#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import dash_table
import math
import numpy as np
import pandas as pd
from pyproj import Transformer
import io
from xhtml2pdf import pisa

# ------------------- Utility Functions -------------------

def dms_to_decimal(degrees, minutes, seconds):
    """Convert DMS to decimal degrees."""
    if None in [degrees, minutes, seconds]:
        return None
    try:
        degrees = float(degrees)
        minutes = float(minutes)
        seconds = float(seconds)
        # Basic range checks
        if not (0 <= degrees <= 360 and 0 <= minutes < 60 and 0 <= seconds < 60):
            return None
        return degrees + minutes / 60 + seconds / 3600
    except ValueError:
        return None

def decimal_to_dms(decimal_degrees):
    """Convert decimal degrees to DMS string."""
    if decimal_degrees is None:
        return "N/A"
    deg = int(decimal_degrees)
    mn = int((decimal_degrees - deg) * 60)
    sc = round(((decimal_degrees - deg) * 60 - mn) * 60, 2)
    return f"{deg:03d}° {mn:02d}' {sc:05.2f}''"

def normalize_angle(angle):
    """Normalize angle to be within [-180, 180]."""
    while angle > 180:
        angle -= 360
    while angle < -180:
        angle += 360
    return angle

def normalize_0_360(angle):
    """Normalize angle to be within [0, 360)."""
    return angle % 360

def dms_input(id_prefix):
    """
    Helper function to build 3 input fields for DMS (Degrees, Minutes, Seconds).
    """
    return dbc.Row([
        dbc.Col(
            dbc.Input(
                id=f"{id_prefix}-degrees",
                placeholder="Degrees",
                type="number",
                min=0,
                max=360,
                step=1,
                required=True
            ),
            width=4
        ),
        dbc.Col(
            dbc.Input(
                id=f"{id_prefix}-minutes",
                placeholder="Minutes",
                type="number",
                min=0,
                max=59,
                step=1,
                required=True
            ),
            width=4
        ),
        dbc.Col(
            dbc.Input(
                id=f"{id_prefix}-seconds",
                placeholder="Seconds",
                type="number",
                min=0,
                max=59.99,
                step=0.01,
                required=True
            ),
            width=4
        ),
    ], className="mb-2")

def calculate_grid_values(e1, n1, e2, n2, hemisphere):
    """
    Calculate both:
      - Grid Bearing (decimal degrees, plus DMS)
      - Grid Convergence (decimal degrees, plus DMS)
    """
    if None in [e1, n1, e2, n2, hemisphere]:
        return None, "N/A", None, "N/A"

    try:
        # Determine UTM zone from e1
        utm_zone = int((e1 // 1000000) % 60) + 1
        south = (hemisphere.lower() == "south")

        # Transformer from UTM->Geographic
        transformer_utm_to_geo = Transformer.from_crs(
            f"+proj=utm +zone={utm_zone} +{'south' if south else ''} +ellps=WGS84",
            "epsg:4326",
            always_xy=True
        )

        lon1, lat1 = transformer_utm_to_geo.transform(e1, n1)
        lon2, lat2 = transformer_utm_to_geo.transform(e2, n2)

        # Grid Bearing
        delta_e = e2 - e1
        delta_n = n2 - n1
        bearing_rad = math.atan2(delta_e, delta_n)
        grid_bearing_dd = (math.degrees(bearing_rad) + 360) % 360

        # Grid Convergence
        central_meridian = utm_zone * 6 - 183
        grid_convergence_rad = math.atan(
            math.tan(math.radians(lon1 - central_meridian)) * math.sin(math.radians(lat1))
        )
        grid_convergence_dd = math.degrees(grid_convergence_rad)

        return (
            grid_bearing_dd, 
            decimal_to_dms(grid_bearing_dd), 
            grid_convergence_dd, 
            decimal_to_dms(grid_convergence_dd)
        )
    except:
        return None, "N/A", None, "N/A"

# ------------------- Dash App & Layout -------------------

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.title = "Heading Sensor Calibration"

# Prepare 19 fixed rows for the Gyrocompass table
default_gyro_data = []
for i in range(1, 21):  # 19 rows
    point_value = "Bow" if i % 2 != 0 else "Stem"
    default_gyro_data.append({
        "utc_time": "",
        "obs_point": point_value,
        "direction_deg": None,
        "direction_min": None,
        "direction_sec": None,
        "obs_distance": None,
        "obs_true_heading": None
    })

app.layout = dbc.Container([
    # --- Header with Title, Device Name Input and Logo ---
    dbc.Row([
        dbc.Col([ 
            html.H1("Heading Sensor Calibration", className="text-center"),
            dbc.Input(id="device-name", placeholder="Enter Device Name", type="text", style={"margin-top": "10px"})
        ], width=8),
        dbc.Col([
            # Example logo; replace as needed
            html.Img(src="https://via.placeholder.com/150", style={"max-width": "100%"})
        ], width=4, className="text-center")
    ], className="my-4"),

    # --- Job Details ---
    dbc.Card([
        dbc.CardHeader(html.H5("Job Details")),
        dbc.CardBody([
            dbc.Row([
                dbc.Col([html.H6("Job Number"), dbc.Input(id="job-number", required=True)], width=3),
                dbc.Col([html.H6("Job Description"), dbc.Input(id="job-description", required=True)], width=3),
                dbc.Col([html.H6("Client"), dbc.Input(id="client", required=True)], width=3),
                dbc.Col([html.H6("Party Chief"), dbc.Input(id="party-chief", required=True)], width=3),
            ]),
            dbc.Row([
                dbc.Col([html.H6("Surveyor"), dbc.Input(id="surveyor", required=True)], width=3),
                dbc.Col([html.H6("Wharf"), dbc.Input(id="wharf", required=True)], width=3),
                dbc.Col([html.H6("Vessel"), dbc.Input(id="vessel", required=True)], width=3),
                dbc.Col([html.H6("Date"), dcc.DatePickerSingle(id="date", placeholder="Select a date", display_format='YYYY-MM-DD')], width=3),
            ]),
            dbc.Row([
                dbc.Col([
                    html.H6("Hemisphere"), 
                    dbc.Select(
                        id="hemisphere", 
                        options=[{"label": "North", "value": "North"}, {"label": "South", "value": "South"}], 
                        placeholder="Select Hemisphere", 
                        required=True
                    )
                ], width=3),
            ])
        ])
    ], className="mb-4"),

    # --- Instrument Station ---
    dbc.Card([
        dbc.CardHeader(html.H5("Instrument Station")),
        dbc.CardBody([
            dbc.Row([
                dbc.Col([html.H6("Instrument Station Name"), dbc.Input(id="instrument-station", required=True)], width=4),
                dbc.Col([html.H6("Easting (m)"), dbc.Input(id="easting", type="number", min=0, required=True)], width=4),
                dbc.Col([html.H6("Northing (m)"), dbc.Input(id="northing", type="number", min=0, required=True)], width=4),
            ]),
            dbc.Row([dbc.Col([html.H6("AHD Height (m)"), dbc.Input(id="AHD-height", type="number", required=True)], width=4)]),
        ])
    ], className="mb-4"),

    # --- Backsight Station ---
    dbc.Card([
        dbc.CardHeader(html.H5("Backsight Station")),
        dbc.CardBody([
            dbc.Row([
                dbc.Col([html.H6("Backsight Station Name"), dbc.Input(id="backsight-station", required=True)], width=4),
                dbc.Col([html.H6("Backsight Easting (m)"), dbc.Input(id="backsight-easting", type="number", min=0, required=True)], width=4),
                dbc.Col([html.H6("Backsight Northing (m)"), dbc.Input(id="backsight-northing", type="number", min=0, required=True)], width=4),
            ]),
            dbc.Row([dbc.Col([html.H6("Backsight AHD Height (m)"), dbc.Input(id="backsight-AHD-height", type="number", required=True)], width=4)]),
        ])
    ], className="mb-4"),

    # --- Calculated Grid Bearing (DMS) ---
    dbc.Card([
        dbc.CardHeader(html.H5("Calculated Grid Bearing (DMS)")),
        dbc.CardBody([html.Div(id="grid-bearing-output", style={"fontWeight": "bold", "fontSize": "16px"})])
    ], className="mb-4"),

    # --- Calculated Grid Convergence (DMS) ---
    dbc.Card([
        dbc.CardHeader(html.H5("Calculated Grid Convergence (DMS)")),
        dbc.CardBody([html.Div(id="grid-convergence-output", style={"fontWeight": "bold", "fontSize": "16px"})])
    ], className="mb-4"),

    # --- Backsight Observation (DMS) ---
    dbc.Card([
        dbc.CardHeader(html.H5("Backsight Observation (DMS)")),
        dbc.CardBody([dms_input("backsight-observation")])
    ], className="mb-4"),

    # --- Gyrocompass Observations: 19 Fixed Rows ---
    dbc.Card([
        dbc.CardHeader(html.H5("Gyrocompass Observations (20 Rows)")),
        dbc.CardBody([
            dash_table.DataTable(
                id='gyro-table',
                columns=[
                    {"name": "UTC Time (hh:mm:ss)",  "id": "utc_time",        "type": "text",    "editable": True},
                    {"name": "Observation Point",    "id": "obs_point",       "type": "text",    "editable": False},
                    {"name": "Direction Deg",        "id": "direction_deg",   "type": "numeric", "editable": True},
                    {"name": "Direction Min",        "id": "direction_min",   "type": "numeric", "editable": True},
                    {"name": "Direction Sec",        "id": "direction_sec",   "type": "numeric", "editable": True},
                    {"name": "Observed Distance (m)","id": "obs_distance",    "type": "numeric", "editable": True},
                    {"name": "Observed True Heading (D.D)", "id": "obs_true_heading", "type": "numeric","editable": True},
                ],
                data=default_gyro_data,
                editable=True,
                row_deletable=False,
                style_cell={
                    'minWidth': '100px', 'width': '150px', 'maxWidth': '200px',
                    'whiteSpace': 'normal',
                    'textAlign': 'left',
                },
                style_header={
                    'backgroundColor': 'rgb(230, 230, 230)',
                    'fontWeight': 'bold'
                },
            ),
        ])
    ], className="mb-4"),

    # --- Submit Button ---
    dbc.Row([
        dbc.Col([
            dbc.Button("Submit", id="submit-button", color="primary", className="mt-4", n_clicks=0)
        ], width=12, className="text-center"),
    ], className="mb-4"),

    # --- Display: Job Details, Results, and Statistical Analysis ---
    dbc.Row([
        dbc.Col([
            html.H4("Job Details"),
            html.Div(id="submitted-job-details")
        ], width=12),
    ], className="mb-4"),

    dbc.Row([
        dbc.Col([
            html.H4("Calculated Results"),
            html.Div(id="results-table")
        ], width=12),
    ], className="mb-4"),

    dbc.Row([
        dbc.Col([
            html.H4("Statistical Analysis"),
            html.Div(id="statistical-analysis", style={"fontWeight": "bold", "fontSize": "18px"})
        ], width=12),
    ], className="mb-4"),

    # --- Export to PDF Button & Download Component ---
    dbc.Row([
        dbc.Col([
            dbc.Button("Export to PDF", id="export-pdf-button", color="secondary", className="mt-4"),
            dcc.Download(id="download-pdf")
        ], width=12, className="text-center")
    ], className="mb-4"),

], fluid=True)

# ------------------- Callbacks -------------------

# 1) Update & display both "Calculated Grid Bearing (DMS)" and "Calculated Grid Convergence (DMS)"
@app.callback(
    [Output("grid-bearing-output", "children"),
     Output("grid-convergence-output", "children")],
    [
        Input("easting", "value"),
        Input("northing", "value"),
        Input("backsight-easting", "value"),
        Input("backsight-northing", "value"),
        Input("hemisphere", "value")
    ]
)
def update_grid_and_convergence(e1, n1, e2, n2, hemisphere):
    gb_dd, gb_dms, gc_dd, gc_dms = calculate_grid_values(e1, n1, e2, n2, hemisphere)
    return gb_dms, gc_dms


@app.callback(
    [
        Output("results-table", "children"),
        Output("statistical-analysis", "children"),
        Output("submitted-job-details", "children"),
    ],
    Input("submit-button", "n_clicks"),
    [
        State("job-number", "value"),
        State("job-description", "value"),
        State("client", "value"),
        State("party-chief", "value"),
        State("surveyor", "value"),
        State("wharf", "value"),
        State("vessel", "value"),
        State("date", "date"),
        State("hemisphere", "value"),
        State("instrument-station", "value"),
        State("easting", "value"),
        State("northing", "value"),
        State("AHD-height", "value"),
        State("backsight-station", "value"),
        State("backsight-easting", "value"),
        State("backsight-northing", "value"),
        State("backsight-AHD-height", "value"),
        # Backsight Observation DMS
        State("backsight-observation-degrees", "value"),
        State("backsight-observation-minutes", "value"),
        State("backsight-observation-seconds", "value"),
        # Table
        State('gyro-table', 'data'),
        State('gyro-table', 'columns'),
        # We also need the displayed Grid Bearing & Convergence
        State("grid-bearing-output", "children"),
        State("grid-convergence-output", "children"),
        # Device name
        State("device-name", "value")
    ]
)
def submit_data(
    n_clicks,
    job_number, job_description, client, party_chief, surveyor,
    wharf, vessel, date, hemisphere,
    instrument_station, e1, n1, ahd1,
    backsight_station, e2, n2, ahd2,
    bs_obs_deg, bs_obs_min, bs_obs_sec,
    gyro_data, gyro_columns,
    grid_bearing_str, grid_convergence_str,
    device_name
):
    # عندما يكون زر الإرسال لم يُضغط بعد
    if n_clicks == 0 or not n_clicks:
        return "", "", ""

    # التحقق من الحقول المطلوبة
    required_inputs = [
        job_number, job_description, client, party_chief, surveyor,
        wharf, vessel, date, hemisphere,
        instrument_station, e1, n1, ahd1,
        backsight_station, e2, n2, ahd2,
        bs_obs_deg, bs_obs_min, bs_obs_sec
    ]
    if any(val is None for val in required_inputs):
        return html.Div("Please fill all required fields."), "", ""

    # تحويل الزاوية DMS للـ Backsight Observation
    backsight_obs_dd = dms_to_decimal(bs_obs_deg, bs_obs_min, bs_obs_sec)
    if backsight_obs_dd is None:
        return html.Div("Backsight Observation (DMS) is invalid."), "", ""

    # حساب Grid Bearing و Convergence
    gb_dd, gb_dms, gc_dd, gc_dms = calculate_grid_values(e1, n1, e2, n2, hemisphere)
    if gb_dd is None or gc_dd is None:
        return html.Div("Could not compute Grid Bearing / Convergence."), "", ""

    # -------------------
    # المرحلة الأولى: نجمع بيانات كل صف ونحسب القيم الأساسية
    # -------------------

    rows_info = []
    total_rows = 0
    successful_rows = 0
    failed_rows = 0

    if not gyro_data:
        gyro_data = []

    for i, row in enumerate(gyro_data):
        time_val = (row.get('utc_time') or "").strip()
        obs_point = (row.get('obs_point') or "").strip()
        direction_deg = row.get('direction_deg')
        direction_min = row.get('direction_min')
        direction_sec = row.get('direction_sec')
        obs_distance  = row.get('obs_distance')
        obs_true_hd   = row.get('obs_true_heading')

        # تحويل قيم النصية إلى float
        def try_float(val):
            try:
                return float(val)
            except:
                return None

        direction_deg = try_float(direction_deg)
        direction_min = try_float(direction_min)
        direction_sec = try_float(direction_sec)
        obs_distance  = try_float(obs_distance)
        obs_true_hd   = try_float(obs_true_hd)

        total_rows += 1

        # إذا كان الصف فارغ تقريبًا فتجاوزه
        if (not time_val and direction_deg is None and 
            direction_min is None and direction_sec is None and 
            obs_distance is None):
            continue

        # محاولة تحويل (deg,min,sec) إلى obs_dir_dd
        obs_dir_dd = dms_to_decimal(direction_deg, direction_min, direction_sec)
        if obs_dir_dd is None or obs_distance is None or obs_distance <= 0:
            # فشل الحساب في هذا الصف
            failed_rows += 1
            # ولكننا ما زلنا نحتفظ بالصف للعرض
            rows_info.append({
                "index": i,
                "time_val": time_val,
                "obs_point": obs_point,
                "obs_dir_dd": None,
                "obs_distance": obs_distance,
                "obs_true_hd": obs_true_hd,
                "plane_bearing_dd": None,
                "calc_e": None,
                "calc_n": None,
                "failed": True
            })
            continue

        # حساب Plane Bearing = gb_dd + (obs_dir_dd - backsight_obs_dd)
        plane_bearing_dd = gb_dd + (obs_dir_dd - backsight_obs_dd)
        plane_bearing_dd = normalize_0_360(plane_bearing_dd)

        # حساب الإحداثيات Calc E/N
        # (هنا معادلة بسيطة كمثال، قد تختلف حسب احتياجاتك)
        plane_bearing_rad = math.radians(plane_bearing_dd)
        calc_e = e1 + obs_distance * math.sin(plane_bearing_rad)
        calc_n = n1 + obs_distance * math.cos(plane_bearing_rad)

        successful_rows += 1

        rows_info.append({
            "index": i,
            "time_val": time_val,
            "obs_point": obs_point,
            "obs_dir_dd": obs_dir_dd,
            "obs_distance": obs_distance,
            "obs_true_hd": obs_true_hd,
            "plane_bearing_dd": plane_bearing_dd,
            "calc_e": calc_e,
            "calc_n": calc_n,
            "failed": False
        })

    # -------------------
    # المرحلة الثانية: حساب Calculated True Heading (و C-O) لصف الـBow
    # -------------------
    # نفترض أن كل Bow في الفهرس (i) يقابله Stem في (i+1) إن كان متاحًا.
    # إن لم يكن متاحًا، تظل القيم N/A.
    co_differences = []

    final_observations = []  # هذه القائمة ستحمل كل صف للعرض النهائي

    # سنبني قاموس (index->data_dict) ليسهل الوصول
    row_dict = { d["index"]: d for d in rows_info }

    # لتسهيل إعادة بناء النتيجة بالترتيب:
    all_indices = sorted(row_dict.keys())

    for i in all_indices:
        row_data = row_dict[i]

        # القيم المشتركة للجميع
        time_val   = row_data["time_val"]
        obs_point  = row_data["obs_point"]
        plane_bear_dd = row_data["plane_bearing_dd"]
        plane_bear_dms = decimal_to_dms(plane_bear_dd) if plane_bear_dd is not None else "N/A"
        dist_str   = f"{row_data['obs_distance']:.2f}" if row_data["obs_distance"] else "N/A"

        calc_e = row_data["calc_e"]
        calc_n = row_data["calc_n"]
        calc_e_str = f"{calc_e:.3f}" if calc_e is not None else "N/A"
        calc_n_str = f"{calc_n:.3f}" if calc_n is not None else "N/A"

        # الحقول الثلاثة الخاصة بالـHeading
        obs_true_heading_str = "N/A"
        calc_true_heading_str = "N/A"
        c_o_str = "N/A"

        # لو هذا الصف فاشل (failed) لا نحسب شيء
        if row_data["failed"]:
            # نضع الصف في النتيجة بدون أي حساب، الكل N/A
            final_observations.append({
                "UTC Time": time_val,
                "Observation Point": obs_point,
                "Plane Bearing (DMS)": "N/A",
                "Obs Dist (m)": "N/A",
                "Calc Easting": "N/A",
                "Calc Northing": "N/A",
                "Obs(O) True Heading (D.D)": "N/A",
                "Calculated True Heading (°)": "N/A",
                "C - O (°)": "N/A"
            })
            continue

        # في حال النجاح، نعرض plane bearing وغيرها
        # إذا obs_point == "bow" نحاول إيجاد الصف التالي إن كان stem
        if obs_point.lower() == "bow":
            # لو أدخل المستخدم Observed True Heading
            if row_data["obs_true_hd"] is not None:
                obs_true_heading_str = f"{row_data['obs_true_hd']:.2f}"

            # نفترض stem هو الصف الذي يليه في الترتيب
            idx_next = i+1
            if idx_next in row_dict:
                stem_data = row_dict[idx_next]
                if stem_data["obs_point"].lower() == "stem" and (not stem_data["failed"]):
                    # نحسب الفرق
                    if row_data["calc_e"] is not None and stem_data["calc_e"] is not None:
                        delta_e = row_data["calc_e"] - stem_data["calc_e"]
                    else:
                        delta_e = None

                    if row_data["calc_n"] is not None and stem_data["calc_n"] is not None:
                        delta_n = row_data["calc_n"] - stem_data["calc_n"]
                    else:
                        delta_n = None

                    if delta_e is not None and delta_n is not None:
                        calc_heading_dd = math.degrees(math.atan2(delta_e, delta_n))
                        # ضبط إلى 0..360
                        if calc_heading_dd < 0:
                            calc_heading_dd += 360
                        # إضافة Grid Convergence
                        if gc_dd is not None:
                            calc_heading_dd = (calc_heading_dd + gc_dd) % 360

                        calc_true_heading_str = f"{calc_heading_dd:.2f}"

                        # لو لدينا Observed True Heading
                        if row_data["obs_true_hd"] is not None:
                            co_val = normalize_angle(calc_heading_dd - row_data["obs_true_hd"])
                            c_o_str = f"{co_val:.2f}"
                            co_differences.append(co_val)

        # بناء صف النتيجة النهائي
        final_observations.append({
            "UTC Time": time_val,
            "Observation Point": obs_point,
            "Plane Bearing (DMS)": plane_bear_dms,
            "Obs Dist (m)": dist_str,
            "Calc Easting": calc_e_str,
            "Calc Northing": calc_n_str,
            "Obs(O) True Heading (D.D)": obs_true_heading_str,
            "Calculated True Heading (°)": calc_true_heading_str,
            "C - O (°)": c_o_str
        })

    # ---------------------------------
    # حساب إحصائيات الـ C - O
    mean_co = None
    std_co = None
    if co_differences:
        mean_co = np.mean(co_differences)
        std_co = np.std(co_differences)

    # ---------------------------------
    # إنشاء جدول تفاصيل الوظيفة
    job_details = {
        "Job Number": job_number,
        "Job Description": job_description,
        "Client": client,
        "Party Chief": party_chief,
        "Surveyor": surveyor,
        "Wharf": wharf,
        "Vessel": vessel,
        "Date": date,
        "Hemisphere": hemisphere,
        "Device Name": device_name if device_name else "N/A",
        "Instrument Station": instrument_station,
        "Instrument Easting": e1,
        "Instrument Northing": n1,
        "Instrument AHD": ahd1,
        "Backsight Station": backsight_station,
        "BS Easting": e2,
        "BS Northing": n2,
        "BS AHD": ahd2,
        "Calculated Grid Bearing (DMS)": gb_dms or "N/A",
        "Calculated Grid Convergence (DMS)": gc_dms or "N/A"
    }
    job_details_table = html.Table([
        html.Tbody([
            html.Tr([html.Th(k), html.Td(str(v))]) for k, v in job_details.items()
        ])
    ], className="table table-bordered table-striped")

    # ---------------------------------
    # إنشاء جدول النتائج (الـ Observations)
    if final_observations:
        columns = [
            "UTC Time",
            "Observation Point",
            "Plane Bearing (DMS)",
            "Obs Dist (m)",
            "Calc Easting",
            "Calc Northing",
            "Obs(O) True Heading (D.D)",
            "Calculated True Heading (°)",
            "C - O (°)"
        ]
        thead = html.Thead([
            html.Tr([html.Th(col) for col in columns])
        ])
        tbody = html.Tbody([
            html.Tr([
                html.Td(row[col]) for col in columns
            ]) for row in final_observations
        ])
        results_table = html.Table([thead, tbody], className="table table-bordered table-striped")
    else:
        results_table = html.Div("No valid Gyrocompass Observations to display.")

    # ---------------------------------
    # التحليل الإحصائي
    stats = [
        f"Total Observations: {total_rows}",
        f"Successful Calculations: {successful_rows}",
        f"Failed Calculations: {failed_rows}"
    ]
    if mean_co is not None and std_co is not None:
        stats.append(f"Mean(C-O): {mean_co:.2f}°")
        stats.append(f"Std Dev(C-O): {std_co:.2f}°")

    stat_div = html.Div([html.Ul([html.Li(s) for s in stats])])

    return results_table, stat_div, job_details_table


if __name__ == "__main__":
    app.run_server(debug=True)

