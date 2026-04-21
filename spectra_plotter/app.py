"""SpectraPlotter - Interactive astronomical spectra visualization.

A Dash web application for loading, processing, and visualizing astronomical
spectra with interactive line identification, inspired by TNS/WISeREP tools.
"""

import base64
import json
import os
from typing import List, Tuple

import numpy as np
import plotly.graph_objs as go

from dash import Dash, dcc, html, callback_context
from dash.dependencies import Input, Output, State, ALL

from .io import list_spectra_files, read_spectrum
from .processing import (
    clean_spectrum, apply_telluric_mask, cosmic_ray_clip,
    apply_savgol, apply_binning, apply_gaussian_smooth, compute_snr,
)
from .line_catalog import (
    SPECTRAL_LINES, LINE_COLORS, TELLURIC_BANDS, DEFAULT_TELLURIC_MASK,
)

UPLOAD_FOLDER = os.path.join(os.path.dirname(__file__), "uploads")
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

DARK_PLOT = dict(
    paper_bgcolor="#0d1117",
    plot_bgcolor="#161b22",
    font=dict(color="#c9d1d9", family="system-ui, -apple-system, sans-serif", size=11),
)

TRACE_COLORS = [
    "#58a6ff", "#f0883e", "#3fb950", "#bc8cff", "#f778ba",
    "#79c0ff", "#d29922", "#56d364", "#d2a8ff", "#ff7b72",
    "#39d353", "#a5d6ff", "#ffa657", "#7ee787", "#e2c5ff",
]

# WISeREP-style: group line IDs into 4 columns for compact display
_LINE_GROUPS_COL1 = [
    "H (Balmer)", "H (Paschen)", "H (Brackett)",
    "He I", "He II",
    "C I", "C II",
    "N II",
]
_LINE_GROUPS_COL2 = [
    "O I", "O II", "[O I]", "[O II]", "[O III]",
    "Si I", "Si II",
    "S I", "S II", "[S III]",
]
_LINE_GROUPS_COL3 = [
    "Na I", "Mg I", "Mg I]", "Mg II",
    "Ca II", "[Ca II]",
    "Fe II", "[Fe II]", "Fe III",
    "Ni II", "[Ni II]",
]
_LINE_GROUPS_COL4 = [
    "[Ar III]", "Sc II", "Ba II", "Ti II", "Cr II",
    "Telluric",
    "Host Galaxy",
    "SLSN-I Pre-peak", "SLSN-I Peak", "SLSN-I Late", "SLSN-I Nebular",
]


def create_app(upload_folder: str = None) -> Dash:
    if upload_folder:
        global UPLOAD_FOLDER
        UPLOAD_FOLDER = upload_folder
        os.makedirs(UPLOAD_FOLDER, exist_ok=True)

    assets = os.path.join(os.path.dirname(__file__), "assets")
    app = Dash(
        __name__,
        assets_folder=assets,
        title="SpectraPlotter",
        update_title="Loading...",
    )

    initial_opts = [
        {"label": os.path.basename(p), "value": p}
        for p in list_spectra_files(UPLOAD_FOLDER)
    ]

    app.layout = html.Div([
        dcc.Store(id="redshift-store", data={}),
        dcc.Store(id="spectrum-metadata", data={}),
        dcc.Download(id="download-data"),

        # Title bar
        html.Div([
            html.Span("SpectraPlotter", style={
                "fontSize": "18px", "fontWeight": "600", "color": "#e6edf3",
            }),
            html.Span(" | Interactive Spectra Viewer", style={
                "fontSize": "12px", "color": "#8b949e", "marginLeft": "8px",
            }),
        ], style={
            "padding": "8px 16px", "borderBottom": "1px solid #30363d",
            "background": "#161b22",
        }),

        # === TOP ROW: spectra selection + data source ===
        html.Div([
            # File source controls
            html.Div([
                dcc.RadioItems(
                    id="file-source",
                    options=[
                        {"label": " Upload", "value": "uploads"},
                        {"label": " Local folder", "value": "folder"},
                    ],
                    value="uploads", inline=True, className="radio-items",
                    style={"fontSize": "12px", "marginBottom": "4px"},
                ),
                html.Div([
                    dcc.Upload(
                        id="upload-data",
                        children=html.Div(["Drop files or ", html.Span("browse", style={"color": "#58a6ff"})]),
                        multiple=True, className="upload-area",
                        style={"height": "36px", "lineHeight": "36px", "textAlign": "center", "fontSize": "12px"},
                    ),
                ], id="upload-controls"),
                html.Div([
                    html.Div([
                        dcc.Input(id="folder-path", type="text", value=os.path.abspath(UPLOAD_FOLDER),
                                  style={"flex": "1", "fontSize": "12px"}),
                        html.Button("Scan", id="scan-folder", n_clicks=0, style={"fontSize": "11px", "padding": "4px 10px"}),
                    ], style={"display": "flex", "gap": "4px"}),
                ], id="folder-controls", style={"display": "none"}),
            ], style={"flex": "0 0 240px"}),

            # File selector
            html.Div([
                dcc.Dropdown(
                    id="file-selector", options=initial_opts, value=[], multi=True,
                    placeholder="Select spectra files...",
                    style={"fontSize": "12px"},
                ),
            ], style={"flex": "1 1 auto", "minWidth": "200px"}),

            # Actions
            html.Div([
                html.Button("Select all", id="select-all-spectra", n_clicks=0,
                            style={"fontSize": "11px", "padding": "4px 8px", "marginRight": "4px"}),
                html.Button("Clear", id="clear-spectra-selection", n_clicks=0,
                            style={"fontSize": "11px", "padding": "4px 8px", "marginRight": "4px"}),
                html.Button("Export CSV", id="export-btn", n_clicks=0,
                            style={"fontSize": "11px", "padding": "4px 8px", "marginRight": "4px"}),
                html.Button("Reload", id="reload-btn", n_clicks=0,
                            style={"fontSize": "11px", "padding": "4px 8px"}),
            ], style={"display": "flex", "alignItems": "center", "flexWrap": "wrap"}),
        ], style={
            "display": "flex", "gap": "10px", "padding": "8px 12px",
            "borderBottom": "1px solid #30363d", "alignItems": "center",
            "flexWrap": "wrap", "background": "#0d1117",
        }),

        # === MAIN AREA: plot left, line IDs right ===
        html.Div([
            # LEFT: Plot + controls below
            html.Div([
                # Plot
                dcc.Graph(
                    id="spectra-plot",
                    config={
                        "displayModeBar": True,
                        "modeBarButtonsToAdd": ["drawline", "eraseshape"],
                        "toImageButtonOptions": {"format": "png", "width": 1600, "height": 800, "scale": 2},
                        "scrollZoom": True,
                    },
                    style={"height": "520px"},
                ),

                # Toolbar below plot (WISeREP-style)
                html.Div([
                    html.Button("Zoom Full", id="reset-zoom-btn", n_clicks=0,
                                style={"fontSize": "11px", "padding": "3px 10px"}),
                    html.Button("Auto Zoom", id="auto-zoom-btn", n_clicks=0,
                                style={"fontSize": "11px", "padding": "3px 10px", "marginLeft": "4px"}),

                    html.Span("Binning:", style={"marginLeft": "12px", "fontSize": "11px", "color": "#8b949e"}),
                    dcc.Input(id="bin-size", type="number", value=1, min=1, step=1,
                              style={"width": "50px", "fontSize": "11px", "marginLeft": "4px"}),

                    html.Span("Smooth:", style={"marginLeft": "12px", "fontSize": "11px", "color": "#8b949e"}),
                    dcc.Dropdown(
                        id="proc-mode",
                        options=[
                            {"label": "None", "value": "none"},
                            {"label": "SavGol", "value": "savgol"},
                            {"label": "Gauss", "value": "gaussian"},
                        ],
                        value="none", clearable=False,
                        style={"width": "90px", "fontSize": "11px", "display": "inline-block", "verticalAlign": "middle"},
                    ),
                    dcc.Input(id="smooth-param", type="number", value=11, min=3, step=2,
                              style={"width": "50px", "fontSize": "11px", "marginLeft": "4px"},
                              placeholder="win"),

                    html.Span("|", style={"margin": "0 8px", "color": "#30363d"}),

                    html.Span("Norm:", style={"fontSize": "11px", "color": "#8b949e"}),
                    dcc.Dropdown(
                        id="norm-mode",
                        options=[
                            {"label": "Off", "value": "off"},
                            {"label": "Each", "value": "each"},
                            {"label": "Common", "value": "common"},
                        ],
                        value="common", clearable=False,
                        style={"width": "85px", "fontSize": "11px", "display": "inline-block", "verticalAlign": "middle"},
                    ),

                    dcc.Checklist(
                        id="cr-clip",
                        options=[{"label": " CR clip", "value": "on"}],
                        value=[], className="checklist",
                        style={"display": "inline-block", "marginLeft": "8px", "fontSize": "11px"},
                    ),

                    dcc.Checklist(
                        id="log-scale",
                        options=[{"label": " Log Y", "value": "log"}],
                        value=[], className="checklist",
                        style={"display": "inline-block", "marginLeft": "4px", "fontSize": "11px"},
                    ),
                ], style={
                    "display": "flex", "alignItems": "center", "flexWrap": "wrap",
                    "padding": "6px 8px", "gap": "2px",
                    "background": "#161b22", "borderTop": "1px solid #30363d",
                    "borderRadius": "0 0 8px 8px",
                }),

                # Redshift inputs for selected spectra
                html.Div(id="redshift-inputs", style={"padding": "4px 8px"}),

                # Status
                html.Div(id="status-bar", style={
                    "padding": "4px 8px", "fontSize": "11px", "color": "#8b949e",
                    "fontFamily": "monospace", "borderTop": "1px solid #21262d",
                }),
            ], style={"flex": "1 1 auto", "minWidth": "500px"}),

            # RIGHT: Line IDs (WISeREP-style compact grid)
            html.Div([
                html.Div([
                    html.Span("Line Identification", style={"fontSize": "13px", "fontWeight": "600", "color": "#e6edf3"}),
                    html.Button("Clear all", id="clear-all-lines", n_clicks=0,
                                style={"marginLeft": "auto", "fontSize": "10px", "padding": "2px 6px"}),
                ], style={"display": "flex", "alignItems": "center", "marginBottom": "6px", "paddingBottom": "4px", "borderBottom": "1px solid #30363d"}),

                # Velocity shift
                html.Div([
                    html.Span("v-shift:", style={"fontSize": "11px", "color": "#8b949e"}),
                    dcc.Input(id="velocity-input", type="number", value=0, step=100,
                              style={"width": "70px", "fontSize": "11px", "marginLeft": "4px"}),
                    html.Span("km/s", style={"fontSize": "10px", "color": "#6e7681", "marginLeft": "2px"}),
                ], style={"marginBottom": "6px", "display": "flex", "alignItems": "center"}),

                # Telluric
                html.Div([
                    html.Span("Telluric:", style={"fontSize": "11px", "color": "#8b949e", "marginRight": "6px"}),
                    dcc.Checklist(
                        id="telluric-mask",
                        options=[{"label": " " + v[0].split("(")[0].strip(), "value": k} for k, v in TELLURIC_BANDS.items()],
                        value=DEFAULT_TELLURIC_MASK,
                        inline=True,
                        className="checklist",
                        style={"fontSize": "10px", "display": "inline"},
                    ),
                ], style={"marginBottom": "6px", "display": "flex", "alignItems": "flex-start", "flexWrap": "wrap"}),

                html.Hr(style={"border": "none", "borderTop": "1px solid #30363d", "margin": "4px 0"}),

                # 4-column line ID grid
                html.Div([
                    _build_line_column(col) for col in [_LINE_GROUPS_COL1, _LINE_GROUPS_COL2, _LINE_GROUPS_COL3, _LINE_GROUPS_COL4]
                ], style={
                    "display": "grid",
                    "gridTemplateColumns": "1fr 1fr 1fr 1fr",
                    "gap": "4px",
                    "maxHeight": "480px",
                    "overflowY": "auto",
                }),
            ], style={
                "flex": "0 0 580px", "padding": "8px 10px",
                "background": "#0d1117", "borderLeft": "1px solid #30363d",
                "overflowY": "auto", "maxHeight": "650px",
            }),
        ], style={
            "display": "flex", "gap": "0",
            "padding": "0", "flexWrap": "wrap",
        }),

        # Footer
        html.Div([
            "SpectraPlotter v1.0 | Harsh Kumar (CfA/IAIFI)",
        ], style={
            "textAlign": "center", "padding": "6px", "color": "#6e7681",
            "fontSize": "11px", "borderTop": "1px solid #30363d",
        }),

        # Hidden elements for savgol polyorder (always 3 default)
        dcc.Store(id="savgol-polyorder-store", data=3),
    ], style={
        "backgroundColor": "#0d1117", "minHeight": "100vh",
        "color": "#e6edf3", "fontFamily": "system-ui, -apple-system, sans-serif",
    })

    _register_callbacks(app)
    return app


def _build_line_column(groups):
    items = []
    for element in groups:
        if element not in SPECTRAL_LINES:
            continue
        lines = SPECTRAL_LINES[element]
        color = LINE_COLORS.get(element, "#e6edf3")
        items.append(html.Div([
            html.Div([
                dcc.Checklist(
                    id={"type": "element-toggle", "index": element},
                    options=[{"label": "", "value": "on"}],
                    value=[],
                    style={"display": "inline-block", "marginRight": "2px"},
                    className="checklist",
                ),
                html.Span(element, style={
                    "color": color, "fontSize": "11px", "fontWeight": "600", "cursor": "pointer",
                }),
            ], style={"display": "flex", "alignItems": "center", "marginBottom": "1px"}),
            dcc.Checklist(
                id={"type": "element-line-selector", "index": element},
                options=[{"label": " " + lbl.split("(")[0].strip(), "value": lbl} for lbl in lines],
                value=[],
                className="checklist",
                style={"paddingLeft": "16px", "fontSize": "10px", "lineHeight": "1.4"},
            ),
        ], style={"marginBottom": "4px"}))
    return html.Div(items)


def _register_callbacks(app):

    @app.callback(
        Output("folder-controls", "style"),
        Output("upload-controls", "style"),
        Input("file-source", "value"),
    )
    def toggle_source(source):
        if source == "folder":
            return {"display": "block"}, {"display": "none"}
        return {"display": "none"}, {"display": "block"}

    @app.callback(
        Output("file-selector", "options"),
        Output("file-selector", "value"),
        Input("upload-data", "contents"),
        State("upload-data", "filename"),
        Input("scan-folder", "n_clicks"),
        Input("reload-btn", "n_clicks"),
        State("folder-path", "value"),
        Input("file-source", "value"),
        State("file-selector", "value"),
        Input("select-all-spectra", "n_clicks"),
        Input("clear-spectra-selection", "n_clicks"),
    )
    def update_file_options(contents, filenames, scan_clicks, reload_clicks,
                            folder_path, source, current,
                            sel_all_clicks, clear_clicks):
        if contents and filenames:
            for content, fname in zip(contents, filenames):
                try:
                    data = content.split(",", 1)[1]
                    path = os.path.join(UPLOAD_FOLDER, os.path.basename(fname))
                    with open(path, "wb") as f:
                        f.write(base64.b64decode(data))
                except Exception:
                    pass

        root = UPLOAD_FOLDER if source == "uploads" else (folder_path or UPLOAD_FOLDER)
        files = list_spectra_files(root)
        options = [{"label": os.path.basename(p), "value": p} for p in files]
        current = current or []

        trig = callback_context.triggered[0]["prop_id"] if callback_context.triggered else ""

        if "clear-spectra-selection" in trig:
            return options, []
        if "select-all-spectra" in trig:
            return options, files

        selection = [p for p in current if p in files]

        if "upload-data" in trig and source == "uploads" and filenames:
            for fn in filenames:
                p = os.path.join(UPLOAD_FOLDER, os.path.basename(fn))
                if p in files and p not in selection:
                    selection.append(p)

        return options, selection

    @app.callback(
        Output("redshift-inputs", "children"),
        Input("file-selector", "value"),
        State("redshift-store", "data"),
    )
    def generate_redshift_inputs(selected, store):
        if not selected:
            return []
        store = store or {}
        children = []
        for fp in selected:
            z = store.get(fp, 0.0)
            name = os.path.basename(fp)
            if len(name) > 30:
                name = name[:27] + "..."
            children.append(html.Div([
                html.Span(f"z=", style={"fontSize": "11px", "color": "#8b949e", "marginRight": "2px"}),
                dcc.Input(
                    id={"type": "redshift-input", "index": fp},
                    type="number", value=float(z), step=0.001,
                    debounce=False, persistence=True, persistence_type="session",
                    style={"width": "80px", "fontSize": "11px"},
                ),
                html.Span(name, style={"fontSize": "11px", "color": "#8b949e", "marginLeft": "6px"}),
            ], style={"display": "inline-flex", "alignItems": "center", "marginRight": "12px", "marginBottom": "2px"}))
        return children

    @app.callback(
        Output("redshift-store", "data"),
        Input({"type": "redshift-input", "index": ALL}, "value"),
        State({"type": "redshift-input", "index": ALL}, "id"),
        State("redshift-store", "data"),
    )
    def update_redshift_store(values, ids, store):
        store = store or {}
        if not ids or not callback_context.triggered:
            return store
        trig = callback_context.triggered[0]["prop_id"].split(".", 1)[0]
        if not trig.startswith("{"):
            return store
        try:
            trig_id = json.loads(trig)
        except Exception:
            return store
        for v, _id in zip(values, ids):
            if _id == trig_id and v is not None:
                store[_id["index"]] = float(v)
                break
        return store

    # Element toggle -> select/deselect all lines in that group
    @app.callback(
        Output({"type": "element-line-selector", "index": ALL}, "value"),
        Input("clear-all-lines", "n_clicks"),
        Input({"type": "element-toggle", "index": ALL}, "value"),
        State({"type": "element-line-selector", "index": ALL}, "options"),
        State({"type": "element-line-selector", "index": ALL}, "value"),
        State({"type": "element-line-selector", "index": ALL}, "id"),
        State({"type": "element-toggle", "index": ALL}, "id"),
    )
    def handle_line_toggles(clear_clicks, toggle_vals, all_options, all_values, all_ids, toggle_ids):
        trig = callback_context.triggered
        if not trig:
            return all_values or [[] for _ in SPECTRAL_LINES]

        prop_id = trig[0]["prop_id"]

        if prop_id == "clear-all-lines.n_clicks":
            return [[] for _ in SPECTRAL_LINES]

        if prop_id.startswith("{"):
            try:
                trig_id = json.loads(prop_id.split(".")[0])
                element = trig_id["index"]
                is_on = "on" in (trig[0].get("value") or trig[0].get("value", []) or [])
            except Exception:
                return all_values or [[] for _ in SPECTRAL_LINES]

            result = []
            for opts, vals, _id in zip(all_options, all_values, all_ids):
                if _id["index"] == element:
                    if is_on:
                        result.append([o["value"] for o in opts])
                    else:
                        result.append([])
                else:
                    result.append(vals or [])
            return result

        return all_values or [[] for _ in SPECTRAL_LINES]

    @app.callback(
        Output("spectra-plot", "figure"),
        Output("status-bar", "children"),
        Output("spectrum-metadata", "data"),
        Input("file-selector", "value"),
        Input("redshift-store", "data"),
        Input({"type": "element-line-selector", "index": ALL}, "value"),
        Input("velocity-input", "value"),
        Input("telluric-mask", "value"),
        Input("proc-mode", "value"),
        Input("smooth-param", "value"),
        Input("bin-size", "value"),
        Input("norm-mode", "value"),
        Input("log-scale", "value"),
        Input("cr-clip", "value"),
        Input("reset-zoom-btn", "n_clicks"),
        Input("auto-zoom-btn", "n_clicks"),
        State("spectra-plot", "relayoutData"),
    )
    def update_plot(
        selected_files, redshift_store, selected_lines_nested,
        velocity, telluric_keys, proc_mode, smooth_param, bin_size,
        norm_mode, log_scale_val, cr_clip_val, reset_clicks, auto_zoom_clicks,
        relayout,
    ):
        y_scale = "log" if ("log" in (log_scale_val or [])) else "linear"

        fig = go.Figure()
        fig.update_layout(
            **DARK_PLOT,
            height=520,
            hovermode="x unified",
            legend=dict(
                bgcolor="rgba(22,27,34,0.9)", bordercolor="#30363d", borderwidth=1,
                font=dict(size=10, color="#e6edf3"),
                orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0,
            ),
            margin=dict(l=55, r=10, t=30, b=45),
            uirevision="main",
        )
        fig.update_xaxes(
            title_text="Rest Wavelength (Angstrom)",
            showgrid=True, zeroline=False,
            gridcolor="#21262d", zerolinecolor="#30363d",
            title_font=dict(size=11, color="#8b949e"),
            tickfont=dict(size=10, color="#8b949e"),
            rangemode="tozero",
            range=[2000, 25000],
        )
        y_label = "Normalized Flux" if norm_mode != "off" else "Flux"
        fig.update_yaxes(
            title_text=y_label,
            showgrid=True, zeroline=False,
            gridcolor="#21262d", zerolinecolor="#30363d",
            title_font=dict(size=11, color="#8b949e"),
            tickfont=dict(size=10, color="#8b949e"),
        )

        if not selected_files:
            fig.add_annotation(
                text="Select or upload spectra to begin",
                x=0.5, y=0.5, xref="paper", yref="paper",
                showarrow=False, font=dict(size=14, color="#8b949e"),
            )
            return fig, "Ready. Select spectra files to plot.", {}

        redshift_store = redshift_store or {}
        telluric_keys = telluric_keys or []
        selected_lines = [ln for sub in (selected_lines_nested or []) for ln in (sub or [])]
        velocity = float(velocity) if velocity is not None else 0.0
        do_cr = "on" in (cr_clip_val or [])
        bin_size = max(1, int(bin_size or 1))
        smooth_param = int(smooth_param or 11)

        spectra: List[Tuple[str, np.ndarray, np.ndarray]] = []
        metadata = {}
        errors = []

        for fp in selected_files:
            try:
                wave_obs, flux = read_spectrum(fp)
                wave_obs, flux = clean_spectrum(wave_obs, flux)
                if wave_obs.size < 2:
                    errors.append(f"{os.path.basename(fp)}: too few points")
                    continue
                flux = apply_telluric_mask(wave_obs, flux, telluric_keys)
                if do_cr:
                    flux = cosmic_ray_clip(flux, kernel=9, sigma=8.0)
                z = float(redshift_store.get(fp, 0.0) or 0.0)
                wave = wave_obs / (1.0 + z)

                if bin_size > 1:
                    wave, flux = apply_binning(wave, flux, bin_size)

                if proc_mode == "savgol":
                    flux = apply_savgol(flux, window=smooth_param, polyorder=3)
                elif proc_mode == "gaussian":
                    flux = apply_gaussian_smooth(flux, sigma=max(1, smooth_param // 2))

                name = os.path.basename(fp)
                snr = compute_snr(flux)
                wmin = float(np.nanmin(wave)) if wave.size > 0 and np.isfinite(np.nanmin(wave)) else 0
                wmax = float(np.nanmax(wave)) if wave.size > 0 and np.isfinite(np.nanmax(wave)) else 0
                metadata[name] = {"z": z, "snr": round(snr, 1), "range": f"{wmin:.0f}-{wmax:.0f}", "pts": int(wave.size)}
                spectra.append((name, wave, flux))
            except Exception as e:
                errors.append(f"{os.path.basename(fp)}: {e}")

        if not spectra:
            msg = "No readable spectra."
            if errors:
                msg += " | " + "; ".join(errors)
            fig.add_annotation(
                text="Could not read selected files", x=0.5, y=0.5,
                xref="paper", yref="paper", showarrow=False,
                font=dict(size=14, color="#f85149"),
            )
            return fig, msg, metadata

        # Normalize
        if norm_mode == "common" and len(spectra) > 1:
            mins = [np.nanmin(w) for _, w, _ in spectra if np.isfinite(w).any()]
            maxs = [np.nanmax(w) for _, w, _ in spectra if np.isfinite(w).any()]
            if mins and maxs:
                lo, hi = float(max(mins)), float(min(maxs))
                if np.isfinite(lo) and np.isfinite(hi) and hi > lo:
                    common_wl = np.linspace(lo, hi, 1500)
                    norms = []
                    for _, w, f in spectra:
                        m = np.isfinite(w) & np.isfinite(f)
                        interp_f = np.interp(common_wl, w[m], f[m]) if m.sum() >= 2 else np.ones_like(common_wl)
                        nrm = np.nanmax(interp_f)
                        norms.append(float(nrm) if np.isfinite(nrm) and nrm != 0 else 1.0)
                    spectra = [(n, w, f / nrm) for (n, w, f), nrm in zip(spectra, norms)]
                else:
                    norm_mode = "each"

        if norm_mode == "each":
            spectra = [
                (n, w, f / (float(np.nanmax(f)) if np.isfinite(np.nanmax(f)) and np.nanmax(f) != 0 else 1.0))
                for n, w, f in spectra
            ]

        if y_scale == "log":
            spectra = [(n, w, np.where(f > 0, f, np.nan)) for n, w, f in spectra]

        # Plot traces
        all_wmin, all_wmax = 1e10, 0
        for i, (name, w, f) in enumerate(spectra):
            color = TRACE_COLORS[i % len(TRACE_COLORS)]
            fig.add_trace(go.Scattergl(
                x=w, y=f, mode="lines", name=name,
                line=dict(width=1.3, color=color),
                hovertemplate="wl=%{x:.1f}A flux=%{y:.4g}<extra></extra>",
                connectgaps=False,
            ))
            valid_w = w[np.isfinite(w) & np.isfinite(f)]
            if valid_w.size > 0:
                all_wmin = min(all_wmin, float(np.nanmin(valid_w)))
                all_wmax = max(all_wmax, float(np.nanmax(valid_w)))

        if y_scale == "log":
            fig.update_yaxes(type="log")

        # Set sensible x-axis range (positive wavelengths only)
        if all_wmin < all_wmax:
            pad = (all_wmax - all_wmin) * 0.02
            fig.update_xaxes(range=[max(0, all_wmin - pad), all_wmax + pad])

        trig = callback_context.triggered[0]["prop_id"] if callback_context.triggered else ""
        is_reset = "reset-zoom-btn" in trig
        is_auto = "auto-zoom-btn" in trig

        if not is_reset and not is_auto and relayout:
            x0 = relayout.get("xaxis.range[0]")
            x1 = relayout.get("xaxis.range[1]")
            y0 = relayout.get("yaxis.range[0]")
            y1 = relayout.get("yaxis.range[1]")
            if x0 is not None and x1 is not None:
                fig.update_xaxes(range=[max(0, x0), x1])
            if y0 is not None and y1 is not None:
                fig.update_yaxes(range=[y0, y1])

        # Line markers
        y_max = 1.0
        try:
            y_max = max(np.nanmax(t.y) for t in fig.data if hasattr(t, "mode") and t.mode == "lines" and t.y is not None)
            if not np.isfinite(y_max) or y_max <= 0:
                y_max = 1.0
        except (ValueError, TypeError):
            y_max = 1.0

        for element, lines in SPECTRAL_LINES.items():
            color = LINE_COLORS.get(element, "#e6edf3")
            for label, w0 in lines.items():
                if label not in selected_lines:
                    continue
                shifted = float(w0) * (1.0 + velocity / 3e5)
                fig.add_trace(go.Scatter(
                    x=[shifted, shifted], y=[0, y_max],
                    mode="lines",
                    line=dict(color=color, dash="dot", width=0.8),
                    name=label, showlegend=False, hoverinfo="skip",
                ))
                fig.add_annotation(
                    x=shifted, y=y_max, text=label.split("(")[0].strip(),
                    showarrow=False, yshift=8,
                    font=dict(size=8, color=color),
                    textangle=-60,
                )

        # Status
        parts = [f"{len(spectra)} spectra"]
        for name, info in metadata.items():
            parts.append(f"{name} z={info['z']:.4f} SNR~{info['snr']} ({info['range']}A, {info['pts']}pts)")
        if errors:
            parts.append("ERR: " + "; ".join(errors))
        status = " | ".join(parts)

        return fig, status, metadata

    @app.callback(
        Output("download-data", "data"),
        Input("export-btn", "n_clicks"),
        State("file-selector", "value"),
        State("redshift-store", "data"),
        prevent_initial_call=True,
    )
    def export_csv(n_clicks, selected, store):
        if not selected or not n_clicks:
            return None
        store = store or {}
        import io as sio
        import csv
        output = sio.StringIO()
        writer = csv.writer(output)
        writer.writerow(["filename", "wavelength_rest_A", "flux"])
        for fp in selected:
            try:
                w, f = read_spectrum(fp)
                w, f = clean_spectrum(w, f)
                z = float(store.get(fp, 0.0) or 0.0)
                w_rest = w / (1.0 + z)
                name = os.path.basename(fp)
                for wi, fi in zip(w_rest, f):
                    if np.isfinite(wi) and np.isfinite(fi):
                        writer.writerow([name, f"{wi:.4f}", f"{fi:.6g}"])
            except Exception:
                pass
        return dict(content=output.getvalue(), filename="spectra_export.csv")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="SpectraPlotter - Astronomical spectra viewer")
    parser.add_argument("--port", type=int, default=8050)
    parser.add_argument("--host", type=str, default="127.0.0.1")
    parser.add_argument("--folder", type=str, default=None, help="Default spectra folder")
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    app = create_app(upload_folder=args.folder)
    print(f"\n  SpectraPlotter running at http://{args.host}:{args.port}\n")
    app.run(host=args.host, port=args.port, debug=args.debug)


if __name__ == "__main__":
    main()
