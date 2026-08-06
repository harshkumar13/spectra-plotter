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

from dash import Dash, dcc, html, callback_context, no_update
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

C_KMS = 299792.458

DARK_PLOT = dict(
    paper_bgcolor="#0D1018",
    plot_bgcolor="#0A0D15",
    font=dict(color="#C9CEDF", family="'SF Mono', ui-monospace, Menlo, monospace", size=11),
)
GRID_COLOR = "#1A2030"
AXIS_TEXT = "#8E96AD"

TRACE_COLORS = [
    "#5CC8FF", "#FFB454", "#4ADE80", "#C4A7FF", "#FF7EB6",
    "#8BD5FF", "#E8C468", "#7EE7A2", "#DFC5FF", "#FF9A8B",
    "#59E0C0", "#B9D5FF", "#FFC98B", "#A5EDB8", "#F0D9FF",
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

        # Header with optical-spectrum signature strip
        html.Div([
            html.Div([
                html.Span("SPECTRA", className="brand-1"),
                html.Span("PLOTTER", className="brand-2"),
            ], className="brand"),
            html.Span("transient spectra console", className="tagline"),
        ], className="app-header"),
        html.Div(className="spectrum-strip"),

        # === TOP ROW: data source + file selection ===
        html.Div([
            html.Div([
                html.Div("Source", className="eyebrow"),
                dcc.RadioItems(
                    id="file-source",
                    options=[
                        {"label": "Upload", "value": "uploads"},
                        {"label": "Local folder", "value": "folder"},
                    ],
                    value="uploads", inline=True, className="radio-seg",
                ),
                html.Div([
                    dcc.Upload(
                        id="upload-data",
                        children=html.Div(["Drop spectra here or ", html.Span("browse", className="upload-link")]),
                        multiple=True, className="upload-area",
                    ),
                ], id="upload-controls"),
                html.Div([
                    html.Div([
                        dcc.Input(id="folder-path", type="text", value=os.path.abspath(UPLOAD_FOLDER),
                                  className="folder-input"),
                        html.Button("Scan", id="scan-folder", n_clicks=0, className="btn-accent"),
                    ], className="folder-row"),
                ], id="folder-controls", style={"display": "none"}),
            ], className="panel source-panel"),

            html.Div([
                html.Div("Spectra", className="eyebrow"),
                dcc.Dropdown(
                    id="file-selector", options=initial_opts, value=[], multi=True,
                    placeholder="Select spectra files...",
                ),
                html.Div([
                    html.Button("Select all", id="select-all-spectra", n_clicks=0),
                    html.Button("Clear", id="clear-spectra-selection", n_clicks=0),
                    html.Button("Reload", id="reload-btn", n_clicks=0),
                    html.Button("Export CSV", id="export-btn", n_clicks=0),
                ], className="btn-row"),
            ], className="panel files-panel"),
        ], className="top-bar"),

        # === MAIN AREA: plot left, line IDs right ===
        html.Div([
            # LEFT: plot + toolbar + redshift + status
            html.Div([
                html.Div([
                    dcc.Loading(
                        dcc.Graph(
                            id="spectra-plot",
                            config={
                                "displayModeBar": True,
                                "displaylogo": False,
                                "modeBarButtonsToAdd": ["drawline", "eraseshape"],
                                "toImageButtonOptions": {"format": "png", "width": 1600, "height": 800, "scale": 2},
                                "scrollZoom": True,
                            },
                            style={"height": "540px"},
                        ),
                        type="circle", color="#FFB454",
                        delay_show=400, delay_hide=100,
                        overlay_style={"visibility": "visible", "opacity": 0.6},
                    ),

                    # Toolbar below plot
                    html.Div([
                        html.Div([
                            html.Span("View", className="tool-label"),
                            html.Button("Zoom full", id="reset-zoom-btn", n_clicks=0),
                            html.Button("Auto zoom", id="auto-zoom-btn", n_clicks=0),
                        ], className="tool-group"),

                        html.Div([
                            html.Span("Process", className="tool-label"),
                            html.Span("bin", className="tool-sub"),
                            dcc.Input(id="bin-size", type="number", value=1, min=1, step=1, className="num-input"),
                            html.Span("smooth", className="tool-sub"),
                            dcc.Dropdown(
                                id="proc-mode",
                                options=[
                                    {"label": "None", "value": "none"},
                                    {"label": "SavGol", "value": "savgol"},
                                    {"label": "Gauss", "value": "gaussian"},
                                ],
                                value="none", clearable=False, className="mini-dropdown",
                            ),
                            dcc.Input(id="smooth-param", type="number", value=11, min=3, step=2,
                                      className="num-input", placeholder="win"),
                            dcc.Checklist(
                                id="cr-clip",
                                options=[{"label": " CR clip", "value": "on"}],
                                value=[], className="checklist inline-check",
                            ),
                        ], className="tool-group"),

                        html.Div([
                            html.Span("Display", className="tool-label"),
                            html.Span("norm", className="tool-sub"),
                            dcc.Dropdown(
                                id="norm-mode",
                                options=[
                                    {"label": "Off", "value": "off"},
                                    {"label": "Each", "value": "each"},
                                    {"label": "Common", "value": "common"},
                                ],
                                value="common", clearable=False, className="mini-dropdown",
                            ),
                            html.Span("offset", className="tool-sub"),
                            dcc.Input(id="trace-offset", type="number", value=0, step=0.1, className="num-input"),
                            dcc.Checklist(
                                id="log-scale",
                                options=[{"label": " Log Y", "value": "log"}],
                                value=[], className="checklist inline-check",
                            ),
                        ], className="tool-group"),
                    ], className="plot-toolbar"),
                ], className="panel plot-panel"),

                # Redshift inputs for selected spectra
                html.Div([
                    html.Div("Redshift", className="eyebrow"),
                    html.Div(id="redshift-inputs", className="z-inputs"),
                    html.Div([
                        dcc.Input(id="z-all-input", type="number", step=0.0001, placeholder="z",
                                  className="num-input z-all-num"),
                        html.Button("Apply to all", id="apply-z-all", n_clicks=0),
                    ], className="z-all-row"),
                ], className="panel z-panel"),

                # Status
                html.Div(id="status-bar", className="status-bar"),
            ], className="main-left"),

            # RIGHT: Line IDs
            html.Div([
                html.Div([
                    html.Span("Line identification", className="eyebrow eyebrow-inline"),
                    html.Span(id="active-lines-badge", className="badge"),
                    html.Button("Clear all", id="clear-all-lines", n_clicks=0, className="btn-small"),
                ], className="line-panel-head"),

                html.Div([
                    html.Div([
                        html.Span("v-shift", className="tool-sub"),
                        dcc.Input(id="velocity-input", type="number", value=0, step=100, className="num-input v-input"),
                        html.Span("km/s (photospheric; skips Telluric & Host)", className="hint"),
                    ], className="vshift-row"),

                    html.Div([
                        html.Span("telluric mask", className="tool-sub"),
                        dcc.Checklist(
                            id="telluric-mask",
                            options=[{"label": " " + k.replace("_", " "), "value": k}
                                     for k in TELLURIC_BANDS],
                            value=DEFAULT_TELLURIC_MASK,
                            inline=True, className="checklist telluric-check",
                        ),
                    ], className="telluric-row"),
                ], className="line-panel-controls"),

                html.Div([
                    _build_line_column(col) for col in
                    [_LINE_GROUPS_COL1, _LINE_GROUPS_COL2, _LINE_GROUPS_COL3, _LINE_GROUPS_COL4]
                ], className="line-grid"),

                html.Div("Click an element name to toggle its whole group.", className="hint line-hint"),
            ], className="panel line-panel"),
        ], className="main-area"),

        # Footer
        html.Div([
            "SpectraPlotter v1.0 · Harsh Kumar (CfA/IAIFI)",
        ], className="app-footer"),
    ], className="app-root")

    _register_callbacks(app)
    return app


def _short_line_label(label: str, element: str, w0: float) -> str:
    """Compact per-line label: 'Hα (6563 Å)' -> 'Hα 6563', 'He I (5876 Å)' -> '5876'."""
    head = label.split("(")[0].strip()
    if head.startswith(element):
        head = head[len(element):].strip()
    return f"{head} {int(round(w0))}".strip()


def _build_line_column(groups):
    items = []
    for element in groups:
        if element not in SPECTRAL_LINES:
            continue
        lines = SPECTRAL_LINES[element]
        color = LINE_COLORS.get(element, "#e6edf3")
        items.append(html.Div([
            html.Div([
                html.Button([
                    html.Span(className="el-dot", style={"background": color}),
                    html.Span(element, style={"color": color}),
                    html.Span(str(len(lines)), className="el-count"),
                ],
                    id={"type": "element-all-btn", "index": element},
                    n_clicks=0, className="el-name-btn",
                    title=f"Toggle all {element} lines",
                ),
            ], className="el-head"),
            dcc.Checklist(
                id={"type": "element-line-selector", "index": element},
                options=[{"label": " " + _short_line_label(lbl, element, w0), "value": lbl}
                         for lbl, w0 in lines.items()],
                value=[],
                className="checklist el-lines",
            ),
        ], className="el-group"))
    return html.Div(items, className="line-col")


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
            return html.Span("Select spectra to set redshifts.", className="hint")
        store = store or {}
        children = []
        for fp in selected:
            z = store.get(fp, 0.0)
            name = os.path.basename(fp)
            short = name if len(name) <= 30 else name[:27] + "..."
            children.append(html.Div([
                html.Span("z=", className="z-label"),
                dcc.Input(
                    id={"type": "redshift-input", "index": fp},
                    type="number", value=float(z), step=0.001,
                    debounce=True, persistence=True, persistence_type="session",
                    className="num-input z-num",
                ),
                html.Span(short, className="z-name", title=name),
            ], className="z-item"))
        return children

    @app.callback(
        Output("redshift-store", "data"),
        Input({"type": "redshift-input", "index": ALL}, "value"),
        State({"type": "redshift-input", "index": ALL}, "id"),
        State("redshift-store", "data"),
    )
    def update_redshift_store(values, ids, store):
        # The displayed inputs are the source of truth: record every current
        # value (file paths contain dots, so parsing prop_id strings is fragile).
        store = store or {}
        for v, _id in zip(values, ids or []):
            if v is not None:
                store[_id["index"]] = float(v)
        return store

    @app.callback(
        Output({"type": "redshift-input", "index": ALL}, "value"),
        Input("apply-z-all", "n_clicks"),
        State("z-all-input", "value"),
        State({"type": "redshift-input", "index": ALL}, "value"),
        prevent_initial_call=True,
    )
    def apply_z_to_all(n_clicks, z, current):
        if not n_clicks or z is None:
            return [no_update] * len(current or [])
        return [float(z)] * len(current or [])

    # Element name button -> toggle all lines in that group; Clear all -> everything off
    @app.callback(
        Output({"type": "element-line-selector", "index": ALL}, "value"),
        Input("clear-all-lines", "n_clicks"),
        Input({"type": "element-all-btn", "index": ALL}, "n_clicks"),
        State({"type": "element-line-selector", "index": ALL}, "value"),
        State({"type": "element-line-selector", "index": ALL}, "options"),
        State({"type": "element-line-selector", "index": ALL}, "id"),
        prevent_initial_call=True,
    )
    def handle_line_toggles(clear_clicks, btn_clicks, all_values, all_options, all_ids):
        n = len(all_ids or [])
        trig = callback_context.triggered
        if not trig:
            return [no_update] * n

        prop_id = trig[0]["prop_id"].rsplit(".", 1)[0]

        if prop_id == "clear-all-lines":
            return [[] for _ in range(n)]

        if not prop_id.startswith("{"):
            return [no_update] * n
        try:
            element = json.loads(prop_id)["index"]
        except Exception:
            return [no_update] * n

        result = []
        for opts, vals, _id in zip(all_options, all_values, all_ids):
            if _id["index"] == element:
                every = [o["value"] for o in opts]
                # toggle: all on if not everything selected, else all off
                result.append([] if set(vals or []) == set(every) else every)
            else:
                result.append(no_update)
        return result

    @app.callback(
        Output("spectra-plot", "figure"),
        Output("status-bar", "children"),
        Output("spectrum-metadata", "data"),
        Output("active-lines-badge", "children"),
        Input("file-selector", "value"),
        Input("redshift-store", "data"),
        Input({"type": "element-line-selector", "index": ALL}, "value"),
        Input("velocity-input", "value"),
        Input("telluric-mask", "value"),
        Input("proc-mode", "value"),
        Input("smooth-param", "value"),
        Input("bin-size", "value"),
        Input("norm-mode", "value"),
        Input("trace-offset", "value"),
        Input("log-scale", "value"),
        Input("cr-clip", "value"),
        Input("reset-zoom-btn", "n_clicks"),
        Input("auto-zoom-btn", "n_clicks"),
        State({"type": "element-line-selector", "index": ALL}, "id"),
        State("spectra-plot", "relayoutData"),
    )
    def update_plot(
        selected_files, redshift_store, selected_lines_nested,
        velocity, telluric_keys, proc_mode, smooth_param, bin_size,
        norm_mode, trace_offset, log_scale_val, cr_clip_val,
        reset_clicks, auto_zoom_clicks, line_ids, relayout,
    ):
        y_scale = "log" if ("log" in (log_scale_val or [])) else "linear"
        redshift_store = redshift_store or {}
        selected_files = selected_files or []

        # View state key: changing the file set or any redshift refits the view,
        # everything else preserves the user's zoom.
        uirev = json.dumps(
            {"f": selected_files, "z": {fp: redshift_store.get(fp, 0.0) for fp in selected_files}},
            sort_keys=True,
        )

        fig = go.Figure()
        fig.update_layout(
            **DARK_PLOT,
            height=540,
            hovermode="x unified",
            hoverlabel=dict(bgcolor="#141927", bordercolor="#2F3850",
                            font=dict(size=11, color="#E8EBF4")),
            legend=dict(
                bgcolor="rgba(13,16,24,0.85)", bordercolor="#232A3D", borderwidth=1,
                font=dict(size=10, color="#C9CEDF"),
                orientation="h", yanchor="bottom", y=1.12, xanchor="left", x=0,
            ),
            margin=dict(l=58, r=12, t=76, b=48),
            uirevision=uirev,
            spikedistance=-1,
        )
        fig.update_xaxes(
            showgrid=True, zeroline=False,
            gridcolor=GRID_COLOR, zerolinecolor=GRID_COLOR,
            title_font=dict(size=11, color=AXIS_TEXT),
            tickfont=dict(size=10, color=AXIS_TEXT),
            range=[2000, 25000],
            showspikes=True, spikemode="across", spikesnap="cursor",
            spikecolor="#3B4560", spikethickness=1, spikedash="solid",
        )
        y_label = "Normalized Flux" if norm_mode != "off" else "Flux"
        fig.update_yaxes(
            title_text=y_label,
            showgrid=True, zeroline=False,
            gridcolor=GRID_COLOR, zerolinecolor=GRID_COLOR,
            title_font=dict(size=11, color=AXIS_TEXT),
            tickfont=dict(size=10, color=AXIS_TEXT),
        )

        if not selected_files:
            fig.update_xaxes(title_text="Wavelength (Å)")
            fig.add_annotation(
                text="Select or upload spectra to begin",
                x=0.5, y=0.5, xref="paper", yref="paper",
                showarrow=False, font=dict(size=14, color="#5A6178"),
            )
            return fig, html.Span("Ready. Select spectra files to plot.", className="hint"), {}, ""

        telluric_keys = telluric_keys or []
        sel_by_el = {i["index"]: set(v or []) for i, v in zip(line_ids or [], selected_lines_nested or [])}
        velocity = float(velocity) if velocity is not None else 0.0
        do_cr = "on" in (cr_clip_val or [])
        bin_size = max(1, int(bin_size or 1))
        smooth_param = int(smooth_param or 11)
        trace_offset = float(trace_offset or 0.0)

        spectra: List[Tuple[str, np.ndarray, np.ndarray]] = []
        applied_z: List[float] = []
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
                applied_z.append(z)
            except Exception as e:
                errors.append(f"{os.path.basename(fp)}: {e}")

        rest_frame = any(abs(z) > 0 for z in applied_z)
        x_title = "Rest Wavelength (Å)" if rest_frame else "Observed Wavelength (Å)"
        fig.update_xaxes(title_text=x_title)

        if not spectra:
            msg = "No readable spectra."
            if errors:
                msg += " | " + "; ".join(errors)
            fig.add_annotation(
                text="Could not read selected files", x=0.5, y=0.5,
                xref="paper", yref="paper", showarrow=False,
                font=dict(size=14, color="#F87171"),
            )
            return fig, html.Span(msg, className="status-error"), metadata, ""

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

        if trace_offset != 0.0:
            spectra = [(n, w, f + i * trace_offset) for i, (n, w, f) in enumerate(spectra)]

        if y_scale == "log":
            spectra = [(n, w, np.where(f > 0, f, np.nan)) for n, w, f in spectra]

        # Plot traces
        all_wmin, all_wmax = 1e10, 0
        for i, (name, w, f) in enumerate(spectra):
            color = TRACE_COLORS[i % len(TRACE_COLORS)]
            fig.add_trace(go.Scattergl(
                x=w, y=f, mode="lines", name=name,
                line=dict(width=1.3, color=color),
                hovertemplate="wl=%{x:.1f}Å flux=%{y:.4g}<extra></extra>",
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
        # Redshift or file-set changes refit the view; other tweaks keep the zoom
        keep_zoom = not (is_reset or is_auto or "redshift-store" in trig or "file-selector" in trig)

        if keep_zoom and relayout:
            x0 = relayout.get("xaxis.range[0]")
            x1 = relayout.get("xaxis.range[1]")
            y0 = relayout.get("yaxis.range[0]")
            y1 = relayout.get("yaxis.range[1]")
            if x0 is not None and x1 is not None:
                fig.update_xaxes(range=[max(0, x0), x1])
            if y0 is not None and y1 is not None:
                fig.update_yaxes(range=[y0, y1])

        if is_auto:
            # Robust y-range: clip outliers so one cosmic ray doesn't set the scale
            allf = np.concatenate([f[np.isfinite(f)] for _, _, f in spectra]) if spectra else np.array([])
            if y_scale == "log":
                allf = allf[allf > 0]
            if allf.size > 10:
                lo, hi = np.nanpercentile(allf, [0.5, 99.7])
                if np.isfinite(lo) and np.isfinite(hi) and hi > lo:
                    pad = (hi - lo) * 0.08
                    if y_scale == "log":
                        fig.update_yaxes(range=[np.log10(max(lo, hi * 1e-6)), np.log10(hi + pad)])
                    else:
                        fig.update_yaxes(range=[lo - pad, hi + pad])

        # Line markers: full-height vlines (log-scale safe, don't affect autoscale)
        z_ref = applied_z[0] if applied_z else 0.0
        n_markers = 0
        for element, lines in SPECTRAL_LINES.items():
            selected = sel_by_el.get(element)
            if not selected:
                continue
            color = LINE_COLORS.get(element, "#e6edf3")
            for label, w0 in lines.items():
                if label not in selected:
                    continue
                x = float(w0)
                note = ""
                if element == "Telluric":
                    # Atmospheric feature: fixed in the observed frame, so map it
                    # onto the (rest-frame) axis with the first spectrum's z
                    x = x / (1.0 + z_ref)
                    if z_ref:
                        note = f" · obs-frame, plotted at z={z_ref:.4f}"
                elif element != "Host Galaxy":
                    x = x * (1.0 + velocity / C_KMS)
                fig.add_shape(
                    type="line", x0=x, x1=x, y0=0, y1=1, yref="y domain",
                    line=dict(color=color, dash="dot", width=1), opacity=0.7,
                )
                fig.add_annotation(
                    x=x, y=1, yref="y domain", yanchor="bottom", yshift=2,
                    text=label.split("(")[0].strip(),
                    showarrow=False, textangle=-55,
                    font=dict(size=9, color=color),
                    hovertext=f"{label}{note}",
                )
                n_markers += 1

        # Status chips (color-matched to traces)
        chips = [html.Span(f"{len(spectra)} spectra", className="chip chip-count")]
        for i, (name, info) in enumerate(metadata.items()):
            color = TRACE_COLORS[i % len(TRACE_COLORS)]
            chips.append(html.Span([
                html.Span(name, className="chip-name"),
                html.Span(f"z={info['z']:.4f}"),
                html.Span(f"SNR~{info['snr']}"),
                html.Span(f"{info['range']}Å"),
                html.Span(f"{info['pts']}pts"),
            ], className="chip", style={"borderLeftColor": color}))
        if errors:
            chips.append(html.Span("ERR: " + "; ".join(errors), className="chip chip-error"))

        badge = f"{n_markers} active" if n_markers else ""
        return fig, chips, metadata, badge

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
    parser.add_argument("--window", action="store_true",
                        help="Open in a local app-style window (no browser tabs/URL bar)")
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    app = create_app(upload_folder=args.folder)
    url = f"http://{args.host}:{args.port}"

    if args.window:
        import subprocess
        import threading
        import time
        import webbrowser

        chrome = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"

        def _open_window():
            time.sleep(1.2)  # let the server come up first
            if os.path.exists(chrome):
                subprocess.Popen([chrome, f"--app={url}", "--window-size=1720,1080"],
                                 stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            else:
                webbrowser.open(url)

        threading.Thread(target=_open_window, daemon=True).start()

    print(f"\n  SpectraPlotter running at {url}\n")
    app.run(host=args.host, port=args.port, debug=args.debug)


if __name__ == "__main__":
    main()
