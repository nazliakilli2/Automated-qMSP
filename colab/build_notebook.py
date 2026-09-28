"""Builds colab/autoqmsp_colab.ipynb. Run: python3 colab/build_notebook.py"""
import json, os

def md(s):
    return {"cell_type": "markdown", "metadata": {}, "source": s.strip("\n").splitlines(True)}

def form(s):
    # cellView "form" hides the code; users see only the title and the fields.
    return {"cell_type": "code", "metadata": {"cellView": "form"}, "execution_count": None,
            "outputs": [], "source": s.strip("\n").splitlines(True)}

INSTALL = r'''
#@title ① Install (run once, takes 1–3 minutes) { display-mode: "form" }
#@markdown Click ▶ on the left. Wait until you see **Ready**.
import json, math, os, re, shutil, subprocess
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from IPython.display import display, HTML

def _run(cmd):
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if p.returncode != 0:
        print(p.stdout[-3000:], p.stderr[-3000:])
        raise SystemExit("Installation failed - see the messages above.")
    return p.stdout

if shutil.which("Rscript") is None:
    print("Installing R ...")
    _run("apt-get -qq update && apt-get -qq install -y r-base > /dev/null")
print("Installing autoqmsp ...")
_run("""Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org"));
  for (p in c("remotes", "writexl", "xml2", "ggplot2"))
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, quiet = TRUE);
  remotes::install_github("nazliakilli2/Automated-qMSP", upgrade = "never", quiet = TRUE);
  library(autoqmsp)'""")

# ---- shared by the next steps ---------------------------------------------
EDS_DIR = "/content/eds"
OUT_DIR = "/content/results"
DATA_DIR = os.path.join(OUT_DIR, "data")
REFERENCE_RE = r"ACTB|B.?ACTIN|BETA.?ACTIN"
NTC_RE = r"NTC|dH2O|dH20|water|blank|^NK"
POSITIVE_RE = r"H460|A549|HT29|positive|^PC\b"
THRESHOLD = 10000.0   # delta Rn where the Ct is read (lab standard)
CT_CUTOFF = 40.0      # methylated when Ct <= this (lab standard)
COLORS = {"detected": "#c0392b", "not detected": "#a7b1bc",
          "positive": "#1f6fb4", "ntc": "#222222"}

def r_str(s):
    return json.dumps(s)  # a JSON string is also a valid R string

def run_r(code, script):
    with open(script, "w") as f:
        f.write(code)
    p = subprocess.run(["Rscript", script], capture_output=True, text=True)
    if p.returncode != 0:
        print(p.stderr[-3000:])
        raise SystemExit("R failed - see the message above.")
    return p.stdout

def read_csv(path, numeric):
    d = pd.read_csv(path, keep_default_na=False, dtype=str)
    for c in numeric:
        if c in d:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    return d

def role_of(sample, task):
    if task.lower() == "ntc" or re.search(NTC_RE, sample, re.I):
        return "ntc"
    return "positive" if re.search(POSITIVE_RE, sample, re.I) else "sample"

def ct_at_threshold(cycles, delta_rn, thr):
    """Cycle where the curve crosses thr for the last time from below
    (log-linear interpolation), None if it ends below. Same rule as R."""
    d = list(delta_rn)
    if not d or thr is None or pd.isna(d[-1]) or d[-1] < thr:
        return None
    below = [i for i, v in enumerate(d) if pd.isna(v) or v < thr]
    if not below:
        return float(cycles[0])
    j = below[-1]
    if pd.isna(d[j]):
        return float(cycles[j + 1])
    step = cycles[j + 1] - cycles[j]
    if d[j] > 0:
        frac = (math.log(thr) - math.log(d[j])) / (math.log(d[j + 1]) - math.log(d[j]))
    else:
        frac = (thr - d[j]) / (d[j + 1] - d[j])
    return float(cycles[j] + step * frac)

def load_data():
    wells = read_csv(os.path.join(DATA_DIR, "wells.csv"),
                     ["well_index", "ct", "amp_status", "cq_conf"])
    curves = read_csv(os.path.join(DATA_DIR, "curves.csv"),
                      ["well_index", "cycle", "rn", "delta_rn"])
    thr = read_csv(os.path.join(DATA_DIR, "thresholds.csv"), ["threshold"])
    wells["role"] = [role_of(s, t) for s, t in zip(wells["sample"], wells["task"])]
    wells["is_reference"] = wells["target"].str.contains(REFERENCE_RE, case=False, regex=True)
    by_well = {k: g.sort_values("cycle") for k, g in
               curves.groupby(["run", "well_index", "target"])}
    pos = curves["delta_rn"][curves["delta_rn"] > 0]
    top = math.log10(pos.max()) if len(pos) else 1.0
    return {"wells": wells, "curves": by_well, "top": top,
            "n_cycles": int(curves["cycle"].max()) if len(curves) else 50,
            "thresholds": {(r, t): v for r, t, v in
                           zip(thr.get("run", []), thr.get("target", []), thr.get("threshold", []))}}

def well_cts(data, run, gene, threshold=THRESHOLD):
    """Ct of every well read at the threshold; the instrument's Ct is kept in instrument_ct."""
    w = data["wells"]
    w = w[(w["run"] == run) & (w["target"] == gene)].copy()
    w["instrument_ct"] = w["ct"]
    cts = []
    for idx in w["well_index"]:
        c = data["curves"].get((run, idx, gene))
        cts.append(None if c is None else
                   ct_at_threshold(list(c["cycle"]), list(c["delta_rn"]), threshold))
    w["ct"] = pd.to_numeric(pd.Series(cts, index=w.index, dtype=object), errors="coerce")
    return w

def add_lines(fig, data, log, row=None, col=None, label=True):
    """Threshold (green, horizontal) and Ct cutoff (red, vertical)."""
    kw = {} if row is None else dict(row=row, col=col)
    fig.add_hline(y=THRESHOLD, line_dash="dot", line_color="#27ae60", line_width=2, **kw)
    fig.add_vrect(x0=CT_CUTOFF, x1=max(data["n_cycles"], CT_CUTOFF + 1), fillcolor="#888",
                  opacity=0.08, line_width=0, **kw)
    fig.add_vline(x=CT_CUTOFF, line_dash="dash", line_color="#c0392b", line_width=2, **kw)
    if label:
        # annotations on a log axis take log10 positions
        fig.add_annotation(x=0.01, xref="paper", y=math.log10(THRESHOLD) if log else THRESHOLD,
                           yref="y", text="threshold ΔRn " + fmt_thr(THRESHOLD), showarrow=False,
                           xanchor="left", yanchor="bottom", font_color="#27ae60")
        fig.add_annotation(x=CT_CUTOFF, y=1, xref="x", yref="paper", text="Ct %g" % CT_CUTOFF,
                           showarrow=False, xanchor="left", yanchor="top", font_color="#c0392b")

def add_curves(fig, data, w, run, gene, cutoff, log, row=None, col=None, legend=False):
    floor = 10 ** (data["top"] - 4)
    for _, r in w.iterrows():
        c = data["curves"].get((run, r["well_index"], gene))
        if c is None:
            continue
        detected = pd.notna(r["ct"]) and r["ct"] <= cutoff
        kind = r["role"] if r["role"] != "sample" else ("detected" if detected else "not detected")
        y = c["delta_rn"].clip(lower=floor) if log else c["delta_rn"]
        ct_txt = "no Ct (does not reach the threshold)" if pd.isna(r["ct"]) else "Ct %.2f" % r["ct"]
        fig.add_trace(go.Scatter(
            x=c["cycle"], y=y, mode="lines", showlegend=False, legendgroup=kind,
            line=dict(color=COLORS[kind], width=2 if kind == "detected" else 1.3,
                      dash={"positive": "dash", "ntc": "dot"}.get(kind, "solid")),
            hovertemplate="<b>%s</b> (%s)<br>%s<br>cycle %%{x}, ΔRn %%{y:.4g}<extra></extra>"
                          % (r["sample"], r["well"], ct_txt)), row=row, col=col)
    if legend:
        labels = {"detected": "Sample, Ct ≤ %g" % cutoff, "not detected": "Sample, no Ct ≤ %g" % cutoff,
                  "positive": "Positive control", "ntc": "No-template control"}
        for kind, label in labels.items():
            fig.add_trace(go.Scatter(x=[None], y=[None], mode="lines", name=label,
                                     legendgroup=kind,
                                     line=dict(color=COLORS[kind],
                                               dash={"positive": "dash", "ntc": "dot"}.get(kind, "solid"))))
    if log:
        fig.update_yaxes(type="log", range=[data["top"] - 4, data["top"] + 0.15], row=row, col=col)

def fmt_thr(v):
    return "{:,.0f}".format(v) if v >= 100 else "%.3g" % v

def show_table(df, max_rows=None):
    display(HTML(df.to_html(index=False, escape=True, na_rep="–", max_rows=max_rows)))

display(HTML("<h3 style='color:#2e8b57'>✔ Ready. Go to step ②.</h3>"))
'''

UPLOAD = r'''
#@title ② Upload your .eds files and look at the curves { display-mode: "form" }
#@markdown Click ▶, then **Choose Files** and select one or more `.eds` files.
#@markdown The amplification curves of every run and gene are drawn below, with the threshold (ΔRn 10,000) and the Ct cutoff (40).
#@markdown Run this step again to replace the files.
import os, shutil
from google.colab import files
from IPython.display import display, HTML

if "run_r" not in globals():
    raise SystemExit("Run step ① first.")
shutil.rmtree(EDS_DIR, ignore_errors=True)
os.makedirs(EDS_DIR)
_uploaded = files.upload()
for _name, _data in _uploaded.items():
    with open(os.path.join(EDS_DIR, _name), "wb") as f:
        f.write(_data)
    os.remove(_name) if os.path.exists(_name) else None
_eds = sorted(f for f in os.listdir(EDS_DIR) if f.lower().endswith(".eds"))
if not _eds:
    raise SystemExit("No .eds files uploaded.")

os.makedirs(DATA_DIR, exist_ok=True)
run_r(f"""library(autoqmsp)
files <- list.files({r_str(EDS_DIR)}, pattern = "\\\\.eds$", full.names = TRUE, ignore.case = TRUE)
runs <- read_eds_files(files, names = tools::file_path_sans_ext(basename(files)))
out <- {r_str(DATA_DIR)}
utils::write.csv(runs$wells, file.path(out, "wells.csv"), row.names = FALSE, na = "")
utils::write.csv(runs$curves, file.path(out, "curves.csv"), row.names = FALSE, na = "")
utils::write.csv(runs$thresholds, file.path(out, "thresholds.csv"), row.names = FALSE, na = "")
""", os.path.join(OUT_DIR, "read_script.R"))

DATA = load_data()
_w = DATA["wells"]
_summary = _w.groupby("run").agg(
    samples=("sample", lambda s: s[_w.loc[s.index, "role"] == "sample"].nunique()),
    controls=("sample", lambda s: ", ".join(sorted(set(s[_w.loc[s.index, "role"] != "sample"])))),
    genes=("target", lambda s: ", ".join(dict.fromkeys(s)))).reset_index()
display(HTML("<h3 style='color:#2e8b57'>✔ %d run(s) loaded</h3>" % len(_summary)))
show_table(_summary)

for _run in dict.fromkeys(_w["run"]):
    _genes = list(dict.fromkeys(_w.loc[_w["run"] == _run, "target"]))
    _ncol = min(3, len(_genes))
    _nrow = math.ceil(len(_genes) / _ncol)
    _fig = make_subplots(rows=_nrow, cols=_ncol, subplot_titles=_genes,
                         horizontal_spacing=0.06, vertical_spacing=0.12 if _nrow > 1 else 0.1)
    for _i, _gene in enumerate(_genes):
        _r, _c = _i // _ncol + 1, _i % _ncol + 1
        add_curves(_fig, DATA, well_cts(DATA, _run, _gene), _run, _gene, CT_CUTOFF,
                   log=True, row=_r, col=_c, legend=(_i == 0))
        add_lines(_fig, DATA, True, row=_r, col=_c, label=False)
    _fig.update_layout(title=dict(text="<b>%s</b>" % _run, font_size=15),
                       height=300 * _nrow + 90, template="plotly_white",
                       margin=dict(l=50, r=20, t=90, b=40),
                       legend=dict(orientation="h", y=1.02, yanchor="bottom", x=1, xanchor="right"))
    _fig.update_xaxes(title_text="Cycle", row=_nrow)
    _fig.show()
display(HTML("<p><span style='color:#27ae60'><b>Green dotted line</b></span>: the threshold, "
             "ΔRn %s. The Ct is read where a curve crosses it. "
             "<span style='color:#c0392b'><b>Red dashed line</b></span>: Ct %g. A gene is methylated "
             "when its curve crosses the threshold at or before it. Hover over a curve to see the sample. "
             "Go to step ③ to look at one gene at a time.</p>" % (fmt_thr(THRESHOLD), CT_CUTOFF)))
'''

EXPLORE = r'''
#@title ③ Look at each gene { display-mode: "form" }
#@markdown Click ▶, then choose a run and a gene. The graph shows the fixed lab rule:
#@markdown - **Threshold** (green, horizontal): ΔRn 10,000. The Ct is the cycle where the curve crosses it.
#@markdown - **Ct cutoff** (red, vertical): 40. A gene is **methylated** when its Ct is 40 or less.
#@markdown
#@markdown The table lists every well with its Ct. You can skip this step.
import ipywidgets as W
from IPython.display import display, HTML

if "DATA" not in globals():
    raise SystemExit("Run step ② first.")

def _viewer(data):
    wells = data["wells"]
    runs = list(dict.fromkeys(wells["run"]))
    genes_of = lambda run: list(dict.fromkeys(wells.loc[wells["run"] == run, "target"]))
    style = {"description_width": "60px"}
    run_dd = W.Dropdown(options=runs, description="Run", layout=W.Layout(width="600px"), style=style)
    gene_dd = W.Dropdown(options=genes_of(runs[0]), description="Gene",
                         layout=W.Layout(width="380px"), style=style)
    log_cb = W.Checkbox(value=True, description="Log scale", indent=False)
    info = W.HTML()
    out = W.Output()
    busy = {"on": False}

    def draw(*_):
        if busy["on"]:
            return
        run, gene = run_dd.value, gene_dd.value
        is_ref = bool(re.search(REFERENCE_RE, gene, re.I))
        w = well_cts(data, run, gene)

        fig = go.Figure()
        add_curves(fig, data, w, run, gene, CT_CUTOFF, log_cb.value, legend=True)
        add_lines(fig, data, log_cb.value)
        fig.update_layout(title="<b>%s</b>%s — %s" % (gene, " (reference gene)" if is_ref else "", run),
                          template="plotly_white", height=480, margin=dict(l=60, r=20, t=60, b=50),
                          xaxis_title="Cycle", yaxis_title="ΔRn",
                          legend=dict(orientation="h", y=-0.18, x=0))

        within = lambda d: d[d["ct"] <= CT_CUTOFF]
        samples = w[w["role"] == "sample"]
        ntc, pc = w[w["role"] == "ntc"], w[w["role"] == "positive"]
        word = "have a reference Ct ≤ %g (DNA OK)" if is_ref else "are methylated (Ct ≤ %g)"
        msgs = ["<b>%d of %d</b> sample wells %s." % (len(within(samples)), len(samples), word % CT_CUTOFF)]
        if len(ntc):
            hit = within(ntc)
            msgs.append("<span style='color:#c0392b'>⚠ No-template control amplified: %s</span>"
                        % ", ".join("%s (Ct %.1f)" % (s, c) for s, c in zip(hit["sample"], hit["ct"]))
                        if len(hit) else "<span style='color:#2e8b57'>✔ No-template controls stay negative.</span>")
        if len(pc):
            msgs.append("<span style='color:#2e8b57'>✔ Positive controls amplify.</span>" if len(within(pc))
                        else "<span style='color:#c0392b'>⚠ No positive control reaches the threshold by Ct %g.</span>" % CT_CUTOFF)
        info.value = "<br>".join(msgs)

        role_name = {"sample": "Sample", "positive": "Positive control", "ntc": "No-template control"}
        tab = pd.DataFrame({
            "Sample": w["sample"], "Well": w["well"], "Type": w["role"].map(role_name),
            "Ct at ΔRn %s" % fmt_thr(THRESHOLD): w["ct"].round(2),
            "Instrument Ct": w["instrument_ct"].where(w["instrument_ct"] < data["n_cycles"]).round(2)})
        tab["Result"] = [("Ct ≤ %g" % CT_CUTOFF) if pd.notna(c) and c <= CT_CUTOFF else
                         ("Ct above %g" % CT_CUTOFF if pd.notna(c) else "does not reach the threshold")
                         for c in w["ct"]]
        tab["Type"] = pd.Categorical(tab["Type"], list(role_name.values()))
        tab = tab.sort_values(["Type", tab.columns[3]], na_position="last")
        with out:
            out.clear_output(wait=True)
            fig.show()
            show_table(tab)

    def on_run(change):
        busy["on"] = True
        genes = genes_of(run_dd.value)
        keep = gene_dd.value if gene_dd.value in genes else genes[0]
        gene_dd.options = genes
        gene_dd.value = keep
        busy["on"] = False
        draw()

    run_dd.observe(on_run, "value")
    gene_dd.observe(draw, "value")
    log_cb.observe(draw, "value")
    draw()
    return W.VBox([run_dd, W.HBox([gene_dd, log_cb]), info, out]), draw

_ui, _draw = _viewer(DATA)
display(_ui)
'''

ANALYSE = r'''
#@title ④ Report { display-mode: "form" }
#@markdown Click ▶. The report and a Ct chart appear below, and the report and Excel file download.
#@markdown
#@markdown **Rule (fixed):** the Ct is read at ΔRn **10,000**; a gene is **methylated** when its Ct is **40 or less**.
#@markdown Every sample is called **Potential cancer**, **Not determined** (repeat) or **Low risk**.

#@markdown ### Cancer risk decision
#@markdown A sample is **Potential cancer** when at least this many genes are methylated:
min_methylated_genes = 1  #@param {type:"integer"}
#@markdown Genes that count (comma separated; empty = all except the reference gene):
panel_genes = ""  #@param {type:"string"}

#@markdown ### Controls (leave empty to detect them from the sample names)
reference_gene = "B ACTIN"  #@param {type:"string"}
no_template_controls = ""  #@param {type:"string"}
positive_controls = ""  #@param {type:"string"}
#@markdown Warn "low DNA input" when the reference gene Ct is above:
reference_ct_low_input = 35  #@param {type:"number"}

#@markdown ---
show_r_code = False  #@param {type:"boolean"}
download_files = True  #@param {type:"boolean"}

import os, re, zlib
from IPython.display import display, HTML

if "run_r" not in globals():
    raise SystemExit("Run step ① first.")

def _names(text):
    return [t.strip() for t in text.split(",") if t.strip()]

def _exact_regex(names):
    return "^(" + "|".join(re.escape(n) for n in names) + ")$"

def build_r_code(eds_dir, out_dir):
    ref = reference_gene.strip()
    args = {
        "reference": r_str(_exact_regex([ref])) if ref else "NULL",
        "threshold": "%g" % THRESHOLD,
        "ct_cutoff": "%g" % CT_CUTOFF,
        "ref_ct_max": "%g" % CT_CUTOFF,
        "ref_ct_warn": reference_ct_low_input,
        "min_methylated_genes": int(min_methylated_genes),
        "panel": "c(%s)" % ", ".join(r_str(g) for g in _names(panel_genes)) if _names(panel_genes) else "NULL",
    }
    if _names(no_template_controls):
        args["ntc"] = r_str(_exact_regex(_names(no_template_controls)))
    if _names(positive_controls):
        args["positive"] = r_str(_exact_regex(_names(positive_controls)))
    arg_text = ",\n".join("  %s = %s" % kv for kv in args.items())
    return f"""library(autoqmsp)
files <- list.files({r_str(eds_dir)}, pattern = "\\\\.eds$", full.names = TRUE, ignore.case = TRUE)
runs <- read_eds_files(files, names = tools::file_path_sans_ext(basename(files)))
res <- analyze_qmsp(
  runs,
{arg_text}
)
dir.create({r_str(out_dir)}, showWarnings = FALSE)
write_report(res, file.path({r_str(out_dir)}, "qmsp_report.html"))
export_results(res, file.path({r_str(out_dir)}, "qmsp_results.xlsx"))
utils::write.csv(res$results, file.path({r_str(out_dir)}, "results.csv"), row.names = FALSE, na = "")
print(res)
"""

def run_analysis(eds_dir=EDS_DIR, out_dir=OUT_DIR):
    if not os.path.isdir(eds_dir) or not any(f.lower().endswith(".eds") for f in os.listdir(eds_dir)):
        raise SystemExit("No .eds files found. Run step ② first.")
    code = build_r_code(eds_dir, out_dir)
    os.makedirs(out_dir, exist_ok=True)
    run_r(code, os.path.join(out_dir, "analysis_script.R"))
    return code

def ct_figure(results):
    r = results[results["role"] == "sample"].copy()
    genes = list(dict.fromkeys(r["target"]))
    top = max(46.0, math.ceil(r["ct"].max() + 1) if r["ct"].notna().any() else 46.0)
    no_ct = top + 2  # row for genes whose curve never reached the threshold
    colors = {"Methylated": "#c0392b", "Unmethylated": "#7fa7c9", "Not determined": "#f0b429"}
    fig = go.Figure()
    for call, color in colors.items():
        d = r[r["call"] == call]
        if not len(d):
            continue
        x = [genes.index(g) + (zlib.crc32(s.encode()) % 1000 / 1000 - 0.5) * 0.5
             for g, s in zip(d["target"], d["sample"] + d["run"])]
        text = ["<b>%s</b><br>%s<br>%s<br>%s" % (s, g, "no Ct" if pd.isna(c) else "Ct %.2f" % c, rn)
                for s, g, c, rn in zip(d["sample"], d["target"], d["ct"], d["run"])]
        fig.add_trace(go.Scatter(x=x, y=d["ct"].fillna(no_ct), mode="markers", name=call, text=text,
                                 hovertemplate="%{text}<extra></extra>",
                                 marker=dict(color=color, size=10, line=dict(color="white", width=1))))
    fig.add_hrect(y0=CT_CUTOFF, y1=no_ct + 1.5, fillcolor="#888", opacity=0.08, line_width=0)
    fig.add_hline(y=CT_CUTOFF, line_dash="dash", line_color="#c0392b", line_width=2,
                  annotation_text="Ct %g: methylated below this line" % CT_CUTOFF,
                  annotation_position="bottom right", annotation_font_color="#c0392b")
    ticks = list(range(15, int(top) + 1, 5))
    fig.update_yaxes(range=[no_ct + 1.5, min(15, (r["ct"].min() - 2) if r["ct"].notna().any() else 15)],
                     title="Ct at ΔRn %s" % fmt_thr(THRESHOLD),
                     tickvals=ticks + [no_ct], ticktext=[str(t) for t in ticks] + ["no Ct"])
    fig.update_xaxes(tickvals=list(range(len(genes))), ticktext=genes, range=[-0.6, len(genes) - 0.4])
    fig.update_layout(title="<b>Ct per gene</b> (each dot is a sample)",
                      template="plotly_white", height=460, legend_title_text="Call",
                      margin=dict(l=60, r=20, t=60, b=40))
    return fig

_code = run_analysis()
if show_r_code:
    print(_code)
_res = read_csv(os.path.join(OUT_DIR, "results.csv"), ["ct"])
with open(os.path.join(OUT_DIR, "qmsp_report.html")) as f:
    display(HTML(f.read()))
ct_figure(_res).show()
if download_files:
    from google.colab import files
    files.download(os.path.join(OUT_DIR, "qmsp_report.html"))
    files.download(os.path.join(OUT_DIR, "qmsp_results.xlsx"))
'''

INTRO = """
# qMSP analysis: upload .eds files and get a report

1. Click ▶ on **step ①** to install the tool. You only need to do this once per session.
2. Click ▶ on **step ②** and choose your QuantStudio `.eds` files. The amplification curves are drawn.
3. Click ▶ on **step ③** to look at one gene at a time, with a table of every well.
4. Click ▶ on **step ④**. The report and a Ct chart appear, and the report (HTML) and the Excel file download.

**The rule:** the Ct is read where the amplification curve crosses **ΔRn = 10,000**. A gene is **methylated** when its **Ct is 40 or less**. Each sample is then called **Potential cancer** (at least one methylated gene), **Not determined** (the result could not be trusted, e.g. the reference gene or a control failed; repeat the sample) or **Low risk**.

To change a setting, edit it and click ▶ on the step again. To see the code behind a step, double-click the step.

*For research use only. The results are not a diagnosis.*
"""

nb = {
    "cells": [md(INTRO), form(INSTALL), form(UPLOAD), form(EXPLORE), form(ANALYSE)],
    "metadata": {"colab": {"provenance": [], "toc_visible": False},
                 "kernelspec": {"name": "python3", "display_name": "Python 3"},
                 "language_info": {"name": "python"}},
    "nbformat": 4, "nbformat_minor": 0,
}
here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, "autoqmsp_colab.ipynb"), "w") as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)
