"""Builds colab/autoqmsp_cutoffs_colab.ipynb. Run: python3 colab/build_cutoff_notebook.py"""
import textwrap

from build_notebook import INSTALL, UPLOAD_LOAD, UPLOAD_CURVES, form, md, notebook, write

UPLOAD = r'''
#@title ② Upload your .eds files { display-mode: "form" }
#@markdown Click ▶, then **Choose Files** and select all the `.eds` files of the samples you know the diagnosis of.
#@markdown Samples with the same name in different runs are pooled.
#@markdown Tick to also draw the amplification curves of every run:
show_curves = False  #@param {type:"boolean"}
''' + UPLOAD_LOAD + "\nif show_curves:\n" + textwrap.indent(
    UPLOAD_CURVES.replace("Go to step ③ to look at one gene at a time.", "Go to step ③."), "    ")

LABEL = r'''
#@title ③ Tell me which samples are which { display-mode: "form" }
#@markdown Click ▶. For every sample, choose its group in the list:
#@markdown - **Cancer**: a sample from a patient with cancer
#@markdown - **Control**: a sample from a person without cancer (healthy / normal)
#@markdown - **Positive control**: methylated DNA, e.g. H460, A549, HT29
#@markdown - **NK**: no-template (water) control
#@markdown - **Leave out**: not used
#@markdown
#@markdown NK and positive controls are pre-filled from their names: check them. Your choices are kept when you run the step again.
#@markdown
#@markdown Tick to fill the groups from a labels file (CSV with the columns `sample` and `group`), for example the `qmsp_labels.csv` saved in step ⑤ before:
upload_labels_file = False  #@param {type:"boolean"}
import ipywidgets as W
from IPython.display import display, HTML

if "DATA" not in globals():
    raise SystemExit("Run step ② first.")

GROUPS = [("— choose —", ""), ("Cancer", "cancer"), ("Control (no cancer)", "control"),
          ("Positive control", "positive"), ("NK (water)", "ntc"), ("Leave out", "exclude")]
GROUP_NAME = {v: k for k, v in GROUPS}

def norm_group(text):
    t = str(text).strip().lower()
    for pattern, group in [(r"^(cancer|tumou?r|case|patient)", "cancer"),
                           (r"^(control|healthy|normal|non.?cancer|benign)", "control"),
                           (r"^(positive|pc|pos)", "positive"),
                           (r"^(ntc|nk|no.?template|water|blank)", "ntc"),
                           (r"^(exclude|leave|skip|none)", "exclude")]:
        if re.match(pattern, t):
            return group
    return ""

LABELS = globals().get("LABELS") or {}
if upload_labels_file:
    from google.colab import files
    _up = files.upload()
    for _n, _b in _up.items():
        _p = os.path.join(OUT_DIR, "uploaded_labels.csv")
        os.makedirs(OUT_DIR, exist_ok=True)
        with open(_p, "wb") as f:
            f.write(_b)
        os.remove(_n) if os.path.exists(_n) else None
        _d = read_csv(_p, [])
        if not {"sample", "group"} <= set(_d.columns):
            raise SystemExit("The labels file needs the columns sample and group.")
        LABELS.update({s: norm_group(g) for s, g in zip(_d["sample"], _d["group"])})

def _sample_table(data):
    w = data["wells"]
    rows = []
    for s in dict.fromkeys(w["sample"]):
        d = w[w["sample"] == s]
        genes = [g for g in dict.fromkeys(d["target"]) if not re.search(REFERENCE_RE, g, re.I)]
        n_runs = d["run"].nunique()
        guess = {"ntc": "ntc", "positive": "positive"}.get(d["role"].iloc[0], "")
        rows.append((s, "%d gene%s · %d run%s" % (len(genes), "" if len(genes) == 1 else "s",
                                                  n_runs, "" if n_runs == 1 else "s"), guess))
    return rows

def _labeller(rows, labels):
    for s, _, guess in rows:
        labels.setdefault(s, guess)
    summary = W.HTML()
    def update(*_):
        n = {}
        for s, _, _ in rows:
            n[labels.get(s, "")] = n.get(labels.get(s, ""), 0) + 1
        parts = ["<b>%s</b>: %d" % (GROUP_NAME[g], n[g]) for g in
                 ["cancer", "control", "positive", "ntc", "exclude"] if n.get(g)]
        if n.get(""):
            parts.append("<b style='color:#c0392b'>no group yet: %d</b>" % n[""])
        summary.value = "<p style='font-size:15px'>%s</p>" % " · ".join(parts)
    boxes = []
    dds = {}
    for s, info, _ in rows:
        dd = W.Dropdown(options=GROUPS, value=labels.get(s, ""), layout=W.Layout(width="210px"))
        def on_change(change, s=s):
            labels[s] = change["new"]
            update()
        dd.observe(on_change, "value")
        dds[s] = dd
        boxes.append(W.HBox([W.HTML("<b>%s</b>" % s, layout=W.Layout(width="300px")),
                             W.HTML("<span style='color:#777'>%s</span>" % info,
                                    layout=W.Layout(width="150px")), dd]))
    fill = W.Dropdown(options=GROUPS[1:], value="control", layout=W.Layout(width="210px"))
    apply = W.Button(description="Apply", layout=W.Layout(width="80px"))
    def on_apply(_):
        for s, dd in dds.items():
            if dd.value == "":
                dd.value = fill.value
    apply.on_click(on_apply)
    update()
    return W.VBox([summary,
                   W.HBox([W.HTML("Give every sample that has no group yet the group:&nbsp;"), fill, apply]),
                   W.HTML("<hr style='margin:6px 0'>")] + boxes)

display(_labeller(_sample_table(DATA), LABELS))
display(HTML("<p>When every sample has a group, go to step ④.</p>"))
'''

FIND = r'''
#@title ④ Find the cutoff for each gene { display-mode: "form" }
#@markdown Click ▶. For every gene, the Ct of the cancer samples is compared with the Ct of the controls (Ct read at ΔRn 10,000).
#@markdown The suggested cutoff is the Ct that separates them best: a sample is **methylated** when its Ct is at or below it.
#@markdown Then pick a gene and move the slider to see what another cutoff would give. Running this step again resets the cutoffs to the suggested ones.

#@markdown ### How to choose the cutoff
choose_cutoff_by = "Best balance (Youden)"  #@param ["Best balance (Youden)", "High specificity: few false cancer results", "High sensitivity: few missed cancers"]
#@markdown Minimum specificity or sensitivity to reach (for the "High …" choices):
minimum = 0.9  #@param {type:"slider", min:0.5, max:1, step:0.05}

#@markdown ### Samples
reference_gene = "B ACTIN"  #@param {type:"string"}
#@markdown Leave out samples whose reference gene has no Ct of 40 or less (failed DNA):
leave_out_failed_dna = True  #@param {type:"boolean"}

#@markdown ### All genes together
#@markdown A sample counts as cancer when at least this many genes are methylated:
min_methylated_genes = 1  #@param {type:"integer"}
show_r_code = False  #@param {type:"boolean"}
import zlib
import ipywidgets as W
from IPython.display import display, HTML

if "LABELS" not in globals():
    raise SystemExit("Run step ③ first.")

METHODS = {"Best balance (Youden)": "youden",
           "High specificity: few false cancer results": "specificity",
           "High sensitivity: few missed cancers": "sensitivity"}
GROUP_STYLE = {"cancer": ("Cancer", "#c0392b"), "control": ("Control", "#1f6fb4"),
               "positive": ("Positive control", "#7d3c98"), "ntc": ("NK", "#555555")}
LABELS_FILE = os.path.join(OUT_DIR, "qmsp_labels.csv")
CUT_DIR = os.path.join(OUT_DIR, "cutoffs")

def write_labels(labels, data):
    present = list(dict.fromkeys(data["wells"]["sample"]))
    os.makedirs(OUT_DIR, exist_ok=True)
    pd.DataFrame({"sample": present,
                  "group": [labels.get(s) or "exclude" for s in present]}).to_csv(LABELS_FILE, index=False)

def find_r_code(args, chosen=None):
    ref = args["reference_gene"].strip()
    ref_re = r_str("^(" + re.escape(ref) + ")$") if ref else "NULL"
    chosen_txt = "NULL"
    if chosen:
        chosen_txt = "c(%s)" % ", ".join("`%s` = %r" % (g.replace("`", ""), float(c))
                                         for g, c in chosen.items() if c is not None and not pd.isna(c))
        chosen_txt = "NULL" if chosen_txt == "c()" else chosen_txt
    return f"""library(autoqmsp)
files <- list.files({r_str(EDS_DIR)}, pattern = "\\\\.eds$", full.names = TRUE, ignore.case = TRUE)
runs <- read_eds_files(files, names = tools::file_path_sans_ext(basename(files)))
groups <- utils::read.csv({r_str(LABELS_FILE)}, stringsAsFactors = FALSE, check.names = FALSE)
cut <- find_cutoffs(
  runs, groups,
  method = "{args['method']}", target = {args['minimum']},
  threshold = {THRESHOLD:g}, reference = {ref_re}, ref_ct_max = {CT_CUTOFF:g},
  exclude_failed_reference = {"TRUE" if args['leave_out_failed_dna'] else "FALSE"},
  cutoffs = {chosen_txt}
)
out <- {r_str(CUT_DIR)}
dir.create(out, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(cut$cutoffs, file.path(out, "cutoffs.csv"), row.names = FALSE, na = "")
utils::write.csv(cut$samples, file.path(out, "samples.csv"), row.names = FALSE, na = "")
utils::write.csv(cut$roc, file.path(out, "roc.csv"), row.names = FALSE, na = "")
print(cut)
"""

def load_found():
    cut = read_csv(os.path.join(CUT_DIR, "cutoffs.csv"),
                   ["cutoff", "suggested_cutoff", "sensitivity", "specificity", "auc", "n_cancer",
                    "n_control", "tp", "fn", "fp", "tn", "positive_ct", "ntc_ct"])
    samp = read_csv(os.path.join(CUT_DIR, "samples.csv"), ["ct", "ref_ct", "n_wells", "n_with_ct"])
    samp["used"] = samp["used"].str.upper() == "TRUE"
    roc = read_csv(os.path.join(CUT_DIR, "roc.csv"), ["cutoff", "sensitivity", "specificity"])
    return cut, samp, roc

def pct(a, b):
    return "–" if not b else "%d of %d (%.0f%%)" % (a, b, 100.0 * a / b)

def gene_stats(samp, gene, cutoff):
    s = samp[(samp["target"] == gene) & samp["used"]]
    hit = (s["ct"].notna() & (s["ct"] <= cutoff)) if cutoff is not None else pd.Series(False, index=s.index)
    ca, co = s["group"] == "cancer", s["group"] == "control"
    return {"tp": int((hit & ca).sum()), "n_ca": int(ca.sum()),
            "fp": int((hit & co).sum()), "n_co": int(co.sum())}

def panel_text(samp, chosen, n_min):
    s = samp[samp["used"] & samp["target"].isin(chosen)].copy()
    s["hit"] = [pd.notna(c) and chosen.get(g) is not None and c <= chosen[g]
                for g, c in zip(s["target"], s["ct"])]
    per = s.groupby(["sample", "group"]).agg(
        n=("hit", "sum"), genes=("target", lambda t: ", ".join(t[s.loc[t.index, "hit"]]))).reset_index()
    per["called"] = per["n"] >= n_min
    ca, co = per[per["group"] == "cancer"], per[per["group"] == "control"]
    missed = ca[~ca["called"]]["sample"].tolist()
    false = ["%s (%s)" % (a, b) for a, b in zip(co[co["called"]]["sample"], co[co["called"]]["genes"])]
    return ("<div style='background:#f4f6f8;padding:10px 12px;border-radius:6px;margin-top:10px'>"
            "<b>All genes together</b>, with the cutoffs as they are now; a sample counts as cancer when at least "
            "<b>%d</b> gene%s %s methylated:<br>"
            "Cancers found: <b>%s</b> (sensitivity) · Controls correctly negative: <b>%s</b> (specificity)<br>"
            "<span style='color:#777'>Cancers missed: %s<br>Controls called cancer: %s</span></div>"
            % (n_min, "" if n_min == 1 else "s", "is" if n_min == 1 else "are",
               pct(int(ca["called"].sum()), len(ca)), pct(int((~co["called"]).sum()), len(co)),
               ", ".join(missed) or "none", ", ".join(false) or "none"))

def cutoff_figure(samp, roc, gene, cutoff, auc_value):
    s = samp[(samp["target"] == gene) & (samp["group"].isin(GROUP_STYLE))]
    top = max(46.0, math.ceil(s["ct"].max() + 1) if s["ct"].notna().any() else 46.0)
    no_ct = top + 2
    groups = [g for g in GROUP_STYLE if (s["group"] == g).any()]
    fig = make_subplots(rows=1, cols=2, column_widths=[0.62, 0.38], horizontal_spacing=0.1,
                        subplot_titles=("Ct of every sample", "ROC curve (AUC %s)" %
                                        ("–" if pd.isna(auc_value) else "%.2f" % auc_value)))
    for i, g in enumerate(groups):
        d = s[s["group"] == g]
        label, color = GROUP_STYLE[g]
        for used in (True, False):
            dd = d[d["used"] == used] if g in ("cancer", "control") else (d if used else d.iloc[0:0])
            if not len(dd):
                continue
            x = [i + (zlib.crc32(n.encode()) % 1000 / 1000 - 0.5) * 0.5 for n in dd["sample"]]
            text = ["<b>%s</b><br>%s<br>%s%s" % (n, label, "no Ct" if pd.isna(c) else "Ct %.2f" % c,
                                                 "" if used else "<br>left out: failed DNA")
                    for n, c in zip(dd["sample"], dd["ct"])]
            fig.add_trace(go.Scatter(
                x=x, y=dd["ct"].fillna(no_ct), mode="markers", text=text, showlegend=False,
                hovertemplate="%{text}<extra></extra>",
                marker=dict(color="white" if not used else color, size=10,
                            line=dict(color=color, width=1.5 if not used else 1),
                            symbol="circle" if used else "circle-open")), row=1, col=1)
    if cutoff is not None:
        fig.add_shape(type="rect", x0=-0.5, x1=len(groups) - 0.5, y0=10, y1=cutoff, fillcolor="#c0392b",
                      opacity=0.06, line_width=0, row=1, col=1)
        fig.add_shape(type="line", x0=-0.5, x1=len(groups) - 0.5, y0=cutoff, y1=cutoff,
                      line=dict(color="#c0392b", dash="dash", width=2), row=1, col=1)
        fig.add_annotation(x=len(groups) - 0.5, y=cutoff, text="cutoff %.1f" % cutoff, showarrow=False,
                           xanchor="right", yanchor="bottom", font_color="#c0392b", row=1, col=1)
    ticks = list(range(15, int(top) + 1, 5))
    lo = min(15, (s["ct"].min() - 2) if s["ct"].notna().any() else 15)
    fig.update_yaxes(range=[no_ct + 1.5, lo], tickvals=ticks + [no_ct],
                     ticktext=[str(t) for t in ticks] + ["no Ct"], title="Ct at ΔRn %s" % fmt_thr(THRESHOLD),
                     row=1, col=1)
    fig.update_xaxes(tickvals=list(range(len(groups))), ticktext=[GROUP_STYLE[g][0] for g in groups],
                     range=[-0.5, len(groups) - 0.5], row=1, col=1)

    r = roc[roc["target"] == gene]
    if len(r):
        fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines", showlegend=False, hoverinfo="skip",
                                 line=dict(color="#bbb", dash="dot")), row=1, col=2)
        fig.add_trace(go.Scatter(
            x=list(1 - r["specificity"]) + [1], y=list(r["sensitivity"]) + [1], mode="lines+markers",
            line=dict(color="#555", shape="hv"), marker=dict(size=5, color="#555"), showlegend=False,
            text=["Ct ≤ %s" % ("none" if pd.isna(c) else "%.1f" % c) for c in r["cutoff"]] + ["all"],
            hovertemplate="%{text}<br>sensitivity %{y:.0%}<br>1 − specificity %{x:.0%}<extra></extra>"),
            row=1, col=2)
        st = gene_stats(samp, gene, cutoff)
        if st["n_ca"] and st["n_co"]:
            fig.add_trace(go.Scatter(x=[st["fp"] / st["n_co"]], y=[st["tp"] / st["n_ca"]], mode="markers",
                                     marker=dict(color="#c0392b", size=14), showlegend=False,
                                     hovertemplate="your cutoff<extra></extra>"), row=1, col=2)
    fig.update_xaxes(range=[-0.03, 1.03], title="1 − specificity (controls called methylated)",
                     tickformat=".0%", row=1, col=2)
    fig.update_yaxes(range=[-0.03, 1.03], title="sensitivity (cancers called methylated)",
                     tickformat=".0%", row=1, col=2)
    fig.update_layout(title="<b>%s</b>" % gene, template="plotly_white", height=470,
                      margin=dict(l=60, r=20, t=80, b=50))
    return fig

def _explorer(cut, samp, roc, chosen, n_min):
    genes = list(cut["target"])
    gene_dd = W.Dropdown(options=genes, description="Gene", layout=W.Layout(width="360px"),
                         style={"description_width": "50px"})
    slider = W.FloatSlider(min=15, max=max(DATA["n_cycles"], 41), step=0.1, description="Cutoff",
                           continuous_update=False, readout_format=".1f", layout=W.Layout(width="520px"),
                           style={"description_width": "50px"})
    reset = W.Button(description="Use suggested", icon="undo", layout=W.Layout(width="150px"))
    info, out, panel = W.HTML(), W.Output(), W.HTML()
    busy = {"on": False}

    def row(gene):
        return cut[cut["target"] == gene].iloc[0]

    def draw(*_):
        gene = gene_dd.value
        c = chosen.get(gene)
        r = row(gene)
        st = gene_stats(samp, gene, c)
        sug = r["suggested_cutoff"]
        msgs = []
        if not st["n_ca"] or not st["n_co"]:
            msgs.append("<span style='color:#c0392b'>⚠ This gene needs at least one cancer and one "
                        "control sample.</span>")
        else:
            msgs.append("At <b>Ct ≤ %s</b>: cancers methylated <b>%s</b> (sensitivity) · controls not "
                        "methylated <b>%s</b> (specificity)."
                        % ("–" if c is None else "%.1f" % c, pct(st["tp"], st["n_ca"]),
                           pct(st["n_co"] - st["fp"], st["n_co"])))
            msgs.append("<span style='color:#777'>Suggested cutoff: %s%s</span>" % (
                "none (no cutoff separates the groups)" if pd.isna(sug) else "%.1f" % sug,
                "" if c is None or pd.isna(sug) or abs(c - sug) < 0.05 else " · you moved it"))
        if c is not None and pd.notna(r["positive_ct"]) and r["positive_ct"] > c:
            msgs.append("<span style='color:#c0392b'>⚠ The positive control (Ct %.1f) is above this "
                        "cutoff.</span>" % r["positive_ct"])
        if c is not None and pd.notna(r["ntc_ct"]) and r["ntc_ct"] <= c:
            msgs.append("<span style='color:#c0392b'>⚠ The NK control (Ct %.1f) is at or below this "
                        "cutoff.</span>" % r["ntc_ct"])
        if r["note"] and "few samples" in r["note"]:
            msgs.append("<span style='color:#b9770e'>⚠ %s.</span>" % r["note"].split(";")[0].capitalize())
        info.value = "<br>".join(msgs)
        panel.value = panel_text(samp, chosen, n_min)
        with out:
            out.clear_output(wait=True)
            cutoff_figure(samp, roc, gene, c, r["auc"]).show()

    def sync():
        busy["on"] = True
        c = chosen.get(gene_dd.value)
        slider.value = c if c is not None else slider.min
        busy["on"] = False

    def on_gene(_):
        sync()
        draw()

    def on_slide(change):
        if busy["on"]:
            return
        chosen[gene_dd.value] = round(float(change["new"]), 1)
        draw()

    def on_reset(_):
        sug = row(gene_dd.value)["suggested_cutoff"]
        chosen[gene_dd.value] = None if pd.isna(sug) else float(sug)
        sync()
        draw()

    gene_dd.observe(on_gene, "value")
    slider.observe(on_slide, "value")
    reset.on_click(on_reset)
    sync()
    draw()
    return W.VBox([gene_dd, W.HBox([slider, reset]), info, out, panel])

write_labels(LABELS, DATA)
FIND_ARGS = {"method": METHODS[choose_cutoff_by], "minimum": minimum, "reference_gene": reference_gene,
             "leave_out_failed_dna": leave_out_failed_dna, "min_methylated_genes": int(min_methylated_genes)}
_n = pd.Series([LABELS.get(s, "") for s in dict.fromkeys(DATA["wells"]["sample"])]).value_counts()
if not _n.get("cancer") or not _n.get("control"):
    raise SystemExit("Give at least one sample the group Cancer and one the group Control in step ③.")
_code = find_r_code(FIND_ARGS)
if show_r_code:
    print(_code)
run_r(_code, os.path.join(OUT_DIR, "cutoff_script.R"))
CUT, SAMP, ROC = load_found()
CHOSEN = {g: (None if pd.isna(c) else float(c)) for g, c in zip(CUT["target"], CUT["suggested_cutoff"])}

_tab = pd.DataFrame({
    "Gene": CUT["target"],
    "Suggested cutoff (Ct ≤)": CUT["suggested_cutoff"].round(1),
    "Sensitivity": [pct(a, a + b) if pd.notna(a) else "–" for a, b in zip(CUT["tp"], CUT["fn"])],
    "Specificity": [pct(b, a + b) if pd.notna(a) else "–" for a, b in zip(CUT["fp"], CUT["tn"])],
    "AUC": CUT["auc"].round(2),
    "Cancer / control samples": ["%d / %d" % (a, b) for a, b in zip(CUT["n_cancer"], CUT["n_control"])],
    "Note": CUT["note"]})
display(HTML("<h3>Suggested cutoff per gene</h3><p style='color:#777'>Sensitivity: cancers with Ct at or "
             "below the cutoff. Specificity: controls above it or with no Ct. AUC: how well the gene "
             "separates cancer from control (1 = perfectly, 0.5 = not at all).</p>"))
show_table(_tab)
display(HTML("<h3>Look at each gene</h3>"))
display(_explorer(CUT, SAMP, ROC, CHOSEN, int(min_methylated_genes)))
display(HTML("<p>When you are happy with the cutoffs, go to step ⑤ to save them.</p>"))
'''

SAVE = r'''
#@title ⑤ Save the cutoffs { display-mode: "form" }
#@markdown Click ▶. Three files download:
#@markdown - **qmsp_cutoffs.csv**: the cutoff of every gene. Load it in step ④ of the main qMSP notebook (tick *use_cutoffs_file*).
#@markdown - **qmsp_cutoffs.xlsx**: all the numbers: the cutoffs, every sample's Ct and every cutoff tried.
#@markdown - **qmsp_labels.csv**: the groups you chose, to load in step ③ next time.
from IPython.display import display, HTML

if "CHOSEN" not in globals():
    raise SystemExit("Run step ④ first.")
write_labels(LABELS, DATA)
_sug = dict(zip(CUT["target"], CUT["suggested_cutoff"]))
_moved = {g: c for g, c in CHOSEN.items()
          if c is not None and (pd.isna(_sug.get(g)) or abs(c - _sug[g]) >= 0.05)}
_code = find_r_code(FIND_ARGS, _moved) + f"""
write_cutoffs(cut, file.path(out, "qmsp_cutoffs.csv"))
export_cutoffs(cut, file.path(out, "qmsp_cutoffs.xlsx"))
"""
run_r(_code, os.path.join(OUT_DIR, "save_script.R"))
_final = read_csv(os.path.join(CUT_DIR, "qmsp_cutoffs.csv"),
                  ["ct_cutoff", "suggested_cutoff", "sensitivity", "specificity", "auc"])
show_table(pd.DataFrame({
    "Gene": _final["gene"], "Cutoff (Ct ≤)": _final["ct_cutoff"].round(1),
    "Suggested": _final["suggested_cutoff"].round(1),
    "Sensitivity": ["–" if pd.isna(v) else "%.0f%%" % (100 * v) for v in _final["sensitivity"]],
    "Specificity": ["–" if pd.isna(v) else "%.0f%%" % (100 * v) for v in _final["specificity"]],
    "AUC": _final["auc"].round(2)}))
display(HTML(panel_text(SAMP, CHOSEN, FIND_ARGS["min_methylated_genes"])))
from google.colab import files
for _f in ["qmsp_cutoffs.csv", "qmsp_cutoffs.xlsx"]:
    files.download(os.path.join(CUT_DIR, _f))
files.download(LABELS_FILE)
'''

INTRO = """
# qMSP: find the Ct cutoff of each gene

Use this notebook with samples whose diagnosis you know. It finds, for every gene, the Ct cutoff that best separates **cancer** samples from **controls** (people without cancer). The Ct is read where the amplification curve crosses **ΔRn = 10,000**, as in the main notebook.

1. Click ▶ on **step ①** to install the tool (once per session).
2. Click ▶ on **step ②** and choose your `.eds` files.
3. Click ▶ on **step ③** and choose the group of every sample: cancer, control, positive control or NK.
4. Click ▶ on **step ④**. You get a suggested cutoff for each gene, with its sensitivity and specificity. Pick a gene and move the slider to try another cutoff.
5. Click ▶ on **step ⑤** to download the cutoffs. Load `qmsp_cutoffs.csv` in step ④ of the main qMSP notebook.

The more samples you have in each group, the more reliable the cutoffs. With fewer than 10 cancers or 10 controls a gene gets a warning. Check new cutoffs on samples that were not used to choose them.

*For research use only. The results are not a diagnosis.*
"""

if __name__ == "__main__":
    write(notebook([md(INTRO), form(INSTALL), form(UPLOAD), form(LABEL), form(FIND), form(SAVE)]),
          "autoqmsp_cutoffs_colab.ipynb")
