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
import subprocess, shutil, sys
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
display(HTML("<h3 style='color:#2e8b57'>✔ Ready. Go to step ②.</h3>"))
'''

UPLOAD = r'''
#@title ② Upload your .eds files { display-mode: "form" }
#@markdown Click ▶, then **Choose Files** and select one or more `.eds` files.
#@markdown Run this step again to replace the files.
import os, shutil
from google.colab import files
from IPython.display import display, HTML

EDS_DIR = "/content/eds"
shutil.rmtree(EDS_DIR, ignore_errors=True)
os.makedirs(EDS_DIR)
_uploaded = files.upload()
for _name, _data in _uploaded.items():
    with open(os.path.join(EDS_DIR, _name), "wb") as f:
        f.write(_data)
    os.remove(_name) if os.path.exists(_name) else None
_eds = sorted(f for f in os.listdir(EDS_DIR) if f.lower().endswith(".eds"))
if _eds:
    display(HTML("<h3 style='color:#2e8b57'>✔ %d file(s) ready: %s</h3><p>Go to step ③.</p>"
                 % (len(_eds), ", ".join(_eds))))
else:
    display(HTML("<h3 style='color:#c0392b'>No .eds files uploaded.</h3>"))
'''

ANALYSE = r'''
#@title ③ Settings and report { display-mode: "form" }
#@markdown Change the settings if you need to, then click ▶. The report appears below, and the report and Excel file download automatically.

#@markdown ### When is a gene methylated?
ct_cutoff = 40  #@param {type:"slider", min:25, max:50, step:0.5}
beta_cutoff = 0  #@param {type:"slider", min:0, max:1, step:0.01}
#@markdown *Beta = methylation level from 0 (none) to 1 (as methylated as the positive control). Keep it at 0 to count any detected methylation.*

#@markdown Different cutoffs for some genes (optional), e.g. `TAC1=38, HOXA7=39`:
ct_cutoff_per_gene = ""  #@param {type:"string"}
beta_cutoff_per_gene = ""  #@param {type:"string"}

#@markdown ### Cancer risk decision
#@markdown A sample is **Potentially cancer** when at least this many genes are methylated:
min_methylated_genes = 1  #@param {type:"integer"}
#@markdown Genes that count (comma separated; empty = all except the reference gene):
panel_genes = ""  #@param {type:"string"}

#@markdown ### Controls (leave empty to detect them from the sample names)
reference_gene = "B ACTIN"  #@param {type:"string"}
no_template_controls = ""  #@param {type:"string"}
positive_controls = ""  #@param {type:"string"}

#@markdown ### Advanced quality settings
reference_ct_max = 40  #@param {type:"number"}
reference_ct_low_input = 35  #@param {type:"number"}
min_cq_confidence = 0.5  #@param {type:"number"}
min_curve_height = 0.2  #@param {type:"number"}

#@markdown ---
show_r_code = False  #@param {type:"boolean"}
download_files = True  #@param {type:"boolean"}

import json, os, re, subprocess
from IPython.display import display, HTML

EDS_DIR = "/content/eds"
OUT_DIR = "/content/results"

def _r_str(s):
    return json.dumps(s)  # a JSON string is also a valid R string

def _names(text):
    return [t.strip() for t in text.split(",") if t.strip()]

def _exact_regex(names):
    return "^(" + "|".join(re.escape(n) for n in names) + ")$"

def _per_gene(default, text):
    parts = ["%s" % default]
    for item in _names(text):
        if "=" not in item:
            raise SystemExit("Per-gene cutoffs must look like GENE=value, got: " + item)
        gene, value = item.rsplit("=", 1)
        parts.append("`%s` = %s" % (gene.strip().replace("`", ""), float(value)))
    return "c(" + ", ".join(parts) + ")"

def build_r_code(eds_dir, out_dir):
    ref = reference_gene.strip()
    args = {
        "reference": _r_str(_exact_regex([ref])) if ref else "NULL",
        "ct_cutoff": _per_gene(ct_cutoff, ct_cutoff_per_gene),
        "beta_cutoff": _per_gene(beta_cutoff, beta_cutoff_per_gene),
        "ref_ct_max": reference_ct_max,
        "ref_ct_warn": reference_ct_low_input,
        "min_cq_conf": min_cq_confidence,
        "min_plateau": min_curve_height,
        "min_methylated_genes": int(min_methylated_genes),
        "panel": "c(%s)" % ", ".join(_r_str(g) for g in _names(panel_genes)) if _names(panel_genes) else "NULL",
    }
    if _names(no_template_controls):
        args["ntc"] = _r_str(_exact_regex(_names(no_template_controls)))
    if _names(positive_controls):
        args["positive"] = _r_str(_exact_regex(_names(positive_controls)))
    arg_text = ",\n".join("  %s = %s" % kv for kv in args.items())
    return f"""library(autoqmsp)
files <- list.files({_r_str(eds_dir)}, pattern = "\\\\.eds$", full.names = TRUE, ignore.case = TRUE)
runs <- read_eds_files(files, names = tools::file_path_sans_ext(basename(files)))
res <- analyze_qmsp(
  runs,
{arg_text}
)
dir.create({_r_str(out_dir)}, showWarnings = FALSE)
write_report(res, file.path({_r_str(out_dir)}, "qmsp_report.html"))
export_results(res, file.path({_r_str(out_dir)}, "qmsp_results.xlsx"))
print(res)
"""

def run_analysis(eds_dir=EDS_DIR, out_dir=OUT_DIR):
    if not os.path.isdir(eds_dir) or not any(f.lower().endswith(".eds") for f in os.listdir(eds_dir)):
        raise SystemExit("No .eds files found. Run step ② first.")
    code = build_r_code(eds_dir, out_dir)
    script = os.path.join(out_dir + "_script.R")
    with open(script, "w") as f:
        f.write(code)
    p = subprocess.run(["Rscript", script], capture_output=True, text=True)
    if p.returncode != 0:
        print(p.stderr[-3000:])
        raise SystemExit("The analysis failed - see the message above.")
    return code

_code = run_analysis()
if show_r_code:
    print(_code)
with open(os.path.join(OUT_DIR, "qmsp_report.html")) as f:
    display(HTML(f.read()))
if download_files:
    from google.colab import files
    files.download(os.path.join(OUT_DIR, "qmsp_report.html"))
    files.download(os.path.join(OUT_DIR, "qmsp_results.xlsx"))
'''

INTRO = """
# qMSP analysis: upload .eds files and get a report

1. Click ▶ on **step ①** to install the tool. You only need to do this once per session.
2. Click ▶ on **step ②** and choose your QuantStudio `.eds` files.
3. Adjust the settings in **step ③** if you need to, then click ▶. The report appears below it, and the report (HTML) and the Excel file download.

To change a setting, edit it and click ▶ on step ③ again. To see the code behind a step, double-click the step.

*For research use only. The results are not a diagnosis.*
"""

nb = {
    "cells": [md(INTRO), form(INSTALL), form(UPLOAD), form(ANALYSE)],
    "metadata": {"colab": {"provenance": [], "toc_visible": False},
                 "kernelspec": {"name": "python3", "display_name": "Python 3"},
                 "language_info": {"name": "python"}},
    "nbformat": 4, "nbformat_minor": 0,
}
here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, "autoqmsp_colab.ipynb"), "w") as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)
