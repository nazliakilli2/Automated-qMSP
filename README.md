# autoqmsp: automated qMSP analysis

Upload QuantStudio / Applied Biosystems run files (`.eds`) and get a report
that says, for every sample, whether it is **Potentially cancer**,
**Inconclusive** (repeat) or **Not risky**, together with a methylation level
(**beta**, 0–1) for every gene. You don't need to export anything manually or
work in Excel.

> For research use only. The results are not a diagnosis.

## Option A: Google Colab (free, nothing to install, no code shown)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/nazliakilli2/Automated-qMSP/blob/claude/gracious-ptolemy-t3d3d8/colab/autoqmsp_colab.ipynb)

The notebook has three steps. Each one is a button with a form, and the code
stays hidden:

1. **Install**: click ▶ and wait for *Ready*.
2. **Upload**: click ▶, then *Choose Files* and pick your `.eds` files.
3. **Settings and report**: adjust the sliders and fields if needed, then click ▶.
   The report appears in the notebook, and the report (HTML) and Excel file
   download.

To see the code behind a step, double-click the step. Tick *show_r_code* to
print the R code of the analysis.

## Option B: the app (point and click)

```r
install.packages(c("remotes", "shiny", "writexl"))
remotes::install_github("nazliakilli2/Automated-qMSP")
autoqmsp::run_app()
```

| Tab | What you do |
|---|---|
| **1. Upload** | Choose the `.eds` files. The runs, samples and genes found are listed. |
| **2. Settings** | Ct and beta cutoffs (also per gene), how many methylated genes make a sample *Potentially cancer*, which genes count, the reference gene, and the control samples (detected automatically, editable). Advanced quality settings are hidden behind a checkbox. |
| **3. Report** | The verdict for every sample, beta values, control checks and wells to check. Download the report (HTML, printable to PDF) or Excel. The R code is shown only if you tick *Show the R code*. |
| **4. Details** | Heatmap, all results, amplification curves, plate layout. |

To give the whole lab one link, deploy the `inst/shiny` folder to
[shinyapps.io](https://www.shinyapps.io) or
[Posit Connect Cloud](https://connect.posit.cloud). Both have free tiers.

## How the results are calculated

**1. Quality control of each well.** An amplified well is sent to *Review*
instead of being trusted when:

- the instrument's Cq confidence is low, or
- the instrument says there was no amplification or it was inconclusive, or
- its curve is much lower than the positive control's (a flat "drift" curve
  rather than a real amplification).

**2. Controls.** For each run and gene:

- No-template (water) controls must stay negative.
- Positive controls must amplify.
- The reference gene (e.g. β-actin) must amplify in every sample. A failure
  makes the sample *Invalid*; a high Ct flags *low DNA input*.

**3. Methylation level (beta, 0–1).**

- beta = PMR / 100, capped at 1.
- PMR (percentage of methylated reference) = 100 × 2^-ΔCt(sample) /
  2^-ΔCt(positive control), where ΔCt = Ct(gene) − Ct(reference gene).
- 0 means no methylation detected. 1 means as methylated as the fully
  methylated positive control.
- In runs without a reference gene, beta = 2^-(Ct(sample) − Ct(positive control)).

**4. Methylated or not.** A gene is **Methylated** when all of these hold:

- Ct ≤ the Ct cutoff
- the wells pass QC
- beta ≥ the beta cutoff

You can set both cutoffs per gene.

**5. Verdict per sample.**

| Verdict | Rule |
|---|---|
| **Potentially cancer** | at least *N* of the panel genes are methylated (default *N* = 1) |
| **Inconclusive** | fewer than *N*, but genes that need review or failed could change that; repeat the sample |
| **Not risky** | fewer than *N* genes methylated, even counting the uncertain ones |

### Default settings (all adjustable)

| Setting | Default |
|---|---|
| Ct cutoff | 40 |
| Beta cutoff | 0 (any detected methylation counts) |
| Methylated genes for *Potentially cancer* | 1 |
| Reference gene invalid / low-input Ct | > 40 / > 35 |
| Minimum Cq confidence | 0.5 |
| Minimum curve height | 20 % of positive control |

## Option C: R code

```r
library(autoqmsp)
runs <- read_eds_files("folder/with/eds/files")
res  <- analyze_qmsp(runs, reference = "B ACTIN", ct_cutoff = c(40, TAC1 = 38),
                     beta_cutoff = 0.05, min_methylated_genes = 2)
res$report                           # verdict per sample
res$results                          # Ct, ΔCt, PMR, beta, call per gene
write_report(res, "qmsp_report.html")
export_results(res, "qmsp_results.xlsx")
plot_methylation(res, "beta")
plot_amplification(res, run = "my run")
```

## Supported files

- `.eds` files from QuantStudio 3/5/7 (Design & Analysis / QuantStudio
  software v1.x), analysed and saved after the run.
- `.edt` template files hold no results and cannot be analysed.
- A sample's gene wells and reference-gene well are matched by **sample
  name**, so use exactly the same name in the plate setup.

Raw instrument files are ignored by git (`.gitignore`) so that patient data is
not committed by accident.
