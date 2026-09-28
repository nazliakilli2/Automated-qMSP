# autoqmsp: automated qMSP analysis

`autoqmsp` is an R package that turns QuantStudio / Applied Biosystems run files
(`.eds`) into quantitative methylation-specific PCR (qMSP) results. It needs no
manual export or Excel work. It reads the `.eds` file directly and gives you:

- Ct values for every well, with *Undetermined* handled correctly
- **Quality control**:
  - no-template controls (NTC) and positive controls, per run and gene
  - reference gene (e.g. β-actin) checks for each sample
  - automatic flags for suspicious amplification: low instrument Cq confidence, the instrument itself saying "no amplification", or a low, linear "drift" curve compared with the positive control
- **Methylation values**: ΔCt, ratio (2^-ΔCt), PMR (percentage of methylated reference)
- **A call** per sample and gene: `Methylated`, `Unmethylated`, `Review` or `Invalid`
- Plots: amplification curves, plate map, methylation heatmap
- Excel export (summary, PMR, full results, controls, wells, settings)

You can use it three ways, all free.

## 1. Google Colab (no installation)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/nazliakilli2/Automated-qMSP/blob/claude/gracious-ptolemy-t3d3d8/colab/autoqmsp_colab.ipynb)

1. Open the notebook with the badge above. It runs R in your browser.
2. Run the first cell to install the package.
3. Click the folder icon on the left and drag in your `.eds` files.
4. Run the remaining cells, then download `qmsp_results.xlsx`.

## 2. Point-and-click app (Shiny)

```r
install.packages(c("remotes", "shiny", "writexl"))
remotes::install_github("nazliakilli2/Automated-qMSP")
autoqmsp::run_app()
```

Upload the `.eds` files, adjust the settings in the sidebar, and download the
Excel report. To share the app with the whole lab through a link, deploy the
`inst/shiny` folder to [shinyapps.io](https://www.shinyapps.io) or
[Posit Connect Cloud](https://connect.posit.cloud). Both have free tiers.

## 3. R code

```r
library(autoqmsp)

runs <- read_eds_files(c("run1.eds", "run2.eds"))   # or a folder of .eds files
res  <- analyze_qmsp(
  runs,
  reference   = "B ACTIN",        # reference gene; NULL if the run has none
  ntc         = "NTC|dH2O|dH20|^NK",
  positive    = "H460|A549|HT29", # methylated positive control samples
  ct_cutoff   = 40                # methylated if gene Ct <= 40
)

res                          # short summary + control problems
results_wide(res, "call")    # sample x gene table
res$results                  # everything: Ct, ref Ct, ΔCt, ratio, PMR, notes
res$wells[res$wells$result == "Review", ]   # wells to look at by eye

plot_methylation(res)
plot_amplification(res, run = "run1")
plot_plate(res, run = "run1")
export_results(res, "qmsp_results.xlsx")
```

## How the calls are made

| Step | Rule (defaults, all adjustable) |
|---|---|
| Well role | Sample name matches `ntc` → NTC; matches `positive` → positive control; otherwise sample |
| Well result | **Negative** if Ct is undetermined or above `ct_cutoff` (40). **Review** if amplified but Cq confidence < `min_cq_conf` (0.5), the instrument says no / inconclusive amplification, or the final ΔRn is < `min_plateau` (20 %) of the positive controls' curves for that gene and run. Otherwise **Positive**. |
| NTC check | *Fail* if any NTC well is Positive, *Review* if any NTC well is Review, for that gene and run |
| Reference gene | *Failed* if Ct is undetermined or > `ref_ct_max` (40) → sample is **Invalid**. *Low input* if > `ref_ct_warn` (35). |
| Call | **Invalid** if the reference failed; **Unmethylated** if the gene is negative; **Review** if the gene well needs review or the NTC for that gene is not clean; otherwise **Methylated** |
| ΔCt | Ct(gene) − Ct(reference) |
| Ratio | 2^-ΔCt (0 for unmethylated samples) |
| PMR | 100 × ratio(sample) / mean ratio(positive controls), same gene and run |

Replicate wells of the same sample and gene (in the same run) are averaged.
The gene and reference wells of a sample are matched by **sample name**, so
use exactly the same name in the plate setup.

## Supported files

- `.eds` files from QuantStudio 3/5/7 (Design & Analysis / QuantStudio software
  v1.x), analysed and saved after the run.
- An unzipped `.eds` folder, or its `apldbio/sds/analysis_result.txt`.
- `.edt` template files hold no results and cannot be analysed.

Raw instrument files are ignored by git (`.gitignore`) so that patient data
is not committed by accident.
