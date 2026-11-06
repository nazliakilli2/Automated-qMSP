# autoqmsp: automated qMSP analysis

Upload QuantStudio / Applied Biosystems run files (`.eds`) and get a report
that says, for every sample, whether it is **Potential cancer**,
**Not determined** (repeat) or **Low risk**. You don't need to export
anything manually or work in Excel.

**The rule:** the Ct is read where the amplification curve crosses
**ΔRn = 10,000**, and a gene is **methylated** when its **Ct is 40 or less**.
Instead of 40 you can use a cutoff per gene, found from samples whose
diagnosis you know with the [cutoff notebook](#finding-a-ct-cutoff-per-gene).

> For research use only. The results are not a diagnosis.

## Option A: Google Colab (free, nothing to install, no code shown)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/nazliakilli2/Automated-qMSP/blob/claude/gracious-ptolemy-t3d3d8/colab/autoqmsp_colab.ipynb)

The notebook has four steps. Each one is a button with a form, and the code
stays hidden:

1. **Install**: click ▶ and wait for *Ready*.
2. **Upload**: click ▶, then *Choose Files* and pick your `.eds` files. The
   amplification curves of every run and gene are drawn, with the threshold
   (green, ΔRn 10,000) and the Ct cutoff (red, 40). Curves of samples that
   cross the threshold by Ct 40 are red. Hover over a curve to see the sample
   and its Ct.
3. **Look at each gene**: pick a run and gene to see its curves, a table of
   every well's Ct, and control checks (for example "⚠ No-template control
   amplified").
4. **Report**: click ▶. The report and a chart of the Ct per gene appear in
   the notebook, and the report (HTML) and Excel file download. Tick
   *use_cutoffs_file* to upload a `qmsp_cutoffs.csv` from the cutoff notebook;
   its genes then use their own cutoff instead of 40.

To see the code behind a step, double-click the step. Tick *show_r_code* to
print the R code of the analysis.

## Finding a Ct cutoff per gene

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/nazliakilli2/Automated-qMSP/blob/claude/gracious-ptolemy-t3d3d8/colab/autoqmsp_cutoffs_colab.ipynb)

A second notebook uses samples whose diagnosis you know to choose, for every
gene, the Ct cutoff that best separates cancer from control samples.

1. **Install**, as above.
2. **Upload** the `.eds` files. Samples with the same name in different runs
   are pooled.
3. **Groups**: for every sample choose *Cancer*, *Control* (no cancer),
   *Positive control*, *NK* (water) or *Leave out*. NK and positive controls
   are pre-filled from their names. A labels file (`sample`, `group`) saved
   earlier can fill the list.
4. **Find the cutoffs**: you get a suggested cutoff per gene with its
   sensitivity, specificity and AUC. Pick a gene to see every sample's Ct by
   group and the ROC curve. Move the slider to try another cutoff: the
   numbers, the graphs and the result for all genes together update.
5. **Save**: downloads `qmsp_cutoffs.csv` (for the main notebook), an Excel
   file with all the numbers, and your labels.

How the cutoff is chosen:

- Each sample gets one Ct per gene at ΔRn 10,000: the mean of its wells, or
  no Ct if fewer than half of its wells reach the threshold. Samples whose
  reference gene has no Ct ≤ 40 (failed DNA) are left out.
- A sample counts as methylated when its Ct is at or below the cutoff. For
  every possible cutoff, **sensitivity** = cancers at or below it and
  **specificity** = controls above it or with no Ct.
- **Best balance (Youden)** picks the highest sensitivity + specificity. You
  can instead require a minimum specificity (few false cancer results) or
  sensitivity (few missed cancers). The cutoff sits halfway between the two
  nearest sample Cts.
- A gene gets a warning with fewer than 10 cancers or 10 controls, when the
  positive control is above the cutoff, or when NK is at or below it. Check
  the cutoffs on new samples that were not used to choose them.

In R:

```r
runs   <- read_eds_files("training/eds/files")
groups <- read.csv("labels.csv")          # columns: sample, group
cut    <- find_cutoffs(runs, groups)      # or method = "specificity", target = 0.9
cut$cutoffs
plot_cutoffs(cut)
write_cutoffs(cut, "qmsp_cutoffs.csv")
res <- analyze_qmsp(new_runs, ct_cutoff = read_cutoffs("qmsp_cutoffs.csv"))
```

## Option B: the app (point and click)

```r
install.packages(c("remotes", "shiny", "writexl"))
remotes::install_github("nazliakilli2/Automated-qMSP")
autoqmsp::run_app()
```

| Tab | What you do |
|---|---|
| **1. Upload** | Choose the `.eds` files. The runs, samples and genes found are listed. |
| **2. Settings** | Shows the fixed rule (ΔRn 10,000, Ct ≤ 40). You choose how many methylated genes make a sample *Potential cancer*, which genes count, the reference gene, and the control samples (detected automatically, editable). |
| **3. Report** | The verdict for every sample, the Ct of every gene, and the control checks. Download the report (HTML, printable to PDF) or Excel. The R code is shown only if you tick *Show the R code*. |
| **4. Details** | Heatmap, all results, amplification curves, plate layout. |

To give the whole lab one link, deploy the `inst/shiny` folder to
[shinyapps.io](https://www.shinyapps.io) or
[Posit Connect Cloud](https://connect.posit.cloud). Both have free tiers.

## How the results are calculated

**1. Ct.** For every well, the Ct is read from its amplification curve at
**ΔRn = 10,000** (all genes, including the reference gene). The Ct is the
cycle where the curve crosses 10,000 for the last time and stays above it,
interpolated between cycles. A curve that never reaches 10,000 has no Ct. The
instrument's own Ct is kept in the Excel file for comparison.

**2. Methylated or not.** A gene is **methylated** when its Ct is **40 or
less**. Otherwise it is unmethylated, unless the result cannot be trusted,
in which case it is **not determined**:

- the sample's reference gene (e.g. β-actin) has no Ct of 40 or less (DNA
  failed), or
- the gene's no-template (water) control amplified, for a methylated result,
  or
- the gene's positive control did not amplify, for an unmethylated result,
  or
- fewer than half of the replicates agree.

**3. Verdict per sample.**

| Verdict | Rule |
|---|---|
| **Potential cancer** | at least *N* genes of the panel are methylated (default *N* = 1) |
| **Not determined** | fewer than *N*, but genes that could not be determined could change that; repeat the sample |
| **Low risk** | fewer than *N* genes methylated, even counting the ones not determined |

ΔCt, PMR and beta are still calculated and listed in the Excel file for
information. They do not change the result.

### Settings

| Setting | Value |
|---|---|
| Threshold (ΔRn) | 10,000 (fixed) |
| Ct cutoff | 40, or per gene from a cutoffs file |
| Reference gene Ct max / low-input warning | 40 / 35 |
| Methylated genes for *Potential cancer* | 1 |

## Option C: R code

```r
library(autoqmsp)
runs <- read_eds_files("folder/with/eds/files")
res  <- analyze_qmsp(runs, reference = "B ACTIN")   # ΔRn 10,000, Ct <= 40
res$report                           # verdict per sample
res$results                          # Ct and call per gene
write_report(res, "qmsp_report.html")
export_results(res, "qmsp_results.xlsx")
plot_methylation(res)
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
