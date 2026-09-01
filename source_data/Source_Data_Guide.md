# Source Data Guide for ArchLink

This guide is the editor-facing companion to `source_data/Source_Data.xlsx`. It explains what the source-data file covers, how it differs from Zenodo and Supplementary Data, and how to keep manuscript wording consistent with the actual deposited materials.

## What the Source Data file is

For Nature Portfolio-style handling, `Source_Data.xlsx` is the workbook that contains the exact plotted or tabulated numerical values underlying the manuscript's quantitative figure panels.

It is not the same thing as:

- raw sequencing files such as `FASTQ`, `BAM`, or long-read alignments
- the broader processed-results archive in Zenodo
- Supplementary Information text
- the rendered figure files themselves such as `PDF`, `PNG`, or `SVG`

Operationally:

- bar plots should map to the values behind each bar
- scatter plots should map to the plotted coordinates
- violin or box plots should map to the observation-level values used to draw each distribution
- heatmaps should map to the underlying matrix or long-format values

## Relationship to Zenodo and Supplementary Data

These three deliverables serve different roles and should be kept distinct:

- `Source_Data.xlsx` supports the paper's quantitative figure panels
- `supplementary_data\` provides extended paired tables, dataset summaries and bin-resolved support files
- the reduced Zenodo package at `F:\zenodo_package_reduced` is the broader processed-results archive for reviewer or post-publication access

Because the reduced Zenodo package was intentionally slimmed down, the source-data workbook should still be distributed as its own submission artifact even if many upstream tables also exist inside the analysis workspace.

## Current ArchLink Source Data package

The current workbook already exists at:

- `source_data/Source_Data.xlsx`

Nature Biotechnology's figure-formatting guidance asks for statistics source data as one Excel file per relevant figure. Submission-ready figure-specific copies are therefore also provided under:

- the journal submission's figure-specific Source Data exports

The split files cover Fig. 2, Fig. 3, Fig. 4, Fig. 5, Fig. 6 and the currently tabulated Supplementary Fig. 11b validation table. The combined workbook remains the master cross-figure source-data record.

Its sheet-by-sheet contents are documented in:

- `source_data/Source_Data_Contents.md`

That file should be treated as the detailed sheet crosswalk, while the present guide should be treated as the higher-level submission note.

## What the workbook currently covers

### Figure 2

- `Fig2a_ARI_ACC`
- `Fig2b_CAMI_bins`
- `Fig2c_real_MAGs`

These sheets capture the assembly-level and MAG-recovery benchmarking values used in the main comparative plots.

### Figure 3

- `Fig3a_precision_TP`
- `Fig3b_longread_main`
- `Fig3c_N50_raw`
- `Fig3c_N50_summary`

These sheets cover scaffold-join precision, primary long-read validation metrics and paired N50 changes.

### Figure 4

- `Fig4a_tax_rank_counts`
- `Fig4a_unique_taxa`
- `Fig4b_align_frac_raw`
- `Fig4b_align_frac_sum`
- `Fig4c_marker_summary`
- `Fig4c_marker_per_bin`
- `Fig4d_DRAM_raw`
- `Fig4d_DRAM_summary`

Important note: the plotted Fig. 4c panel is dataset-level, but the workbook now also includes the supporting bin-level paired marker table so that the net bar values remain traceable to specific bins.

### Figure 5

- `Fig5a_BGC_counts`
- `Fig5b_BGC_lengths`
- `Fig5b_violin_input`
- `Fig5c_BGC_class_resolved`
- `Fig5c_BGC_paired`
- `Fig5c_BGC_per_bin`
- `Fig5_summary`

These sheets collectively cover end-to-end BGC counts, per-call length distributions, the class-resolved Fig. 5c panel table, paired before/after BGC summaries and per-bin boundary-completeness support tables.

### Figure 6

- `Fig6a_scaffold_summary`
- `Fig6d_PAE_matrix`
- `Fig6e_depth_profile`

### Supplementary Figures

- `SuppFig1c_gap_support`

## Panels that are intentionally not flattened into numeric sheets

Some panels are primarily schematic, structural, tree-rendered, gene-cluster-cartoon or software-generated visual outputs. For these, the workbook relies on:

- `Panel_Map`
- `NonTabular_Notes`

inside `Source_Data.xlsx`, plus the external explanation in `Source_Data_Contents.md`.

This is a normal and defensible approach as long as the manuscript does not claim that every non-tabular panel has a standalone numeric source table.

## Practical submission rule

If an editor or reviewer asks, "What exact numbers were used to draw this quantitative panel?", the answer should be findable in `Source_Data.xlsx`.

If they ask, "Where are the raw reads, full MAG FASTA files, antiSMASH outputs, or other processed archives?", the answer should point to SRA, Zenodo, or the supplementary-data package instead.

## Recommended manuscript wording

If the workbook is included in the submission package, the manuscript can safely retain wording such as:

`Source data underlying quantitative figure panels are provided with the paper.`

If the workbook is temporarily omitted from a transfer package, that sentence should be revised before submission.

## Final pre-submission check

Before uploading the final package, confirm all four items below:

1. `Source_Data.xlsx` is present in the submission materials.
2. `Source_Data_Contents.md` still matches the actual workbook sheets.
3. The manuscript figure legends and `Data availability` statement do not promise a source-data artifact that is missing from the upload.
4. The current Zenodo package wording does not imply that the source-data workbook is embedded there unless it truly has been added.
