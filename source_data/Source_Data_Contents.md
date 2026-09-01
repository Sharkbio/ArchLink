# Source Data Contents

This note is intended as an editor-facing guide to `Source_Data.xlsx`. It explains, sheet by sheet, which manuscript figure panel each table supports and what kind of source data are included.

For upload, the same sheets are also separated into one workbook per relevant display item under `source_data_by_figure`, following Nature Biotechnology's figure-specific source-data guidance. `Source_Data.xlsx` remains the combined master workbook.

`Source_Data.xlsx` is organized around the quantitative figure panels. Panels that are conceptual, schematic, externally rendered, or otherwise non-tabular are documented in the workbook sheets `Panel_Map` and `NonTabular_Notes` rather than being forced into artificial numeric tables.

Two scope notes should still be kept in mind:

- `Fig4c_marker_per_bin` now provides full bin-level paired SCG and duplicated-marker rows, but the published panel itself still visualizes dataset-level net changes.
- `Fig5c_BGC_paired` and `Fig5c_BGC_per_bin` now provide bin-resolved paired BGC support tables. The current revised source set is traceable at the `BinID` plus within-bin pair-index level; it does not expose native antiSMASH region identifiers as explicit columns in the exported summary sheets.

## Workbook Navigation Sheets

`README`
: Workbook-level scope statement. Explains that the file covers quantitative figure panels and records exclusions for non-tabular panels.

`Panel_Map`
: Crosswalk between manuscript figure panels and their plotting scripts, rendered subfigure files, and primary raw-data sources. This is a provenance sheet rather than a plotted-data sheet.

## Figure 2

`Fig2a_ARI_ACC`
: Supports Fig. 2a. Long-format table of ARI and ACC values by dataset and tool.

`Fig2b_CAMI_bins`
: Supports Fig. 2b. Recovered bin counts across completeness thresholds for the CAMI II datasets.

`Fig2c_real_MAGs`
: Supports Fig. 2c. Recovered MAG counts across completeness thresholds for the real-world datasets.

## Figure 3

`Fig3a_precision_TP`
: Supports Fig. 3a. Precision and true-positive contig-join counts for ArchLink ablations and MetaCarvel comparison modes.

`Fig3b_longread_main`
: Supports Fig. 3b. Primary long-read validation metrics in the Waste Water dataset: validation rate, depth consistency, structural purity, and anti-chimeric safety.

`Fig3c_N50_raw`
: Supports Fig. 3c. Per-bin before/after N50 values and the resulting N50 gains used to draw the distribution panel.

`Fig3c_N50_summary`
: Supports Fig. 3c. Dataset-level summary statistics and effect-size values used for the accompanying summary display.

## Figure 4

`Fig4a_tax_rank_counts`
: Supports Fig. 4a. Rank-level GTDB count table summarizing shared and workflow-specific classifications across COMEBin, ArchLink bin-only, and ArchLink scaffolded outputs. This sheet supports the rank-summary/alluvial-style panel; it is not a MAG-by-MAG transition table.

`Fig4a_unique_taxa`
: Supports Fig. 4a. Lists taxa that were unique to individual workflows at each rank and provides the specific labels underlying the rank-summary panel.

`Fig4b_align_frac_raw`
: Supports Fig. 4b. Raw paired entries for alignable-reference-fraction comparisons, including ANI, aligned fragments, total fragments, and aligned fraction.

`Fig4b_align_frac_sum`
: Supports Fig. 4b. Paired summary of ANI changes and aligned-fragment changes before and after scaffolding.

`Fig4c_marker_summary`
: Supports Fig. 4c. Dataset-level paired summary of SCG gains and duplicated-marker reductions after scaffolding. This sheet matches the net-change values displayed in the panel. Upstream source folder: `E:\binning_paper\graph\SCG`.

`Fig4c_marker_per_bin`
: Supports Fig. 4c. Bin-level paired SCG and duplicated-marker table derived from the pre/post GTDB-Tk marker summaries. This sheet provides the detailed support table for which bins changed and which marker sets were gained or lost. Upstream source folder: `E:\binning_paper\graph\SCG`.

`Fig4d_DRAM_raw`
: Supports Fig. 4d. Raw DRAM module-completeness values for the selected 14 paired MAGs analyzed before and after scaffolding.

`Fig4d_DRAM_summary`
: Supports Fig. 4d. Dataset-level summary of DRAM completeness changes within the selected 14-MAG paired subset.

## Figure 5

`Fig5a_BGC_counts`
: Supports Fig. 5a. Total BGC counts by dataset, tool, and ArchLink recovery category.

`Fig5b_BGC_lengths`
: Supports Fig. 5b. Per-call BGC total length and core-length values used for the violin plots, including `BinID`, state and boundary-completeness status.

`Fig5b_violin_input`
: Supports Fig. 5b. Helper input table used for the violin-length plotting workflow.

`Fig5c_BGC_class_resolved`
: Supports Fig. 5c directly. Class-resolved scaffold-associated BGC summary by dataset, including median positive total-length gain, the number of BGCs becoming boundary-complete after scaffolding, and total class counts. Upstream source table: `F:\Binlink-main\BinLINK\BinLINK\bgc_class_detail.csv`; rendered figure script: `E:\binning_paper\graph\fig5c_bgc_class\BGC_stats_classify_graph.py`.

`Fig5c_BGC_paired`
: Supports Fig. 5c and also underlies the paired summaries discussed for Fig. 5b. Contains bin-resolved paired before/after BGC length rows, within-bin pair indices and linked boundary-completeness counts. Upstream source folder: `E:\binning_paper\graph\bgc_length`.

`Fig5c_BGC_per_bin`
: Supports Fig. 5c. Per-bin counts of boundary-complete BGCs before and after scaffolding. Upstream source folder: `E:\binning_paper\graph\bgc_length`.

`Fig5_summary`
: Figure 5 support sheet. Dataset-level summary of BGC counts, median lengths, median core lengths, and boundary-completeness changes.

## Figure 6

`Fig6a_scaffold_summary`
: Supports Fig. 6a. Scaffold component lengths, contig orientations, junction identity values and the before/after annotated lengths of the NAGGN and terpene/RiPP-like BGC regions.

`Fig6d_PAE_matrix`
: Supports Fig. 6d. Full predicted aligned error matrix used for the protein-structure panel.

`Fig6e_depth_profile`
: Supports Fig. 6e. Per-position short-read depth and clipped-read counts across the `connect_00055` fusion-gene region.

## Supplementary Figures

`SuppFig1c_gap_support`
: Supports Supplementary Fig. 11b. Gap-level case-study validation values including spanning long-read counts, depth ratios and junction identity values.

## Non-Tabular Panels

`NonTabular_Notes`
: Records panels whose underlying source is not most naturally represented as a flat quantitative table. This includes conceptual schematics, styled tree renders, gene-cluster cartoons, molecular structure renders, and software-generated comparative-genomics visualizations. The current notes cover Fig. 1a-d, Fig. 6a-c and Supplementary Fig. 10, 11a and 12, with source artifacts cross-referenced to `Panel_Map`.
