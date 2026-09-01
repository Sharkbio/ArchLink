# Figure-generation entry points

The scripts in this directory regenerate the data-driven manuscript panels from
the versioned `source_data/Source_Data.xlsx` workbook.  They use a fixed random
seed (`20260831`), explicit sheet names, fixed dataset/tool ordering, a shared
color map, and a default output resolution of 300 dpi.  Every run writes
`figure_build_manifest.json` with the workbook SHA-256, Git commit, software
versions, parameters and output paths.

## One-command build

From the repository root:

```bash
python figures/make_all_figures.py \
  --source-data source_data/Source_Data.xlsx \
  --output figures/output \
  --formats pdf
```

Install the pinned figure-only dependencies with
`pip install -r figures/requirements.txt`, or use the repository conda
environment after its figure dependencies have been installed.

## Panel map

| Figure | Source Data sheet(s) | Script | Default output |
|---|---|---|---|
| Fig. 2a | `Fig2a_ARI_ACC` | `fig2_binning_benchmark.py` | `Fig2a.pdf` |
| Fig. 2b | `Fig2b_CAMI_bins` | `fig2_binning_benchmark.py` | `Fig2b.pdf` |
| Fig. 2c | `Fig2c_real_MAGs` | `fig2_binning_benchmark.py` | `Fig2c.pdf` |
| Fig. 3a | `Fig3a_precision_TP` | `fig3_link_validation.py` | `Fig3a.pdf` |
| Fig. 3b | `Fig3b_longread_main` | `fig3_link_validation.py` | `Fig3b.pdf` |
| Fig. 3c | `Fig3c_N50_raw`, `Fig3c_N50_summary` | `fig3_link_validation.py` | `Fig3c.pdf` |
| Fig. 4a | `Fig4a_tax_rank_counts` | `fig4_taxonomy_metabolism.py` | `Fig4a.pdf` |
| Fig. 4b | `Fig4b_align_frac_raw` | `fig4_taxonomy_metabolism.py` | `Fig4b.pdf` |
| Fig. 4c | `Fig4c_marker_summary` | `fig4_taxonomy_metabolism.py` | `Fig4c.pdf` |
| Fig. 4d | `Fig4d_DRAM_raw` | `fig4_taxonomy_metabolism.py` | `Fig4d.pdf` |
| Fig. 5a | `Fig5a_BGC_counts` | `fig5_bgc_analysis.py` | `Fig5a.pdf` |
| Fig. 5b | `Fig5b_BGC_lengths` | `fig5_bgc_analysis.py` | `Fig5b.pdf` |
| Fig. 5c | `Fig5c_BGC_class_resolved` | `fig5_bgc_analysis.py` | `Fig5c.pdf` |
| Fig. 6a | `Fig6a_scaffold_summary` | `fig6_quantitative.py` | `Fig6a.pdf` |
| Fig. 6d | `Fig6d_PAE_matrix` | `fig6_quantitative.py` | `Fig6d.pdf` |
| Fig. 6e | `Fig6e_depth_profile` | `fig6_quantitative.py` | `Fig6e.pdf` |

## Statistical and display conventions

- Fig. 2b filters the source table at precision >= 90%; Fig. 2c uses precision
  >= 95% and completeness thresholds >= 50%.
- Fig. 3a displays the across-dataset mean and standard deviation for precision
  and true-positive joins.
- Fig. 3c displays every bin-level N50 gain, source-workbook mean/95% CI values,
  and the source-workbook approximate effect size on a symmetric log axis.
- Fig. 4d aggregates post-scaffolding numeric DRAM module completeness into the
  documented ETC, carbon-fixation and central-carbon categories.
- Fig. 5b uses the observation-level total BGC lengths on a logarithmic axis;
  the violin median is drawn directly from those observations.

## Non-tabular and assembled panels

Fig. 1 is a manually assembled conceptual schematic and is not claimed to be
fully script generated.  Fig. 4a's workbook contains rank-level counts rather
than the MAG-level transition table required to recreate the original alluvial
geometry; the supplied script therefore regenerates the quantitative rank-count
summary, not the Illustrator layout.  Fig. 6b is generated from antiSMASH/GenBank
context files and Fig. 6c from an AlphaFold/PyMOL structure render; their final
composition remains an external/manual assembly step.  The quantitative Fig. 6
panels (a, d and e) have explicit scripts here.
