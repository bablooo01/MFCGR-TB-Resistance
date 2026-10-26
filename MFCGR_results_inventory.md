# MFCGR Project — Results Inventory

## Project
Methylation-aware Fractal Chaos Game Representation for Characterizing
Drug-Resistant Mycobacterium tuberculosis Genomes

## Results inventory

- [PRESENT] `selected_120_isolates.csv` — 120 data rows
- [PRESENT] `methyltransferase_gene_variants_120.csv` — 504 data rows
- [PRESENT] `methyltransferase_variant_consequences.csv` — 504 data rows
- [PRESENT] `cohort_top_variant_validation.csv` — 240 data rows
- [PRESENT] `top_variant_fisher_tests.csv` — 6 data rows
- [PRESENT] `candidate_annotation_check.csv` — 8 data rows
- [PRESENT] `candidate_allele_group_summary.csv` — 8 data rows
- [PRESENT] `preliminary_candidate_gene_ranking.csv` — 2 data rows
- [PRESENT] `circos_candidate_data.csv` — 16 data rows
- [PRESENT] `candidate_gene_sequences.fasta`
- [PRESENT] `candidate_sequence_summary.csv` — 4 data rows

## Generated figures

- [PRESENT] `figures/CGR_Rv2024c.png`
- [PRESENT] `figures/CGR_Rv2756c.png`
- [PRESENT] `figures/CGR_Rv2761c.png`
- [PRESENT] `figures/CGR_Rv3263.png`
- [PRESENT] `figures/candidate_loci_circular_overview.png`

## Current exploratory findings

- 120 isolates were selected, with 30 in each resistance category.
- Four candidate loci were examined: Rv2024c, Rv2756c, Rv2761c and Rv3263.
- Two candidate positions were assessed using read-level evidence across the cohort.
- Fisher's exact tests and FDR-adjusted q-values were calculated for the two positions.
- Reference-sequence CGR plots and a circular candidate-locus overview were generated.

## Interpretation and limitations

- Association results are exploratory and do not establish causation.
- Population structure and lineage effects require further assessment.
- Coding consequences require independent validation before biological conclusions.
- Rv2024c is annotated as a hypothetical protein in the reference annotation examined.
- The current CGR plots use reference CDS sequences, not isolate-specific sequences.
- The circular overview is a Circos-style plot, not a full Circos track analysis.
- Illumina sequence variants do not directly measure DNA methylation.
- A validated methylation-aware analysis requires appropriate methylation measurements
  or a clearly defined and independently justified alternative.

## Pending before final submission

1. Independently validate variant calls and coding consequences.
2. Verify gene identities and functions against authoritative annotations and literature.
3. Assess lineage/population structure and robustness of resistance associations.
4. Obtain compatible methylation evidence or narrow the methylation claims.
5. Extend sequence analysis to isolate-specific sequences where appropriate.
6. Finalize gene ranking using transparent criteria and validated evidence.
7. Generate final figures and write the methods, results, discussion and limitations.
