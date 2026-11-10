# Comparing Gene Expression in Childhood B-Cell and T-Cell Leukemia

## Research Question

Which genes show different expression levels between B-cell and
T-cell acute lymphoblastic leukemia (ALL) in children?

## Dataset

This project uses the public NCBI GEO dataset [GSE13425](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13425).
I will check the sample labels to identify the two groups.

## Analysis Plan

I will import the data into R using GEOquery, check data quality,
and use limma to compare gene expression between the two groups.
I will use a heatmap to show the genes with the strongest evidence
of a difference.
