# Task list
--------------------------------------------------------------------------------
## Important biology
- [X] Principle component plot
- [X] Timecourse scatterplot (modified volcano)
- [X] Facet categories from REACTOME GSEA plot
- [X] GSEA REACTOME heatmap with NES at 0, 12, 24 hours p.i.
- [X] correlate viral load data with gene expression data
   - [X] dot plot of viral replication markers
   - [X] add to makefile 

### Comments from NAMSED meeting 5-23-18
- [ ] Plot PC1,2,3 vs. time, viral titer to show association of PCs with other experimental variables
  - [ ] extract eigenvalues to show major genes driving variance
- [ ] Annotate GO pathways with parent-child relationships
- [ ] Broad overview figure
  - schematic showing experimental design
  - Timecourse scatterplot
  - tables with number of DE genes
  - Principle component analysis
- [ ] Broad activation of antiviral response figure
  - GSEA REACTOME heatmap with NES at 0, 12, 24 hours p.i.
  - correlate viral load data with gene expression data

## Coding and database maintenance
- [ ] write make for kallisto alignment and DESeq2 count export
- [ ] write script to check host machine for installed R dependencies and install any missing packages
- [X] Test function of 'test = "Wald"' in DESeq2::results()
- [ ] Create archive of raw data prerequisites

## Input from collaborators
- [ ] Edit list of representative pathways for GSEA figure
  - GSEA by time ('img/gsea-reactome-heatmap.png')
  - GSEA by viral titer ('results/GSEA/GSEA_titer-corr_REACTOME.csv')

