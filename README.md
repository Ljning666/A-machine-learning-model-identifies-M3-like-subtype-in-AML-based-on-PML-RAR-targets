# A-machine-learning-model-identifies-M3-like-subtype-in-AML-based-on-PML-RAR-targets
787_GSEAscore.R conclude the code of machine learning model.
The input data is a list of genes ranked from highest to lowest FC for each AML patient, calculated by comparing their expression levels with the average expression levels in normal samples.
The gene set activated by PML/RARα or genes set repressed by PML/RARα was used as a gene set to perform GSEA.
The output data is the M3-LS' score for each AML patient. The M3-LS' scores for all patients were further normalized in the range from 0 to 1.
