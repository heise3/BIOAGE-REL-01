# Post hoc comparison of conventional indices and PhenoAgeAccel

The v1.0 contract was frozen after the primary results were known. The unchanged numerical program is in code/run.R. The frozen/ directory holds aggregate results and gate records from the single study run. The two tableS7 TSV files are the source of manuscript Supplementary Table S19, panels A–C.

To independently reproduce the comparison, first run the repository's reproduce.sh in an empty analysis workspace. From a separate, empty output directory, run:

    Rscript /absolute/path/to/repository/S7/code/run.R /absolute/path/to/analysis-workspace

The program reads the primary workflow's local NHANES inputs and intermediate files, writes its numerical TSVs into the current directory and checks the reference ICC and mortality estimates before the comparisons. Existing completed output files prevent an accidental second run there. Compare the TSVs against frozen/. A1_variables.tsv documents the source fields and laboratory methods. Missing or ineligible results retain their original status.
