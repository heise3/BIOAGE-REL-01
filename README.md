# BIOAGE-REL-01
Analysis code for **Reliable ranking, unreliable change: short-term retest variability of blood-chemistry PhenoAge in repeat NHANES blood draws**.

This release accompanies manuscript v0.3.3. It contains the frozen v1.0 screening contract, its implementation lock, the public input inventory with SHA-256 checksums, the primary and replication analysis, S1–S4 and S6 quantitative supplement scripts, selected aggregate results and the five figure renderers. S5 was a dated literature review. Its historical decision is supplied for the S1 entry check. The strict-subset and other cross-project results used in the figures are supplied as published aggregate inputs. This release covers REL-01 and does not reproduce the separate RCV project in full.

S7/ contains the post hoc comparison of conventional blood indices and PhenoAgeAccel with its original code, frozen contract, calibration and aggregate tables. See S7/README.md for its entry point.

## Reproduce
Use Python 3 (standard library) and R with haven, readr, dplyr, jsonlite, survey, survival, ggplot2, patchwork and ragg. See core/sessionInfo.txt for the archived analysis environment. Arial, Cairo graphics and a working TIFF device are required for figure export. Network access to CDC/NCHS and the pinned BioAge source repository is required on the first run.

In an empty working directory, run `bash /absolute/path/to/repository/reproduce.sh`. The scripts download NHANES III and NHANES 2001–2002 public-use files, validate the frozen core input hashes, generate local person-level intermediate files, and build aggregate outputs. If a frozen input changes, inspect the mismatch before using new results. For figure-only reproduction after S1, run `bash /absolute/path/to/repository/render.sh work/S/s1_results.rds figures-build`. The two individual-level plotting inputs are reconstructed on the reader's machine. All five figures are written under figures-build/figures/.

Primary NHANES III analyses used 1,801 participants. The replication included 407 paired participants. The estimates concern short-term retest variation in these analytic samples. The repository archives an existing analysis and does not establish prospective public registration.
