# PRISM Repurposing Public 24Q2

##############
## Overview ##
##############

This data release contains two most recent PRISM Repurposing screens: Repurposing-1M and Repurposing-300. All Repurposing-1M [REP1M] and Repurposing-300 [REP300] compounds (1514 = 1280 REP1M + 234 REP300) were screened in the PRISM assay at a dose of 2.5 μM with a 5-day treatment against 906 cancer cell lines (859 of them passed quality checks -QC- for all tested compound with two high quality replicates). Two PRISM cell line collections were used in the assay: PR500A, which includes only adherent cell lines, and PR500B, which has adherent and suspension cell lines. Together, these PRISM cell line collections form the PR1000 cell line collection. All compounds were run in triplicate, and each plate contained positive (Bortezomib, 20μM) and negative (DMSO) controls. The screen can be considered an extension of the PRISM Repurposing Primary Screen, with PR500A mainly covering the existing cell line panel, while PR500B extends the cell line collection with new subtypes and lineages. For the assay details, please refer to Corsello et al., 2020 (doi.org/10.1038/s43018-019-0018-6) and https://www.theprismlab.org.


###########
## Files ##
###########

### Repurposing_Public_24Q2_Cell_Line_Meta_Data.csv

Contains the meta data about the  cell lines  screened in the assay, a list of columns with simple descriptions are given below:

* ccle_name / depmap_id : CCLE and arxspan identifiers for each cell line, please use these to join with other depmap datasets.
* pool_id: Pool identifier of each cell line, please note that some cell lines are included in more than one pool.
* culture: The cell set for each cell line. Takes values PR500A or PR500B. Please note, some cell lines are included in both cell sets.
* row_id: The unique identifier for each row of the raw LMFI values (Repurposing_Public_24Q2_LMFI_matrix.csv), please use this column to join with other data tables in this release.
* screen: Takes the values of "REP1M" and "REP300" to represent the particular run that the data is generated.  


### Repurposing_Public_24Q2_Treatment_Meta_Data.csv

Contains the meta data about the treatment/perturbation in each well of each plate. A list of column names with simple descriptions are given below:
* profile_id: The unique identifier for each prism_replicate (detection plate), perturbation well, and screen. Please use this column to join with the other data tables. This column also represents the column names of Repurposing_Public_24Q2_LMFI_matrix.csv.
* culture: The cell set for each prism_replicate. Takes values PR500A or PR500B. 
* compound_plate: Unique identifier for each compound plate.
* prism_replicate: Unique identifier for each detection plate (3 for each compound plate x culture pair)
* rep: Replicate code for each prism_replicate, takes values of X1, X2 and X3.
* perturbation_type: Takes one of the following 3 values: ctl_vehicle (negative control/DMSO), trt_cp (test compound/treatment), trt_poscon (positive control/Bortezomib) 
* perturbation_well: Physical location of each perturbation on each detection plate (prism replicate)
* dose: Dose of the treatment (in uM).
* broad_id: Unique identifier for each compound. Please note, there might be multiple broad_id's corresponding to each compound but there is only a single compound for each broad_id. 
* name: Name of the compound treatment, according to the Drug Repurposing Hub (https://clue.io/repurposing)
* screen: Takes the values of "REP1M" and "REP300" to represent the particular run that the data is generated.


### Repurposing_Public_24Q2_LMFI_matrix.csv

Contains the log2 median fluorescence intensity values for each cell line (or control barcode) and perturbation well. This is the raw data coming from the luminex detection.


### Repurposing_Public_24Q2_LMFI_NORMALIZED.csv

Contains the log2 median fluorescence intensity values before and after normalizing them to spike in barcodes. List of column names follows: 
* row_id: The unique identifier for each row of the raw LMFI values (Repurposing_Public_24Q2_LMFI_matrix.csv), please use this column to join with other data tables in this release.
* profile_id: The unique identifier for each prism_replicate (detection plate), perturbation well, and screen. Please use this column to join with the other data tables. 
* LMFI: log2 median fluorescence intensity values for each cell line perturbation pair.
* LMFI.ref: Reference values for the control barcodes used for normalization. These values are computed as the centered median LMFI values among negative control wells for each spike-in barcode in each detection plate (prism replicate).
* LMFI.normalized: LMFI values after normalizing them using the control barcodes (spike-in's) to make all wells on the same plate are in the same scale.
* LMFI.normalized.se: The standard errors for the normalization step, these rough values are not used in the rest of the pipeline but can be examined for more detailed troubleshooting.


### Repurposing_Public_24Q2_QC_table.csv

Contains the QC metrics for each cell line x detection plate pairs. A list of columns and brief explanations are given below.
* row_id: The unique identifier for each row of the raw LMFI values (Repurposing_Public_24Q2_LMFI_matrix.csv), please use this column to join with other data tables in this release.
* culture: The cell set for each prism_replicate. Takes values PR500A or PR500B. 
* compound_plate: Unique identifier for each compound plate.
* prism_replicate: Unique identifier for each detection plate (3 for each compound plate x culture pair)
* screen: Takes the values of "REP1M" and "REP300" to represent the particular run that the data is generated.
* error_rate : Bayes error while classifying negative and positive control wells for each cell line x detection plate pair.
* NC.median / PC.median : The median normalized LMFI values for negative/positive control wells in a given detection plate (prism_replicate) 
* NC.mad / PC.mad : The median absolude deviation (MAD) of normalized LMFI values for negative/positive control wells in a given detection plate (prism_replicate) 
* SSMD : Strictily standardized mean difference, defined as (ctl_vehicle_md-trt_poscon_md)/sqrt(ctl_vehicle_mad^2+trt_poscon_mad^2)
* DR : Dynamic range of the assay, defined as the median LFC of high dose bortezomib, i.e. ctl_vehicle_md - trt_poscon_md
* PASS : TRUE iff error_rate < .05 and DR > 2. Flags a QC failure if it is FALSE.
* n.PASS : How many replicates of the corresponding compound_plate x row_name pair passes the QC. Cell lines with n.PASS < 2 are filtered out for the rest of the analysis.


### Repurposing_Public_24Q2_LFC.csv

Contains the log2-fold changes between the perturbations and the median of the negative control wells in each plate. A brief description of the column names is given below.
* row_id: The unique identifier for each row of the raw LMFI values (Repurposing_Public_24Q2_LMFI_matrix.csv), please use this column to join with other data tables in this release.
* profile_id: The unique identifier for each prism_replicate (detection plate), perturbation well, and screen. Please use this column to join with the other data tables. 
* LFC : log2-fold-change values between the normalized LMFI values for each treatment and the median of negative control wells in each detection plate. 
* LFC_cb : LFC values after applying a batch correction step (ComBat) for potential pooling artifacts. 
* PASS : Binary flag to decide if this value will be used in the final table/matrix or not. (See Repurposing_Public_24Q2_QC_table.csv)
          


### Repurposing_Public_24Q2_LFC_COLLAPSED.csv

Log2-fold-change viability values (LFC_cb) computed in Repurposing_Public_24Q2_LFC.csv file are filtered based on the QC results and median collapsed across replicates. The list of the columns in thes data table are described below.
* row_id: The unique identifier for each row of the raw LMFI values (Repurposing_Public_24Q2_LMFI_matrix.csv), please use this column to join with other data tables in this release.
* broad_id: Unique identifier for each compound. Please note, there might be multiple broad_id's corresponding to each compound but there is only a single compound for each broad_id. 
* dose: Dose of the treatment (in uM).
* compound_plate: Unique identifier for each compound plate.
* screen: Takes the values of "REP1M" and "REP300" to represent the particular run that the data is generated.
* culture: The cell set for each prism_replicate. Takes values PR500A or PR500B. 
* LFC: Median collapsed LFC_cb values from Repurposing_Public_24Q2_LFC.csv (excluding rows with PASS = FALSE).


### Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv

The final LFC values in Repurposing_Public_24Q2_LFC_collapsed.csv table are cast into a matrix with rows corresponding to individual treatments and columns are depmap_id's. This final matrix also contains the compounds from the original PRISM Repurposing Primary screen for convenience. Please use Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv to get the meta-data for each row.


### Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv

Metadata for the perturbations represented by each row of Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv. Columns of this table are given below:
* IDs: Unique identifier to join with the rows of Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv
* screen: Identifier for in which Repurposing screen this particular profile is created.
* dose: Dose of the treatment (in uM).
* Drug.Name: Human read-able name for the drug.       
* repurposing_target: Annotated target of the drug.
* MOA: Annotated mechanism of action for the drug.
* Synonyms: Known synonyms of the drug name.


### Repurposing_Public_24Q2_Readme.txt

Description of all files

