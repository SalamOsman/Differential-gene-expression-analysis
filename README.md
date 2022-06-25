# Diffrential gene expression analysis of TCGA datasets.
A comprehensive analysis of differential gene expression matrix from The Cancer Genome Atlas (TCGA) repository.

# Prerequisite for the analysis.
1. The installation of R software and packages for counting raw reads 
2. (Optional) Use the MS windows generic text editor "Notepad" and if you want to use more efficient text editor then install either a "TextPad" or "Notepad++".   

# Source for RNA-seq raw counts datasets and clinical files
* Broad GDAC firehose (https://gdac.broadinstitute.org/)
* Link of the TCGA mRNA preprocessed file used for this tutorial (http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/HNSC/20160128/gdac.broadinstitute.org_HNSC.mRNAseq_Preprocess.Level_3.2016012800.0.0.tar.gz)

# Preparing the datasets
1. Before starting the analysis, always make a separate directory so that they may not replace the files of the previous analysis i.e. "C:/Users/User_name/Desktop/GitHub/HNSCC_DEGs/" on the Desktop.   
2. In this turoial, you will be using the mRNAseq_preprocess archived file of Head and Neck squamous cell carcinoma obtained from the Broad GDAC firehose. Decompress the downloaded file and import the file names as "HNSC.mRNAseq_raw_counts.txt" in R.   
3. For assigning features to the raw counts, also retrieve the clinical archived file that is termed as "Clinical_Pick_Tier1". After decompressing the file, restructure the file named as "All_CDEs.txt" (Use the methodology described in "Creating-a-meta-file-from-TCGA-raw-clinical-file" session). For convenienece, the "meta_data.txt" has been added to the current repository. 
# Running the script
1. Change the directory path to your dataset location. 
2. Customize the name of input and output file names as you like. 
3. Run the script step by step.
4. Let's see the score using statistical cutoff of each gene in the resulting list i.e. try adjusted p-value (padj) <0.05 and <0.02 and also absolute log2FoldChange >1 and >2.
5. Also, get a background knowledge for reading the plots and look at interesting observations.  

# If you need any help, ask me!
abdussalam@precisionmedicine.pk
salam.nurzai@gmail.com


