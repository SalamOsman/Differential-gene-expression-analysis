# Diffrential gene expression analysis of TCGA datasets.
A comprehensive analysis of differential gene expression matrix from The Cancer Genome Atlas (TCGA) repository.

# Prerequisite for the analysis.
1. The installation of R software and packages for counting raw reads 
2. Use the MS windows generic text editor "Notepad" and if you want to use more efficient text editor then install either a "TextPad" or "Notepad++".   

# Source for RNA-seq raw counts datasets and clinical files
* Broad GDAC firehose (https://gdac.broadinstitute.org/)

# Preparing the datasets
1. Before starting the analysis, always make a separate directory so that they may not replace the files of the previous analysis i.e. "C:/Users/User_name/Desktop/GitHub/HNSCC_DEGs/" on the Desktop.   
2. In this turoial, you will be using the mRNAseq_preprocess archived file of Head and Neck squamous cell carcinoma obtained from the Broad GDAC firehose. Decompress the       downloaded file and use the "HNSC.mRNAseq_raw_counts.txt" file as a source for "raw_counts.txt".    
3. For assigning features to the raw counts, also retrieve the clinical archived file that is termed as "Clinical_Pick_Tier1". After decompressing the file, import "All_CDEs.txt" into MS Excel application and transpose the data frame of the sheet. Also remember that you should look and remove    
4. Prepare a meta_data.txt from the clinical 

3. Download mRNAseq_preprocess dataset. In this turoial, we will be using the mRNAseq_preprocess file of Head and Neck squamous cell carcinoma.   
4. In addition, we will also need to download clinical dataset which includes information related to the age, sample type, gender, cancer stage etc of the sample IDs.
5. Create a metafile from the clinical dataset    

# Running the script
1. Change the directory path to your dataset location. 
2. Customize the name of input and output file names as you like. 
3. Run the script step by step.
4. Let's see the score using statistical cutoff of each gene in the resulting list. 

# If you need any help and support ask me!
salam.nurzai@gmail.com


