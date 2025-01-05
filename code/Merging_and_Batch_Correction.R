# Merge and Batch correction the datasets to correct the technical variation between the datasets

# 1- Merging the metadata 
library(dplyr)

# Load all metadata files
MM <- read.csv("GSE19784_metadata.csv")
Lymphoma <- read.csv("GSE132929_metadata.csv")
Leukemia_AML_1 <- read.csv("GSE34860_AML_metadata.csv")
Leukemia_AML_2 <- read.csv("GSE107465_metadata.csv")
Leukemia_CLL <- read.csv("GSE21029_metadata.csv")
Leukemia <- read.csv("GSE51082_metadata.csv")
Leukemia_ALL <- read.csv("GSE28497_metadata.csv")

# Combine all metadata files row-wise
All_metadata <- bind_rows(MM, Lymphoma, Leukemia_AML_1,
                          Leukemia_AML_2, Leukemia_CLL, Leukemia, Leukemia_ALL)

# Save the metadata file: All
write.csv(All_metadata, "Merged_Metadata.csv")

#-------------------------------------------------------------------------------
# 2- Merging the count matrix 

# Import the count matrix of each dataset
MM <- read.csv("MM_GSE19784_normalized.csv", row.names = 1)  
lymphoma <- read.csv("Lymphoma_GSE132929_normalized.csv", row.names = 1)
Leukemia <- read.csv("Leukemia_GSE51082_normalized.csv",  row.names = 1)
Leukemia_ALL <- read.csv("Leukemia_ALL_GSE28497_normalized.csv", row.names = 1)
Leukemia_CLL <- read.csv("Leukemia_CLL_GSE21029_normalized.csv")   
Leukemia_AML_1 <- read.csv("Leukemia_AML1_GSE107465_normalized.csv")
Leukemia_AML_2 <- read.csv("Leukemia_AML2_GSE34860_normalized.csv")

# rename the first column in Leukemia_CLL /Leukemia_AML_2/ Leukemia_AML_1 datasets
colnames(Leukemia_CLL)[1] <- "Sample.ID"
colnames(Leukemia_AML_2)[1] <- "Sample.ID"
colnames(Leukemia_AML_1)[1] <- "Sample.ID"

# Rows: sample IDs.
# Columns: gene expression values.

# Create a list of column names from all datasets
datasets <- list(
  colnames(Leukemia),
  colnames(Leukemia_AML_1),
  colnames(Leukemia_AML_2),
  colnames(Leukemia_ALL),
  colnames(Leukemia_CLL),
  colnames(lymphoma),
  colnames(MM)
)

print(colnames(lymphoma)[0:10])

# Identify the common column names (genes) across all datasets
comm_genes <- Reduce(intersect, datasets)
length(comm_genes)


# Subset each dataset to include 'only' the common genes
filtered_datasets <- lapply(datasets, function(df) df[, comm_genes])


# Merge all datasets by row-binding
Merged_CountMatrix <- do.call(rbind, filtered_datasets)


# Save the merged dataset
write.csv(Merged_CountMatrix, "Final Files/All_datasets_merged_data.csv", row.names = TRUE)


#-------------------------------------------------------------------------------------------------
###    Batch correction: to handel the issue of technical variation between the datasets

# BiocManager::install("sva")
library(sva)
library(limma)

all_dataset <- read.csv('All_datasets_merged_data.csv', row.names = 1)

# Prepare the data for batch correction analysis:
# metadata$Sample.ID for batch correction
metadata <- read.csv("Microarray data preparation/Final Files/Metadata.csv")
batch_info <- metadata$Sample.ID  # Extract batch info


# 'ComBat' method: genes as rows and samples as columns.
expression_data <- t(all_dataset)
expression_data <- as.data.frame(expression_data)  
colnames(expression_data) <- rownames(all_dataset)

# Assign the first row as column names
colnames(expression_data) <- expression_data[1, ]
# Remove the first row
expression_data <- expression_data[-1, ]

# Perform Batch Correction
corrected_data <- ComBat(
  dat = expression_matrix,
  batch = metadata$Sample.ID,
  mod = mod1,
  par.prior=TRUE,
)

# Combine the metadat with the df
corrected_data_transposed <- t(expression_data)
corrected_data_transposed <- as.data.frame(corrected_data_transposed)

corrected_data_transposed <- merge(corrected_data_transposed, 
                                   metadata, by = "Sample.ID")


# Save corrected data
write.csv(corrected_data_transposed, "Final_CountMatrix_after_BatchCorrection.csv")












