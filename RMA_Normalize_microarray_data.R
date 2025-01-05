#    Microarray gene expression data analysis 
# RMA normalization: Prepares normalized expression data from raw microarray data
# CEL file: Cell intensity file, probe level PM and MM values

# BiocManager::install("GEOquery")  # Get data from NCBI Gene Expression Omnibus (GEO)
# BiocManager::install("affy")      # Methods for Affymetrix Oligonucleotide Arrays
# install.packages("tibble")

#-------------------------------------------------------------------------------
### Note ###
# These process was performed separately for each dataset to ensure individual normalization before merging
#-------------------------------------------------------------------------------

# Load necessary libraries
library(tidyverse)
library(GEOquery)
library(affy)
library(dplyr)
library(tibble)

# get dataset from GEO: Supplementary files
getGEOSuppFiles("GSE107465")         # Take DatasetID
# untar files
untar("GSE148537/GSE148537_RAW.tar", exdir='data/')      # path
# Read raw data from CEL files
raw_data <- ReadAffy(celfile.path = "Leukemia/GSE51082/GSE51082_RAW/")

# Perform RMA normalization
normalized_data  <- rma(raw_data)
str(normalized_data)
normalized_data <- read.csv("../GSE51082_normalized_expression.csv")      # path for normalized dataset

# Extract expression values
normalized_expression <-as.data.frame(exprs(normalized_data))

# map probe IDs to gene symbols 
gse <- getGEO("GSE107465", GSEMatrix = TRUE)          # Dataset ID

# fetch feature data to get ID - gene symbol mapping
feature_data <- gse$GSE107465_series_matrix.txt.gz@featureData@data       # retrieve from GEO database
# subset ID and gene symbol columns
feature_data <- feature_data[, c(1,11)]           

# Merge the the data with gene symbol
normalized_expression <- normalized_expression %>%
  rownames_to_column(var = 'ID') %>%
  inner_join(., feature_data, by = 'ID') 


# Remove the first two columns
normalized_expression <- normalized_expression[, -c(1)]
normalized_expression <- t(normalized_expression)
normalized_expression <- as.data.frame(normalized_expression)

# Extract the last row as the header
header <- as.character(normalized_expression[nrow(normalized_expression), ])
# Remove the last row
normalized_expression <- normalized_expression[-nrow(normalized_expression), ]
colnames(normalized_expression) <- header
# Convert row names to a column named "Sample ID" : "as in metadata"
normalized_expression$'Sample.ID' <- rownames(normalized_expression)
# Reorder the column to the first position (Sample ID)
normalized_expression <- normalized_expression[, c(ncol(normalized_expression), 1:(ncol(normalized_expression) - 1))]
rownames(normalized_expression) <- NULL
# Rename the first column to "Sample ID"
colnames(normalized_expression)[1] <- "Sample.ID"


# Remove a specific substring ("_suffix")
normalized_expression$'Sample.ID'  <- gsub("_.*", "", normalized_expression$'Sample.ID')


# Save the normalized data
write.csv(normalized_expression, "BloodType_normalized_expression.csv")
