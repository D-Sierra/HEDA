# HLA Eplet-based Delisting Assistant (HEDA) is a program designed to 
# facilitate the delisting process of prohibited HLA alleles in highly 
# sensitized patients using eplet incompatibility as the base criterion.
# Copyright (C) 2024  Daniel Álvarez-Sierra and David San Segundo Arribas.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

HLA_dataframe <- function(receptor_HLA, donor1_HLA, donor2_HLA, donor3_HLA, donor4_HLA, donor5_HLA){
  samples <- list(receptor_HLA, donor1_HLA, donor2_HLA, donor3_HLA, donor4_HLA, donor5_HLA)
  names(samples) <- c("receptor", "donor_1", "donor_2", "donor_3", "donor_4", "donor_5")
  #Remove NA samples
  samples <- samples[which(!is.na(samples))]
  #We create a df that will store the quality control of the typings
  HLA_df <- data.frame(matrix(NA, nrow = sum(!is.na(samples)), ncol = 15))
  column_names <- c("Sample", "A1", "A2", "B1", "B2", "C1", "C2", "DRB11", "DRB12", "DRB3451", "DRB3452", "DQB11", "DQB12", "DQA11", "DQA12")
  names(HLA_df) <- column_names
  HLA_df$Sample <- names(which(!is.na(samples)))
  

  #We populate the table obtained with the alleles provided for each sample
  for (i in 1:length(samples)){
    #We iterate the non-NA samples to obtain all the alleles contained in each genotype as a vector
    alleles <- gsub("\\\\t|\\s+", "", strsplit(samples[i][[1]], ",")[[1]])
    for (allele in alleles) {
      #Find the corresponding column for the allele
      if (grepl("^A\\*", allele)){
        column_name <- "A1"
      } else if (grepl("^B\\*", allele)){
        column_name <- "B1"
      } else if (grepl("^C\\*", allele)){
        column_name <- "C1"
      } else if (grepl("^DRB1\\*", allele)){
        column_name <- "DRB11"
      } else if (grepl("^DRB[345]\\*", allele)){
        column_name <- "DRB3451"
      } else if (grepl("^DQB1\\*", allele)){
        column_name <- "DQB11"
      }
      
     #Assign the allele to the first available column corresponding to the HLA type.
     if (is.na(HLA_df[i, column_name])) {
       HLA_df[i, column_name] <- allele
     } else {
         #If the column is already occupied, search for the next available column.
         next_column <- sub("1$", "2", column_name)
         HLA_df[i, next_column] <- allele
      }
    }
  }

  
  for (i in 1:nrow(HLA_df)) {
    for (column in names(HLA_df)[grepl("2$", names(HLA_df))]) {
      if (is.na(HLA_df[i, column])) {
        #Assign the value of the column ending in ‘1’ if the column ending in ‘2’ is blank
        HLA_df[i, column] <- HLA_df[i, sub("2$", "1", column)]
      }
    }
  }
  
  #DQA1 imputation
  HLA_df <- DQA1_imputation(HLA_df)
  

  return(HLA_df)
}
