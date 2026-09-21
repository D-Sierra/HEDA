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

#Function to obtain the missmatch at allele level between recipient and previous donors
HLA_MM <- function(receptor_HLA, donor1_HLA = NA, donor2_HLA = NA, donor3_HLA = NA, donor4_HLA = NA, donor5_HLA = NA){
  samples <- list(receptor_HLA, donor1_HLA, donor2_HLA, donor3_HLA, donor4_HLA, donor5_HLA)
  names(samples) <- c("receptor", "donor_1", "donor_2", "donor_3", "donor_4", "donor_5")
  #We create a df that will store the quality control of the typings
  check_df <- data.frame(matrix(0, nrow = sum(!is.na(samples)), ncol = 9))
  check_columns <- c("A", "B", "C", "DRB1", "DRB345", "DQB1", "DQA1", "unknown_alleles", "PASS")
  names(check_df) <- check_columns
  rownames(check_df) <- names(which(!is.na(samples)))
  check_df[,"PASS"] <- FALSE
  
  #REGEX list to assess alleles
  check_regex <- c("A\\*\\d{2}:\\d{2}", "B\\*\\d{2}:\\d{2}", "C\\*\\d{2}:\\d{2}", "DRB1\\*\\d{2}:\\d{2,3}", 
                   "DRB[345]\\*\\d{2}:\\d{2}", "DQB1\\*\\d{2}:\\d{2,3}", "DQA1\\*\\d{2}:\\d{2}")
  #Loop to assess whether all non-NA genotypes follow some of the corresponding patterns with A, B, C, DRB1, DRB345, DQB1, DQA1 alleles.
  for (i in 1:length(samples)){
    if (!is.na(samples[i])){
      alleles <- gsub("\\\\t|\\s+", "", strsplit(samples[i][[1]], ",")[[1]])
      for (allele in alleles){
        matched <- FALSE
        for (j in 1:length(check_regex)){
          if (grepl(check_regex[j], allele)){
            check_df[i,j] <- check_df[i,j] + 1 #If an allele follows a pattern present in check_regex, 1 is added to the value of the corresponding row and column
            matched <- TRUE
          }
        }
        if (!matched){
          check_df[i, "unknown_alleles"] <- check_df[i, "unknown_alleles"] + 1
        }
      }
    }
  }
  #Loop to check whether each provided typing contains a logical number of alleles for each of the loci
  for (i in 1:nrow(check_df)){
    if (check_df[i, "A"] == 2 & 
        check_df[i, "B"] == 2 & 
        check_df[i, "C"] == 2 & 
        check_df[i, "DRB1"] == 2 & 
        check_df[i, "DRB345"] >= 0 & check_df[i, "DRB345"] <= 2 & #DRB345 alleles need not be associated with a haplotype
        check_df[i, "DQB1"] == 2 & 
        check_df[i, "DQA1"] >= 0 & check_df[i, "DQA1"] <= 2 &
        check_df[i, "unknown_alleles"] == 0){ #We accept absent DQA1 alleles
      check_df[i, "PASS"] <- TRUE
    }
  }
  
  #QC checking that all typing is correct
  QC <- all(check_df[, "PASS"])
  
  
  #We obtain the unique alleles for the recipient and all previous donors
  receptor <- unique(gsub("\\\\t|\\s+", "", strsplit(receptor_HLA, ",")[[1]]))
  donor <- unique(gsub("\\\\t|\\s+", "", unlist(lapply(samples[which(grepl("^donor", names(samples)) & !is.na(samples))], function(x) strsplit(x, ",")[[1]]))))
  #Calculation of prohibited alleles that have been in previous donors and not in the recipient
  forbidden_alleles <- paste(sort(setdiff(donor, receptor)), collapse = ", ")
  
  return(list(forbidden_alleles, QC))
}
