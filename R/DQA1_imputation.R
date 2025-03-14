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

DQA1_imputation <- function(df){
  #Table of equivalences for DQB1, DRB1, DQA1 alleles
  equivalence_table <- data.frame(
    DQB1 = c(".*\\*02:01.*", ".*\\*02:02.*", ".*\\*02:02.*", ".*\\*03:01.*", ".*\\*03:01.*",".*\\*03:01.*", ".*\\*03:02.*", ".*\\*03:03.*", ".*\\*03:03.*", ".*\\*03:19.*", ".*\\*04.*", ".*\\*04.*",
             ".*\\*05:01.*", ".*\\*05:02.*", ".*\\*05:03.*", ".*\\*06:01.*", ".*\\*06:02.*", ".*\\*06:03.*", ".*\\*06:04.*", ".*\\*06:08.*", ".*\\*06:09.*"),
    DRB1 = c(NA, ".*\\*04:05.*", "^(?!.*\\*04:05).*$", ".*\\*04:01.*", ".*\\*08:01.*","^(?!.*(?:\\*04:01|\\*08:01)).*$", NA, ".*\\*07:01.*", "^(?!.*\\*07:01).*$", NA,  ".*\\*08.*",
             "^(?!.*\\*08).*$", NA, NA,NA,NA,NA,NA,NA,NA,NA),
    DQA1 = c("DQA1*05:01", "DQA1*03:01", "DQA1*02:01", "DQA1*03:01", "DQA1*06:01","DQA1*05:01", "DQA1*03:01", "DQA1*02:01","DQA1*03:01", "DQA1*05:05", "DQA1*04:01", "DQA1*03:01", "DQA1*01:01", "DQA1*01:02",
             "DQA1*01:01", "DQA1*01:01", "DQA1*01:02", "DQA1*01:03", "DQA1*01:02", "DQA1*01:01", "DQA1*01:01")
  )
  
  #DQA1 assignment based on DQB1 and DRB1 typing
  for (i in 1:nrow(df)) {
    for (j in 1:2) {
      dqb_value <- ifelse(j == 1, df[i, "DQB11"], df[i, "DQB12"]) #For each row, both columns of the DQB1 typing are evaluated separately
      drb_value <- paste(df[i, "DRB11"], ", ", df[i, "DRB12"])   #For each row, the value of the two DRB1 alleles are concatenated for joint evaluation
      
      #For each DQB1 allele, the rows containing that allele are listed in the equivalence table
      index <- which(sapply(equivalence_table[, "DQB1"], function(x) any(grepl(x, dqb_value, perl = TRUE))))
      #If only a single index is obtained, there is no ambiguity in the assignment of DQA1
      if (length(index) == 1) {
        df[i, ifelse(j == 1, "DQA11", "DQA12")] <- equivalence_table[index, "DQA1"]
        #If the index obtained has more than one value then the DQB1 allele can be assigned to more than one DQA1 allele depending on the presence or absence of a particular DRB1 allele
        #By means of regex, the presence of DRB1 alleles is evaluated to resolve the ambiguity in the equivalence table
        #A single index of the equivalence table is obtained and the corresponding value of DQA1 is assigned to the column of the dataframe
      } else if (length(index) > 1) {
        index2 <- sapply(equivalence_table[index, "DRB1"], function(x) any(grepl(x, drb_value, perl = TRUE)))
        if (any(index2)) {
          df[i, ifelse(j == 1, "DQA11", "DQA12")] <- equivalence_table[index[index2], "DQA1"]
        } else {
          next
        }
      } else {
        next
      }
    }
  }
  return(df)
}
