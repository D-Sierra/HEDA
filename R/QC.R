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

###########################################################################
#QC function for alleles introduced as recipient and donor types
check_alleles <- function(typing) {
  #We convert the typing into a vector by removing spaces and tabs.
  alleles <- gsub("\\\\t|\\s+", "", strsplit(typing[[1]], ",")[[1]])
  #Patterns that must be met by any valid allele
  patterns <- c("^A\\*[0-9]{2}:[0-9]{2,3}$",
                "^B\\*[0-9]{2}:[0-9]{2,3}$",
                "^C\\*[0-9]{2}:[0-9]{2,3}$",
                "^DRB1\\*[0-9]{2}:[0-9]{2,3}$",
                "^(DRB3|DRB4|DRB5)\\*[0-9]{2}:[0-9]{2,3}$",
                "^DQB1\\*[0-9]{2}:[0-9]{2,3}$",
                "^DQA1\\*[0-9]{2}:[0-9]{2,3}$")
  
  #Logical vector that will store the result of checking whether each allele of the typing follows a valid pattern
  checks <- logical(length(alleles))
  #Loop to iterate each allele
  for (i in seq_along(alleles)) {
    match_found <- FALSE
    #Loop to iterate through each pattern to check if it is found on the allele
    for (j in seq_along(patterns)) {
      if (grepl(patterns[j], alleles[i], perl = TRUE)) {
        match_found <- TRUE #If pattern j is found on the allele, match_found change to TRUE
        break
      }
    }
    #Store the result in the vector checks
    checks[i] <- match_found
  }
  
  #Output
  if (all(checks)) {
    return(NULL) #If all alleles are valid the result is null
  } else {
    #If any allele does not follow a valid pattern it is stored in a vector and concatenated to be displayed as an error message.
    non_matching_alleles <- alleles[!checks]
    return(paste("The following alleles entered are not valid: ", paste(non_matching_alleles, collapse = ", ")))
  }
}

###################################################################
#Quality control function for alleles in the table for DL1
DL1_QC_alleles <- function(df_luminex = df_luminex){
  colnames(df_luminex) <- df_luminex[1,]
  df_luminex <- df_luminex[-1,]
  colnames(df_luminex)[1] <- "alleles"
  pattern <- "(A|B|C|DRB1|DRB3|DRB4|DRB5|DQB1|DQA1|DPB1|DPA1)\\*[0-9]{2,3}:[0-9]{2,3}"
  #We identify which rows of the original table do not contain valid allele names in the first column.
  unr_alleles_index <- which(!grepl(pattern, df_luminex$alleles)) #Row indexes
  unrecognised_alleles <- df_luminex$alleles[unr_alleles_index] #Row values
  #We modify the values by creating a vector indicating the original unrecognised values and their original row to make it easier to search for them
  for (i in seq_along(unrecognised_alleles)){
    unrecognised_alleles[i] <- paste0(unrecognised_alleles[i], " (fila ", unr_alleles_index[i], ")")
  }
  if (length(unrecognised_alleles) == 0){
    resultado <- NULL
  } else {
    resultado <- paste("Values not recognised as alleles in the first column: ",paste(unrecognised_alleles, collapse = ", "))
  }
  return(resultado)
}

###################################################################
#Quality control function for alleles in the table for DL1
DL1_QC_dates <- function(df_luminex = df_luminex){
  colnames(df_luminex) <- df_luminex[1,]
  df_luminex <- df_luminex[-1,]
  pattern <- ".*?(\\d{1,2}/\\d{1,2}/\\d{2,4}).*"
  column_names <- colnames(df_luminex[-1])
  #We identify which column names, excluding the first one, do not contain a correct date format.
  unr_dates_index <- which(!grepl(pattern, column_names)) #Column indexes
  unrecognised_dates <- column_names[unr_dates_index] #Column values
  #We modify the values by creating a vector indicating the original unrecognised values and their original row to make it easier to search for them
  for (i in seq_along(unrecognised_dates)){
    unrecognised_dates[i] <- paste0(unrecognised_dates[i], " (column ", unr_dates_index[i]+1, ")")
  }
  if (length(unrecognised_dates) == 0){
    resultado <- NULL
  } else {
    resultado <- paste("Column names not recognised as dates: ",paste(unrecognised_dates, collapse = ", "))
  }
  return(resultado)
}
