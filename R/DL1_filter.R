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

##################################################################################
#Function for filtering and modifying the format of the table to calculate the DL1
DL1_filter <- function(df_luminex = df_luminex,
                          allele_MM = NA,
                          allele_compl = NA,
                          cut_off_positive = 5000, #Minimum MFI that must have presented a specificity at some point in time to be considered as truly positive
                          cut_off_negative = 5000, #MFI cut-off to consider a specificity negative 
                          time_range = 24, #Number of months during which a specificity must have been maintained below cut_off_negative MFI value to be considered as DL1
                          ignore_DP = TRUE){
  
  ##################################################################
  #Code to read and format input files correctly
  #Assign the first row as column name and delete it
  colnames(df_luminex) <- df_luminex[1,]
  df_luminex <- df_luminex[-1,]
  #We extract from the column names, excluding the first one, the patterns matching date xx/xx/xx or xx/xx/xxxx
  pattern <- ".*?(\\d{1,2}/\\d{1,2}/\\d{2,4}).*"
  column_names <- colnames(df_luminex[-1])
  if(length(which(!grepl(pattern, column_names))) > 0){df_luminex <- df_luminex[-(which(!grepl(pattern, column_names)) + 1)]}  #We remove columns where a date is not recognised in the name
  colnames(df_luminex) <- names(df_luminex) %>% gsub(pattern, "\\1", .) %>% as.Date(., format = "%d/%m/%Y")
  colnames(df_luminex)[1] <- "alleles" #The previous line modifies the first column name, a fixed one is given

  #We dissociate alleles found in pairs in the same row (e.g. ‘DPB1*01:01, DPA1*02:01’).
  delimitadores <- c("; ", ", ", ",") #Possible delimiters between allele pairs that will be recognised
  df_luminex <- df_luminex %>% 
    mutate(alleles = str_replace_all(alleles, paste(delimitadores, collapse = "|"), ";")) %>% #Replace all known delimiters with ‘;’.
    separate_longer_delim(alleles, delim = ";")#Pair allele splitting
  
  #If there are any columns whose cells are empty (they are NAs), they will be deleted
  if (length(which(colSums(is.na(df_luminex)) == nrow(df_luminex))) != 0){
    df_luminex <- df_luminex %>% subset(select = -which(colSums(is.na(.)) == nrow(.)))#Remove empty columns
  }
    
  
  #Extract the allele names from the first column by removing unnecessary characters
  #We identify rows containing in the first column a pattern compatible with an HLA allele
  pattern <- "(A|B|C|DRB1|DRB3|DRB4|DRB5|DQB1|DQA1|DPB1|DPA1)\\*[0-9]{2,3}:[0-9]{2,3}"
  matches <- grepl(pattern, df_luminex$alleles)
  #We replace for each cell of df_luminex$alleles[matches] its value by the pattern found
  df_luminex$alleles[matches] <- gsub(paste0(".*(", pattern, ").*"), "\\1", df_luminex$alleles[matches])
  df_luminex <- df_luminex[matches,]  #We select only rows containing allele-matched patterns
  
  #Extract the numeric value inside each cell to remove annotation
  extract_numeric <- function(x) {
    as.numeric(gsub(".*?(-?\\d+\\.?\\d*).*", "\\1", x)) #Locate a pattern that matches any number that can contain decimals separated by ‘.’
  }
  df_luminex <- df_luminex %>% mutate(across(-alleles, extract_numeric)) #We apply f(x) to all cells except the first column
  

  ################################################################
  #Code for filtering DL1 candidate specificities
  
  #We mark the rows according to whether they have been positive at some point: MFI > cut_off_positive
  filas_con_positivos <- which(apply(df_luminex[-1], 1, function(x) any(x > cut_off_positive)))
  df_luminex$filter <- ifelse(1:nrow(df_luminex) %in% filas_con_positivos, TRUE, FALSE)
  
  #We create a new dataframe containing for each unique allele the number of times it has associated TRUE/FALSE values in df_luminex$filter
  allele_filter <- aggregate(filter ~ alleles, data = df_luminex, FUN = function(x) sum(x == TRUE, na.rm = TRUE))
  allele_filter$false <- aggregate(filter ~ alleles, data = df_luminex, FUN = function(x) sum(x == FALSE, na.rm = TRUE))[, 2]
  colnames(allele_filter) <- c("alleles", "filter_TRUE_count", "filter_FALSE_count")
  #Message containing alleles that can be positive by association to other strands.
  mensaje_advertencia <- paste("The following alleles ", paste(unlist(strsplit(subset(allele_filter, allele_filter$filter_TRUE_count > 0 & allele_filter$filter_FALSE_count > 0)[,"alleles"], " ")), collapse = ", "), 
                               "have simultaneous positive and negative MFIs. This is considered to be due to a reaction against the alpha/beta chain with which they are associated, and therefore a false positive, so they are ignored for analysis.")
  
  #Single alleles with at least one positive and no complete negative history are filtered out
  allele_filter <- subset(allele_filter, allele_filter$filter_TRUE_count > 0 & allele_filter$filter_FALSE_count == 0)
  #We filter the dataframe to keep only alleles that have ever been positive, while ignoring alleles with ambiguous df_luminex
  df_luminex <- df_luminex[which(df_luminex$alleles %in% allele_filter$alleles), -which(colnames(df_luminex) == "filter")]
  
  #We group by unique alleles maintaining the maximum MFI value in each assay for those alleles repeated in the luminex panel
  df_luminex <- aggregate(df_luminex[,-1], by = list(df_luminex$alleles), FUN = max)
  #Alleles are passed to row names and that column is deleted
  rownames(df_luminex) <- df_luminex[,1]
  df_luminex <- df_luminex[,-1]
  
  #We filter the dataframe based on the parameters of time_range and cut_off_negative
  #We get the columns with the results within the indicated period according to time_range
  fecha_dos_anos <- as.Date(colnames(df_luminex[1])) - months(time_range) 
  columnas_dos_anos <- colnames(df_luminex)[as.Date(colnames(df_luminex)) >= fecha_dos_anos] 
  #Filter the dataframe to get only rows with values less than cut-off_negative for the set number of previous months
  df_luminex <- subset(df_luminex, apply(df_luminex[columnas_dos_anos], 1, function(x) all(x < cut_off_negative, na.rm = FALSE)))
  #Filter to remove alleles corresponding to DPA1 and DPB1 loci
  if(ignore_DP == TRUE){df_luminex <- subset(df_luminex, !grepl("^DP", rownames(df_luminex)))}
  
  #############################################################################################################
  #Code for annotation of the resulting dataframe with complement fixation alleles and MM with previous donors
  
  if (nrow(df_luminex > 0)){ #The conditional is added because if the resulting df has nrow = 0 the app will give an error
    df_luminex <- df_luminex %>% 
      mutate(allele_MM = NA, comp_fixation = NA) %>% 
      select(allele_MM, comp_fixation, everything())
    #We extract the banned alleles as a vector and note whether any of them are DL1 candidates
    if (!is.na(allele_MM)){
      allele_MM <- strsplit(gsub(" " , "", allele_MM), ",")[[1]]
      df_luminex$allele_MM <- ifelse(rownames(df_luminex) %in% allele_MM, "Yes", NA)
    }
    
    #We extract the complement binding alleles as a vector and note whether any of them are DL1 candidates
    if (!is.na(allele_compl)){
      allele_compl <- strsplit(gsub(" " , "", allele_compl), ",")[[1]]
      df_luminex$comp_fixation <- ifelse(rownames(df_luminex) %in% allele_compl, "Yes", NA)
    }
  }
  
  return(df_luminex)
}
