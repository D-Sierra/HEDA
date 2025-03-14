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

sample_eplet_assignment <- function(hres_df, df_eplets){
  #Inicialization of the data frame
  eplet_assignment <- data.frame(Sample = hres_df$Sample,
                                 setNames(data.frame(matrix(NA, nrow = nrow(hres_df), ncol = sum(!(names(hres_df) %in% "Sample")))),
                                          paste0(names(hres_df)[!(names(hres_df) %in% "Sample")], "_eplet")))
  
  #Loop to assign eplets to each cell of eplet_assignment
  for (i in 1:nrow(hres_df)){
    for (j in 2:ncol(hres_df)){
      eplet_assignment[i, j] <- paste(hlapro::lookup_eplets(df_eplets, hres_df[i, j])[[1]], collapse = " ")
    }
  }
  #Convert empty cells to NA
  eplet_assignment[eplet_assignment == ""] <- NA
  
  #Obtain single eplets between class I or class II or all of them
  #Obtain the indices of the columns that correspond to class I or class II
  HLA_I_indices <- which(grepl("^[ABC]", names(eplet_assignment)))
  HLA_II_indices <- which(grepl("^D", names(eplet_assignment)))
  HLA_all_indices <- which(grepl("^[ABCD]", names(eplet_assignment)))
  
  #Function for splitting, obtaining unique values and collapsing
  split_unique_collapse <- function(cell) {
    all_values <- paste(cell, collapse = " ")
    unique_vals <- unique(strsplit(all_values, " ")[[1]])
    non_na_unique_vals <- unique_vals[unique_vals != "NA"]
    collapsed_str <- paste(non_na_unique_vals, collapse = " ")
    return(collapsed_str)
  }
  
  #Apply split_unique_collapse function to columns HLA_I_eplets and HLA_II_eplets
  eplet_assignment$HLA_I_eplets <- apply(eplet_assignment[HLA_I_indices], 1, split_unique_collapse)
  eplet_assignment$HLA_II_eplets <- apply(eplet_assignment[HLA_II_indices], 1, split_unique_collapse)
  eplet_assignment$HLA_all_eplets <- apply(eplet_assignment[HLA_all_indices], 1, split_unique_collapse)
  
  return(eplet_assignment)
}
