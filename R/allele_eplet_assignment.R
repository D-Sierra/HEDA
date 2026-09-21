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

allele_eplet_assignment <- function(allele_df, df_eplets){
  #Initialisation of the data frame
  eplet_assignment <- data.frame(Alleles = allele_df$Alleles,
                                 Eplets = NA)
  
  #Loop to assign eplets to each cell of eplet_assignment
  for (i in 1:nrow(allele_df)){
    eplet_assignment[i, 2] <- paste(hlapro::lookup_eplets(df_eplets, allele_df[i,])[[1]], collapse = " ")
  }
  #Convert all empty cells to NA
  eplet_assignment[eplet_assignment == ""] <- NA
 
  return(eplet_assignment)
}
