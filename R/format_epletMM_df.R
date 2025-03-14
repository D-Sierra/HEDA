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

format_epletMM_df <- function(df_MM, df_eplets){
  #Code to split the calculated prohibited eplets into new columns according to their level of exposure
  #Convert the column df_eplets$exposition to factor
  df_eplets$exposition <- as.factor(df_eplets$exposition)
  #We create new columns in candidates_DL23 according to the levels of exposure and a column for unknown candidates (‘Unknown’)
  for (nivel in c(levels(df_eplets$exposition), "Unknown")) {
    df_MM[[nivel]] <- NA
  }
  
  #Loop that iterates through each row of df_MM and gets the value of the cell epletMM
  for (i in seq_len(nrow(df_MM))) {
    elementos <- unlist(strsplit(as.character(df_MM$MM[i]), " "))
    
    #Loop that iterates over the vector of eplets obtained and looks up the exposure values for that eplet in df_eplets
    for (elemento in elementos) {
      exposicion <- unique(df_eplets$exposition[df_eplets$name == elemento])
      #If the exposure is NA, the value ‘Unknown’ is assigned
      if (is.na(exposicion)) {
        exposicion <- "Unknown"
      }
      #Locates for the row being iterated the column that matches the exposure level and assigns the analysed eplet to it
      if (is.na(df_MM[[exposicion]][i])) {
        df_MM[[exposicion]][i] <- elemento
      } else {
        df_MM[[exposicion]][i] <- paste(df_MM[[exposicion]][i], elemento, sep = " ")
      }
    }
  }
}
