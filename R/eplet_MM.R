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

eplet_MM <- function(eplet_table){
  #Vector with all unique eplets of the receiver
  eplets_receptor <- strsplit(subset(eplet_table, eplet_table$Sample == "receptor")[, "HLA_all_eplets"], " ")[[1]]
  #Vector with unique eplets in the donor set
  filas_donantes <- grep("donor", eplet_table$Sample)
  eplets_donors <- unique(strsplit(paste(eplet_table[filas_donantes, "HLA_all_eplets"], collapse = " "), " ")[[1]])
  #eplet missmatch between recipient and donors, also removes the annotations between []
  eplet_missmatch <- setdiff(eplets_donors, eplets_receptor)
  return(eplet_missmatch)
}