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

#Function calculating the intersection between two eplet vectors
forbidden_eplets <- function(x, eplets){
  x <- strsplit(x[[1]], " ")[[1]]
  resultado_interseccion <- intersect(x, eplets)
  #Check if the result is of size 0
  if (length(resultado_interseccion) == 0) {
    #If there is no coincidence between the two, NA is assigned to the value of the intersection
    resultado_interseccion <- "No eplet missmatch"
  } else {
    #Otherwise it collapses to form a single string
    resultado_interseccion <- paste(resultado_interseccion, collapse = " ")
  }
  return(resultado_interseccion)
}
