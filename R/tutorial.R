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

tutorial <- function(){
  fluidRow(width = 12,
    column(width = 10,
    tags$div(style = "font-size: 16px", 
      tags$p("This tool is intended to assist in the eplet-based delisting process of highly sensitised organ recipients.")),
    tags$h3(tags$b("Data input")),
    tags$h4(tags$b("HLA typing")),
    div(style = "font-size: 16px", 
      tags$p("In this section, the recipient's and previous donors' HLA typings are entered in the corresponding fields. Considerations to be taken into account:"),
      tags$li("Typings should be high resolution, in case low resolution genotyping is available, it should be previously imputed using other tools such as ", 
              tags$a( href = "https://www.haplostats.org", "Haplostats", target = "_blank"), "."),
      tags$li("Although this platform does not perform low to high resolution typing conversion, it does perform DQA1 locus imputation which is not usually provided by other applications."),
      tags$li("Alleles must be separated by commas (,)."),
      tags$li("In the case where only a single allele is provided for any of the loci, the individual is considered homozygous for that locus."),
      tags$li("The application checks if all elements follow a pattern compatible with the designation of an HLA allele, and will report most of the alleles with an incorrect format or with some kind of typing error.
                  However, the user should check if the typing recognition and the imputation of the DQA1 locus are correct in the corresponding table of the Data input - Summary section.")),
    tags$h4(tags$b("Delisting 1 (DL1) file input")),
    tags$div(style = "font-size: 16px", 
      tags$p("This section allows the attachment of a file containing the recipient's historical anti-HLA antibody data for analysis according to the parameters chosen in section Delisting 1.
              The attached file must meet the following requirements:"),
      tags$li("The file provided must be in text format (.txt) separated by tabs."),
      tags$li("The first column should contain the names of the specificities in high resolution.  
              Any set of characters that do not correspond to an HLA allele at high resolution will be removed from the corresponding cell."),
      tags$li("In case a row contains information of two alleles that are cuantified by the same Luminex bead, it will be split into two identical rows with the same MFI values provided that the alleles are separated by the symbols ; or , 
              (e.g. a row with the following value in the first column ‘DPA1*01:05; DPB1*28:01’ will be duplicated and the resulting rows will be assigned to ‘DPA1*01:05’ and ‘DPB1*28:01’ alleles respectively)."),
      tags$li("The remaining columns, one for each time point, should contain MFIs values for the corresponding allele. The numbers should follow an Anglo-Saxon notation and use the dot (.) as decimal separator.
               Any set of characters that do not correspond to a number will be removed from these cells."),
      tags$li("The first row corresponds to column names, which must contain the date of the test in one of the following formats: dd/mm/yyyy or dd/mm/yyyy.
               Any set of characters that do not correspond to a date in one of the previous formats will be removed."),
      tags$li("If for any row or column a valid allele or date value is not identified, these will be ignored for further analysis.
              In both cases, an error message will be displayed indicating the row or column name of the original file that has been ignored and the value of the corresponding cells where a correct allele or date value was not found.")),
    tags$br(),
    tags$div(
      align = "center",
      tags$img(src = "tableDL1_example.png", height = "150px", width = "auto"),
      tags$figcaption("Example of the header of an excel table before exporting it to a tab-delimited text file. 
                      The ‘Ac Anti-HLA’ prefixes in the first column will be removed and only the allele names will be kept. 
                      For the column names only the dates will be kept, in this case in dd/mm/yy format and numbers representing laboratory IDs will be removed. 
                      The remaining cells contain the MFI value together with the qualitative interpretation of the result (P) or (N), these annotations will be removed and only the numeric values will be kept. 
                      The rows containing the joint analysis of alpha and beta chains of a locus will be duplicated generating a row for each individual chain.")
    ),
    tags$br(),
    tags$h4(tags$b("Complement-fixing DSA")),
    div(style = "font-size: 16px",
        tags$p("This section allows manual selection of all specificities with complement fixing capability. 
                This is to be noted during the DL1 process so that this information can be taken into account for the final selection of the specificities that can be delisted.")),
    tags$h3(tags$b("Delisting 1")),
    div(style = "font-size: 16px",
        tags$p("In this section you can configure the parameters to filter the list of DL1 candidate alleles.The programme will evaluate the previously attached file, and if it complies with the established format requirements,
               will return a table with the candidate specificities that meet the selected criteria. To facilitate the interpretation, it will indicate for each allele if it is present in previous donors, or if it is complement fixer if you have also provided this information.
               Finally, a graph will be created showing for each candidate allele the evolution of its MFI over time."),
        tags$p("A PDF report of the analysis can be downloaded using the ‘Download DL1 results’ button. 
               The report contains the parameters used for the filtering of the candidate specificities and the table and plot showing those that have been selected by the application.")),
    tags$h3(tags$b("Delisting 2/3")),
    div(style = "font-size: 16px",
        tags$p("Delisting stages 2/3 based on eplet incompatibility require a more thorough review of the data and the programme does not incorporate a fully automatic analysis tool to select candidate alleles.
               Instead, a form is provided to manually select candidate specificities for these delisting steps. 
               The configuration section allows to select on which subset of eplets we want to restrict our analysis according to whether or not they have been verified by antibodies and according to their degree of exposure.
               The application will assign each candidate specificity to one of the two delisting phases: DL2 when the allele does not share eplets with the list of eplets prohibited for the receptor, or to the DL3 phase when it does."),
        tags$p("A PDF report of the analysis can be downloaded using the ‘Download DL2/3 results’ button.  
               The report contains the parameters used for the filtering of eplets that have been considered for the analysis, as well as the table of results shown in this section."))
  ))
  
}