# VarCon
VarCon is a tool to retrieve genomic sequence information from sequence variants. Upon upload of the respective human reference genome fasta file, it reports sequence surroundings and the potential impact of the variant on splicing regulatory sequence elements like SRP binding sites or splice site sequences. 
<br/><br/>
## Installation and Usage
#### 1. Install dependencies
The VarCon tool was developed as an R package (version ≥ R-3.5.2) using integrated Perl scripts (version ≥ Perl 5.18.1) for MaxEntScan score calculations. Therefore, if necessary please install Perl and R first. In total four R packages are required to execute the tool, namely “Biostrings”, “BSgenome”, "seqinr" and "ggplot2". To execute the attached shiny application, following R packages have to be additionally installed: "shiny", "shinyFiles" and "shinycssloaders".
<br/><br/>
In order to execute the VarCon tool, please download every file of the VarCon repository for installation or use the devtools option to install R packages from github.<br/>
`library(devtools)`<br/>
`install_github("caggtaagtat/VarCon")`
<br/><br/>
#### 2. Upload reference genome fasta file
First, please upload the reference genome fasta file, which will be used to retrieve respective surrounding sequences. <br/>
`prepare_reference_fasta(path_to_reference_genome_fasta_file)`
<br/><br/>
#### 3. Get nucleotide sequences around genomic position
Now get the surrounding nucleotide sequence around the sequence variation for later comparison.<br/>
`results <- get_seq_info_from_variation(reference_dnastringset, transcript_id, variation, nt_window, transcript_table, gene2transcript)`
<br/><br/>
#### 4. Generate HZEI plot
To visualize potential impact of the sequence variant on the SRP-binding potential, the HZEI-score per nucleotide is calculated for the sequence with and without the sequence variant. Additionally the intrinisc strength of splice site sequences within the respective sequences is visualized.<br/>
`get_HEXplorer_plot(results, nt_window)`
<br/><br/>
#### 5. Shiny application
To locally start the shiny application, please enter the shiny sub-directory, open the shiny R-script "app.R", change the working directory to source file location and start the app.
<br/><br/>
