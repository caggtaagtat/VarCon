# VarCon
VarCon is a tool to retrieve genomic sequence information from sequence variants. Upon upload of the respective human reference genome fasta file, it reports sequence surroundings and the potential impact of the variant on splicing regulatory sequence elements like SRP binding sites or splice site sequences. 
<br/><br/>
## Installation and Usage
#### 1. Install dependencies
The VarCon tool was developed as an R package (version ≥ R-4.0) using integrated Perl scripts (version ≥ Perl 5.18.1) for MaxEntScan score calculations. Therefore, if necessary please install Perl and R first.
<br/><br/>
In order to execute the VarCon tool, please download every file of the VarCon repository for installation, use Bioconductor or use the devtools option to install R packages from github. <br/>
`library(devtools)`<br/>
`install_github("caggtaagtat/VarCon")`
<br/><br/>
#### 2. Upload reference genome fasta file and transcript table
First, please upload a reference genome fasta file, which will be used to retrieve respective surrounding sequences. <br/>
`prepareReferenceFasta(path_to_reference_genome_fasta_file)`
<br/><br/> Now upload a transcript table which shows the same columns as the exemplary short transcript table "transCoord"." Ready to use transcript tables for SNV positions refering to the reference genome GRCh37 or GRCh38 can be downloaded from: https://github.com/caggtaagtat/VarConTables
#### 3. Get nucleotide sequences around genomic position
Now get the surrounding nucleotide sequence around the sequence variation for later comparison.<br/>
`results <- getSeqInfoFromVariation(referenceDnaStringSet, transcriptID, variation, ntWindow=20, transcriptTable, gene2transcript=gene2transcript)`
<br/><br/>
#### 4. Generate HZEI plot
To visualize potential impact of the sequence variant on the SRP-binding potential, the HZEI-score per nucleotide is calculated for the sequence with and without the sequence variant. Additionally the intrinisc strength of splice site sequences within the respective sequences is visualized.<br/>
`generateHEXplorerPlot(results, nt_window=20)`
<br/><br/>
#### 5. Shiny application
To locally start the shiny application, please execute the function `startVarConApp()`
<br/><br/>
