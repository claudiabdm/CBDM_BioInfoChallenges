# CBDM_BioInfoChallenges
A code repository for the Bioinformatics Programming Challenges course of the Master's degree in Computional Biology. 

## Assignment 1:

### Description
Program to simulate the planting of a number of seeds in a stock database and save the updated data to a new file. It also reports the genes that are genetically linked.
  - Input files: gene_information.tsv seed_stock_data.tsv cross_data.tsv new_file_data.tsv
  - Output file: new_file_data.tsv
  
### Command

 ```
 ruby process_database.rb gene_information.tsv seed_stock_data.tsv cross_data.tsv new_file_data.tsv
 ```

Last row of each data file was added to test the program.

## Assignment 2:

### Gems needed

```
gem install rest-client
gem install ruby-progressbar
```

### Description

This program will find protein-protein interactions between genes from a gene codes list and will predict possible interaction networks. It will return a report file with all the interaction networks with KEGG and GO functional annotations for the network and for members. 
  - Input files: gene_codes_list.txt name_report_file.txt
  - Output file: name_report_file.txt
  
 ### Command
  
 ```
 ruby main.rb gene_codes_list.txt name_report_file.txt
 ```

## Assignment 3:

### Gems needed

```
gem install 'rest-client'
gem install 'bio'
gem install 'ruby-progressbar'
```

### Description

This program will find the desired motif in the exons of genes from a list of AGI codes and will execute a GFF3 file recording each repeat. This file can be used to display the repeats in a genome browser.
  - Input files: gene_codes_list.txt
  - Output file: genes_motif.txt genome_motif.txt genes_without_motif_in_any_exon.txt
  
 ### Command
 
 It can be ran in the terminal. There is also a jupyter notebook available (the Assingment3.pynb).
  
 ```
 ruby main.rb gene_codes_list.txt
 ```
