#frameshift_detector
frameshift_detector identifes both annotated and predicted +1 programmed ribosomal frameshifts through the detection of specific 
heptameric sequences in the genome of an organism, and detects frameshift conservation across specified organisms of interest.

##How it works
###main_input.json + species_and_heptamers.csv -> genome download -> frameshift detection -> frameshift conservation

The following two files are provided as inputs for the frameshift_detector to indicate
main_input.json file layout:
```python
{
    "genome_path" : "genbank_files",                # path to store genbank files
    "fs_inputs_path" : "input_files",               # path for individual species input files containing species metadata
    "frame" : 1,                                    # indicates +1 frameshift
    "ustream_limit" : 10000,                        # upstream search limit 
    "dstream_limit" : 10000,                        # downstream search limit
    "stop_codons" : ["TAG","TAA","TGA"],            # stop codon list
    "start_codons" : ["ATG"],                       # start codon list
    "csv_input_file": "species_and_heptamers.csv",  # input file name for species and heptamers of interest
    "results_dir": "results"                        # path to store frameshift csv results
}
```

species_and_heptamers.csv file layout
Name,Heptamers
Saccharomyces cerevisiae S288C,"[CCTAGGC, 10], [CCTAGTT, 10], [CCTTGAC, 10], [CCCAGGC, 10], [CCCAGTT, 10], [CCCTGAC,10], [CCGAGGC, 10], [CCGAGTT, 10], [CCGTGAC, 10]"

Genome download
Species name -> tax id -> assembly accession number -> gbff file download using ncbi datasets
*Genomes are downloaded regardless of completeness in annotation

For each species, an input file is generated at the fs_inputs_path containing metadata for the frameshift detection step.

Frameshift detection
Insert diagrams :)

After the frameshift detection step, a csv file named "frameshifts.csv" is created in the results_dir and contains frameshifts found for all species of interest. 

frameshifts.csv file header layout
Species,Accession,Description,Locus Tag,Protein ID,Known,Product,Strand,Case,Signal,Signal Score,Frameshift Stop Codon,Annotated Gene Location,Frameshift Location,Annotated Gene Product Length,Frameshift Product Length,Annotated Gene Product,Frameshift Product,Spliced Annotated Gene Sequence,Spliced Frameshift Sequence,Orthology Group

A text file is also created separately for each species and saved into the results directory

Frameshift conservation
A blast database is created from the frameshifted proteins using makeblastdb. Each frameshift is blasted against this database and given an orthology group number that is appended to the original frameshifts.csv file.