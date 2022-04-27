'''
Author: Sara Tahir
Date: 04/15/2022
Example usage: python3 frameshift_detector.py input.json --verbose
'''

from contextlib import nullcontext
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import argparse
import json
import time
import os

potential_frameshifts = []
genome_features = []

class GenomeFeature:
    '''
    An object containing information about a specific feature pulled from a Genebank file.
    '''
    def __init__(self, nuc_acc, nuc_desc, feature):
        # Basic feature metadata
        self.accession = nuc_acc
        self.description = nuc_desc
        self.feature = feature
        self.type = feature.type
        self.strand = feature.strand
        self.location = feature.location
        self.locus_tag = feature.qualifiers['locus_tag'][0]
        self.protein_id = feature.qualifiers['protein_id'][0]
        self.product = feature.qualifiers['product'][0]
        # If CDS has multiple exons
        if self.feature.location_operator == 'join':
            self.unspliced_sequence = self.get_unspliced_DNA()
            self.spliced_sequence = self.get_spliced_DNA()[0]
            self.splice_seq_true_positions = self.get_spliced_DNA()[1]
        # If CDS has one exon
        else:
            self.unspliced_sequence = nucrec[self.location.start:self.location.end].seq
            self.spliced_sequence = nucrec[self.location.start:self.location.end].seq
            self.splice_seq_true_positions = list(range(int(self.location.start), int(self.location.end)+1))
    
    def get_unspliced_DNA(self):
        '''
        Returns a sequence including introns for the genome feature
        '''
        first_exon_start = self.location.parts[0].start # start of first exon
        last_exon_end = self.location.parts[-1].end # end of last exon
        return nucrec[first_exon_start:last_exon_end].seq

    def get_spliced_DNA(self):
        '''
        Returns a tuple containing the sequence excluding introns for the genome feature and an
        array of true nucleotide positions for the feature
        '''
        sequence_string = ""
        true_positions = []
        for part in self.location.parts:
            print(part)
            true_positions.extend(list(range(int(part.start), int(part.end)+1)))
            sequence_string += str(nucrec[part.start:part.end].seq)

        print(sequence_string)
        print(true_positions,'\n')
        return (Seq(sequence_string), true_positions)
    

class Frameshift:
    '''
    Holds information for a potential frame shift
    '''
    def __init__(self, fs_case, location, signal_found, genome_feature):
        self.genome_feature = genome_feature
        self.case = fs_case
        self.location = location
        self.signal_found = signal_found
        self.sequence = nucrec[self.start:self.end].seq
    
def download_gbk_files(entrez_email, entrez_api_key, genome_path):
    '''
    Downloads genebank files to genome_path output dir using accession numebrs from input json file
    Parameters:
        entrez_email - entrez email
        entrez_api_key - entrez api key
        genome_path - output directory to save genebank files to
    '''
    Entrez.email=entrez_email
    Entrez.api_key=entrez_api_key
    dir_contents = os.listdir(genome_path)
    #for each chromid in assembly, download GBK, and save using accession number
    print("Downloading GBK files and saving using accession number")
    for chromid_id in params['assembly_chromids']:
        if args.verbose: print("Processing download of: " + chromid_id)
        #check whether the chromid is already in
        if not((chromid_id + '.gb') in dir_contents):
            if args.verbose: print("--> Downloading: " + chromid_id)
            #retrieve genome record
            net_handle = Entrez.efetch(db="nuccore",id=chromid_id,
                                        rettype='gbwithparts', retmode="txt")
            chromid_record=net_handle.read()
            time.sleep(1.5)
            #write the record to file
            out_handle = open(genome_path + '/' + chromid_id + ".gb", "w")
            out_handle.write(chromid_record)
            time.sleep(1.5)
            out_handle.close()
            net_handle.close()

def find_case_one_frameshift(feature):
    '''
    Search for +1 frameshift in given feature
    Parameters:
        feature - search for frameshift for this feature
    '''
    print('Searching for case 1 frameshift')

def find_case_two_frameshift(feature, shift, ustream_limit, stop_codons, signals):
    '''
    Search for -1 frameshift in given feature
    Parameters:
        feature - search for frameshift for this feature
        shift - integer value that determines search shift
        ustream_limit - how far up we should look
        stop_codons - list of stop codons
        signals - list of signals to search for
    '''
    roi_left = 0 # the first +1 destination frame stop codon - needs better initial value
    roi_right = feature.location.end # the stop codon in the source frame
    ustream_count = 0
    # Frameshift + 1 to destination, start at stop + 1, go upstream until we find the first stop codon
    current_pos = feature.location.end + shift
    while ustream_count < ustream_limit and current_pos >= feature.location.start:
        if nucrec[current_pos-3:current_pos].seq in stop_codons:
            roi_left = current_pos # found first +1 destination frame stop codon
            break
        current_pos -= 3
        ustream_count += 3

    # search for signal in roi
    while current_pos < roi_right:
        if str(nucrec[current_pos:current_pos+7].seq) in signals[0][0]:
            if args.verbose: print("Found case 2 frameshift")
            if args.verbose: print(feature)
            if args.verbose: print(nucrec[feature.location.start:feature.location.end].seq)
            frameshift = Frameshift(2, nucrec.name, nucrec.description, current_pos, signals[0][0], feature)
            potential_frameshifts.append(frameshift)
            break
        current_pos += 1

def write_to_txt(output_filename):
    with open(output_filename+'.txt','w') as outfile:
        for fs in potential_frameshifts:
            outfile.write('\n\nAccession: ' + fs.accession)
            outfile.write('\nDescription: ' + fs.description)
            outfile.write('\nCase: ' + str(fs.case))
            outfile.write('\nType: ' + fs.type)
            outfile.write('\nStrand: ' + str(fs.strand))
            outfile.write('\nStart: ' + str(fs.start))
            outfile.write('\nEnd: ' + str(fs.end))
            outfile.write('\nSignal Found: ' + fs.signal_found)
            outfile.write('\nLocus Tag: ' + fs.locus_tag)
            outfile.write('\nProtein ID: ' + fs.protein_id)
            outfile.write('\nProduct: ' + fs.product)
            outfile.write('\nSequence: ' + str(fs.sequence))
        
if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(description='Find some frameshifts')
    # Require input json file
    parser.add_argument('input', nargs='+', help='an input json file with frameshift search parameters')
    # Optional --verbose flag for printing debug statements
    parser.add_argument('--verbose', action='store_true', help='print helpful things')
    # Optional --download_GBK_files flag to download files if not done so already
    parser.add_argument('--download_GBK_files', action='store_true', help='download genebank files if not already in dir')
    parser.set_defaults(verbose=False, download_GBK_files=False)
    args = parser.parse_args()

    print("Opening input File: ", args.input[0])
    with open(args.input[0]) as json_file:  
        params = json.load(json_file)
    
    if args.download_GBK_files:
        download_gbk_files(params['entrez_email'], params['entrez_apikey'], params['genome_path'])

    #for each nucleotide record associated with assembly
    for chromid_id in params['assembly_chromids']:
        nfile = open(params['genome_path'] + '/' + chromid_id + '.gb', "r")
        global nucrec
        nucrec = SeqIO.read(nfile, "genbank")
        nfile.close()

        print("Processing: " + nucrec.id)
        prev_CDS = None
        for i in range(len(nucrec.features)):
            feat = nucrec.features[i]
            if (feat.type == 'CDS'):
                if (feat.strand == 1):
                    genome_features.append(GenomeFeature(nuc_acc=nucrec.name, nuc_desc=nucrec.description, feature=feat))
                    #find_case_two_frameshift(feat, params['frame'], params['ustream_limit'], params['stop_codons'], params['signals'])

    write_to_txt(params['outfile_name'])