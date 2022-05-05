'''
Author: Sara Tahir
Date: 04/15/2022
Purpose: To identify new +1 programmed ribosomal frameshifts through the detection of specific heptameric sequences
Example usage: python3 frameshift_detector.py input.json --verbose 
'''

from termcolor import colored
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import argparse
import json
import csv
import time
import os

detected_frameshifts = [] # Holds all Frameshift objects created after successful heptamer detection
genome_features = [] # Holds all GenomeFeature objects created from features read in from genbank files

class GenomeFeature:
    '''
    Contains information about a specific feature pulled from a Genbank file.
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
            if args.verbose: print('Feature has multiple exons')
            self.unspliced_sequence = self.get_unspliced_DNA()
            self.spliced_sequence = self.get_spliced_DNA()[0]
            self.splice_seq_true_positions = self.get_spliced_DNA()[1]
        # If CDS has one exon
        else:
            self.unspliced_sequence = nucrec[self.location.start:self.location.end].seq
            self.spliced_sequence = nucrec[self.location.start:self.location.end].seq
            self.splice_seq_true_positions = list(range(int(self.location.start), int(self.location.end)+1))
        # Set downstream region (x nucleotides after annotated stop)
        self.downstream_region = self.get_downstream_region()[0]
        self.downstream_region_true_pos = self.get_downstream_region()[1]
    
    def get_unspliced_DNA(self):
        '''
        Remove introns from the DNA sequence for the genome feature
        Return:
            unspliced_seq (seq object): unspliced DNA sequence
        '''
        first_exon_start = self.location.parts[0].start # start of first exon
        last_exon_end = self.location.parts[-1].end # end of last exon
        unspliced_seq = nucrec[first_exon_start:last_exon_end].seq
        return unspliced_seq

    def get_spliced_DNA(self):
        '''
        Combines all exons into one sequence for the genome feature
        Return:
            (seq object, [int]): tuple containing the spliced sequence [0] and an int array of 
            the true positions for the nucrec
        '''
        sequence_string = ""
        true_positions = []
        for part in self.location.parts:
            #print(part)
            true_positions.extend(list(range(int(part.start), int(part.end)+1)))
            sequence_string += str(nucrec[part.start:part.end].seq)
        #print(sequence_string)
        #print(true_positions,'\n')
        return (Seq(sequence_string), true_positions)
    
    def get_downstream_region(self):
        '''
        Adds x downstream nucleotides to an array that can be indexed later
        Return:
            (seq object, [int]): tuple containing the downstream sequence [0] and an int array of 
            the true positions for the nucrec
        '''
        true_positions = []
        last_exon_end = self.location.parts[-1].end # end of last exon
        true_positions.extend(list(range(int(last_exon_end), int(last_exon_end)+10000)))
        return (nucrec[last_exon_end:last_exon_end+10000].seq, true_positions)

    def get_true_pos(self, spliced_sequence_pos):
        '''
        Returns nucrec position given a spliced sequence position
        Parameters:
            spliced_sequence_pos (int): position in the spliced DNA sequence
        Return:
            true_pos (int): nucrec position
        '''
        combined_sequence_pos = self.splice_seq_true_positions + self.downstream_region_true_pos
        return combined_sequence_pos[spliced_sequence_pos]



class Frameshift:
    '''
    Holds information for a detected frameshift
    '''
    def __init__(self, genome_feature, signal_found, fs_case, heptamer_location, start_pos, stop_pos):
        self.genome_feature = genome_feature
        self.signal_found = signal_found
        self.case = fs_case
        self.heptamer_location = heptamer_location
        self.start_pos = start_pos
        self.stop_pos = stop_pos
        self.seq_end = 0
        self.stop_codon = 'None'
        self.get_frameshifted_seq()

    def get_original_seq(self):
        '''
        Adds space between 3 base pairs in original sequence
        Return:
            (string): string for the original sequence with spacing between codons
        '''
        seq_string = str(self.genome_feature.spliced_sequence)
        return ' '.join(seq_string[i:i+3] for i in range(0,len(seq_string),3))

    def get_frameshifted_seq(self):
        '''
        Returns a string of the frameshifted sequence with a space between 3 base paids
        Return:
            (string): string for the frameshifted sequence with spacing between codons
        '''
        frameshifted_seq = ''
        frameshifted_translate = ''
        extended_sequence = str(self.genome_feature.spliced_sequence) + str(self.genome_feature.downstream_region)
        if self.case == 'Downstream':
            cur_pos = self.start_pos
            while extended_sequence[cur_pos:cur_pos+3] not in params['stop_codons'] and (cur_pos-self.start_pos < params['ustream_limit']):
                if cur_pos == self.heptamer_location:
                    frameshifted_seq += str(extended_sequence[cur_pos:cur_pos+3]) + ' '
                    frameshifted_translate += str(Seq(extended_sequence[cur_pos:cur_pos+3]).translate())
                    frameshifted_seq += str(extended_sequence[cur_pos+3:cur_pos+4]) + ' '
                    cur_pos += 4
                else:
                    frameshifted_seq += str(extended_sequence[cur_pos:cur_pos+3]) + ' '
                    frameshifted_translate += str(Seq(extended_sequence[cur_pos:cur_pos+3]).translate())
                    cur_pos += 3
                
            frameshifted_seq += str(extended_sequence[cur_pos:cur_pos+3])
            if str(extended_sequence[cur_pos:cur_pos+3]) in params['stop_codons']:
                self.stop_codon = str(extended_sequence[cur_pos:cur_pos+3])
                self.seq_end = cur_pos + 4
            self.frameshifted_seq = frameshifted_seq
            self.frameshifted_translate = frameshifted_translate



def download_gbk_files(entrez_email, entrez_api_key, genome_path):
    '''
    Downloads genbank files to genome_path output dir using accession numbers from input json file
    Parameters:
        entrez_email (string): entrez email
        entrez_api_key (string) - entrez api key
        genome_path (string) - output directory to save genbank files to
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

def find_heptamer(feature, signals, start_pos, stop_pos):
    '''
    Search for heptamer in geneome feature within the given range
    Parameters:
        feature (GenomeFeature obj): object with the sequence to scan
        signals [("String", Int)]: array of tuples ("Heptamer", Score)
        start_pos (int): start position of scan
        stop_pos (int): stop position of scan
    Return:
        None, prints frameshifted sequence + creates and appends Frameshift object detected_frameshifts array
    '''
    if args.verbose:print('Searching from', start_pos, 'to', stop_pos)
    current_pos = start_pos
    for signal in signals:
        while current_pos < stop_pos-6:
            if str(feature.spliced_sequence[current_pos:current_pos+7]) in signal[0]:
                print('Found Downstream Heptamer')
                frameshift = Frameshift(feature, signal[0], 'Downstream', current_pos, 0, stop_pos)
                detected_frameshifts.append(frameshift)
                print_pos = start_pos
                while print_pos <= stop_pos + 1:
                    if print_pos == current_pos:
                        print(colored(feature.spliced_sequence[print_pos:print_pos+3], 'green'), end=' ')
                        print(colored(feature.spliced_sequence[print_pos+3:print_pos+4], 'green'), end=' ')
                        print(colored(feature.spliced_sequence[print_pos+4:print_pos+7], 'green'), end=' ')
                        print_pos += 7
                    else:
                        print(feature.spliced_sequence[print_pos:print_pos+3], end=' ')
                        print_pos += 3
                print()
                break
            current_pos += 3

def find_upstream_frameshift(feature, shift, dstream_limit, stop_codons, signals):
    '''
    Search for +1 upstream frameshift in given feature - this is the cases where we frameshift
    into the annotated gene
    Parameters:
        feature - search for frameshift for this feature
    '''
    if args.verbose: print('\n********** Searching for uptream frameshift **********')

def find_downstream_frameshift(feature, shift, ustream_limit, stop_codons, signals):
    '''
    Search for +1 frameshifts downstream in given feature - this is the case where we frameshift
    out of the annotated gene.
    Parameters:
        feature (GenomeFeature): search for frameshift for this feature
        shift (int): integer value that determines search shift i.e. + 1
        ustream_limit (int):  how far up we should look
        stop_codons ["",""]: array of strings of stop codons
        signals [("String", Int)]: array of tuples ("Heptamer", Score)
    '''
    if args.verbose: print('\n********** Searching for downstream frameshift **********')
    source_start_codon_pos = 0
    source_stop_codon_pos = len(feature.spliced_sequence) - 3 # the stop codon in the source frame

    if args.verbose:
        print('\nSource Frame: ', end='')
        for x in range(0, len(feature.spliced_sequence)//3):
            if x == 0:
                print(colored(feature.spliced_sequence[x*3:x*3+3], 'green'), end=' ')
            elif x == len(feature.spliced_sequence)//3 - 1:
                print(colored(feature.spliced_sequence[x*3:x*3+3], 'red'))
            else:
                print(feature.spliced_sequence[x*3:x*3+3], end=' ')

    ustream_count = 0
    # Frameshift + 1 to destination, start at stop + 1, go upstream until we find the first stop codon
    destination_stop_codon_pos = []
    current_pos = source_stop_codon_pos + shift

    while ustream_count < ustream_limit and current_pos >= 0:
        if feature.spliced_sequence[current_pos:current_pos+3] in stop_codons:
            roi_left = current_pos # found first +1 destination frame stop codon
            destination_stop_codon_pos.append(current_pos)
        current_pos -= 3
        ustream_count += 3
    
    if args.verbose:
        print('\nDestination Frame: ', end='')
        for x in range(0, len(feature.spliced_sequence)//3):
            if x*3+1 in destination_stop_codon_pos:
                print(colored(feature.spliced_sequence[x*3+1:x*3+1+3], 'red'), end=' ')
            else:
                print(feature.spliced_sequence[x*3+1:x*3+1+3], end=' ')        
        print()


    destination_stop_codon_pos = sorted(destination_stop_codon_pos)

    if args.verbose:
        print('\nDestination stop positions: ', end='')
        for stop_pos in destination_stop_codon_pos:
            print(stop_pos, end=', ')
    
    if args.verbose: print('\n\n...Searching for heptamers...')
    # search for signal in roi's
    for roi_count in range(0,len(destination_stop_codon_pos)):
        if roi_count == 0: 
            # source start to destination first stop
            find_heptamer(feature, signals, source_start_codon_pos, destination_stop_codon_pos[roi_count]-1)
            if len(destination_stop_codon_pos) == 1:
                # destination stop to source stop
                find_heptamer(feature, signals, destination_stop_codon_pos[roi_count]-1, source_stop_codon_pos)
                break
            else:
                # destination first stop to next stop
                find_heptamer(feature, signals, destination_stop_codon_pos[roi_count]-1, destination_stop_codon_pos[roi_count+1]-1)
        elif roi_count == len(destination_stop_codon_pos)-1: 
            # destination stop to source stop
            find_heptamer(feature, signals, destination_stop_codon_pos[roi_count]-1, source_stop_codon_pos)
        else:
            # destination stop to next destination stop
            find_heptamer(feature, signals, destination_stop_codon_pos[roi_count]-1, destination_stop_codon_pos[roi_count+1]-1)

def write_to_txt(output_filename):
    '''
    Create and write to text file with frameshift information
    '''
    with open(output_filename+'.txt','w') as outfile:
        for fs in detected_frameshifts:
            if fs.stop_codon != 'None':
                outfile.write('\n\nAccession: ' + fs.genome_feature.accession)
                outfile.write('\nDescription: ' + fs.genome_feature.description)
                outfile.write('\nLocus Tag: ' + fs.genome_feature.locus_tag)
                outfile.write('\nProtein ID: ' + fs.genome_feature.protein_id)
                outfile.write('\nProduct: ' + fs.genome_feature.product)
                outfile.write('\nStrand: ' + str(fs.genome_feature.strand))
                outfile.write('\nCase: ' + str(fs.case))
                outfile.write('\nSignal Found: ' + fs.signal_found)
                outfile.write('\nFrameshift Stop Codon: ' + fs.stop_codon)
                outfile.write('\nOriginal Location: [' + str(fs.genome_feature.location))
                outfile.write('\nOriginal Sequence:\n' + fs.get_original_seq())
                outfile.write('\nOriginal Translate:\n' + str(fs.genome_feature.spliced_sequence.translate()))
                outfile.write('\nFrameshift Location: [' + str(fs.genome_feature.get_true_pos(fs.start_pos)) + 
                ':' + str(fs.genome_feature.get_true_pos(fs.seq_end)) + ']')
                outfile.write('\nFrameshifted Sequence:\n' + str(fs.frameshifted_seq))
                outfile.write('\nFrameshifted Translate:\n' + str(fs.frameshifted_translate))

def write_to_csv(output_filename):
    fields = ['Accession', 'Description', 'Locus Tag', 'Protein ID', 'Product', 'Strand', 'Case', 'Signal', 'Stop Codon', 
    'Original Location', 'Frameshift Location', 'Original Product', 'Frameshift Product', 'Original Sequence', 'Frameshift seq']
    with open(output_filename + '.csv', 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        for fs in detected_frameshifts:
            if fs.stop_codon != 'None':
                csvwriter.writerow([fs.genome_feature.accession, fs.genome_feature.description, fs.genome_feature.locus_tag,
                fs.genome_feature.protein_id, fs.genome_feature.product, str(fs.genome_feature.strand), str(fs.case), fs.signal_found,
                fs.stop_codon, str(fs.genome_feature.location), [str(fs.genome_feature.get_true_pos(fs.start_pos)) + ':' + 
                str(fs.genome_feature.get_true_pos(fs.seq_end))], str(fs.genome_feature.spliced_sequence.translate()),
                str(fs.frameshifted_translate), fs.get_original_seq(), str(fs.frameshifted_seq)])

        
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

    print('Opening input File: ', args.input[0])
    global params
    with open(args.input[0]) as json_file:  
        params = json.load(json_file)
    
    if args.download_GBK_files:
        download_gbk_files(params['entrez_email'], params['entrez_apikey'], params['genome_path'])

    # For each nucleotide record associated with assembly
    for chromid_id in params['assembly_chromids']:
        nfile = open(params['genome_path'] + '/' + chromid_id + '.gb', "r")
        global nucrec
        nucrec = SeqIO.read(nfile, "genbank")
        nfile.close()

        print("Processing: " + nucrec.id)
        # For each feature in the genebank file, create a GenomeFeature object and store in global genome_features array
        for i in range(len(nucrec.features)):
            feat = nucrec.features[i]
            if (feat.type == 'CDS'):
                if (feat.strand == 1):
                    genome_features.append(GenomeFeature(nuc_acc=nucrec.name, nuc_desc=nucrec.description, feature=feat))
    
        # For each GenomeFeature in genome_features, search for upstream and downstream frameshifts
        for feature in genome_features:
            find_downstream_frameshift(feature, params['frame'], params['ustream_limit'], params['stop_codons'], params['signals'])

    write_to_txt(params['outfile_name'])
    write_to_csv(params['outfile_name'])