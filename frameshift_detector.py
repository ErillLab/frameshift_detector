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
import time
import os

detected_frameshifts = [] # Holds all Frameshift objects created after successful heptamer detection
genome_features = [] # Holds all GenomeFeature objects created from features read in from genbank files

class GenomeFeature:
    '''
    Contains information about a specific feature pulled from a Genbank file.
    '''
    def __init__(self, nuc_acc, nuc_desc, curr_feat):
        # Basic feature metadata
        self.accession = nuc_acc
        self.description = nuc_desc
        self.curr_feat = curr_feat
        # If CDS has multiple exons
        if self.curr_feat.location_operator == 'join':
            if args.verbose: print('Feature has multiple exons')
            self.spliced_sequence = self.get_spliced_DNA()[0]
            self.splice_seq_true_positions = self.get_spliced_DNA()[1]
        # If CDS has one exon
        else:
            self.spliced_sequence = nucrec[self.curr_feat.location.start:self.curr_feat.location.end].seq
            self.splice_seq_true_positions = list(range(int(self.curr_feat.location.start), int(self.curr_feat.location.end)+1))
        self.spliced_start_pos = 0
        self.spliced_stop_pos = len(self.spliced_sequence)-3
    
    def set_prev_feat(self, feat):
        '''
        Set the prev feature, used primarily in the upstream case
        Parameters:
            feat (GenomeFeature obj): the previous feature
        Return:
            None. Sets prev_feat member variable
        '''
        self.prev_feat = feat
    
    def set_next_feat(self, feat):
        '''
        Set the next feature, used primarily in the downstream case
        Parameters:
            feat (GenomeFeature obj): the next feature
        Return:
            None. Sets next_feat member variable
        '''
        self.next_feat = feat
        if self.next_feat != None: 
            # Append next feature to the end of the current feature
            self.spliced_sequence += feat.spliced_sequence
            # Append next feature's true positions to current feature's true positions
            self.splice_seq_true_positions.extend(feat.splice_seq_true_positions) 

    def get_next_feat_metadata(self):
        '''
        Get next feature's metadata (locus tag, prot id, and product)
        Return:
            output (string): string of next feature's metadata
        '''
        output = ''
        if self.next_feat != None:
            output += '\nNext Feat Locus Tag: ' + self.next_feat.curr_feat.qualifiers['locus_tag'][0]
            output += '\nNext Feat Protein ID: ' + self.next_feat.curr_feat.qualifiers['protein_id'][0]
            output += '\nNext Feat Product: ' + self.next_feat.curr_feat.qualifiers['product'][0]
        return output

    def get_unspliced_DNA(self):
        '''
        Keep introns in the feature's DNA sequence
        Return:
            unspliced_seq (seq object): unspliced DNA sequence
        '''
        first_exon_start = self.curr_feat.location.parts[0].start # start of first exon
        last_exon_end = self.curr_feat.location.parts[-1].end # end of last exon
        unspliced_seq = nucrec[first_exon_start:last_exon_end].seq
        return unspliced_seq

    def get_spliced_DNA(self):
        '''
        Combines all exons into one sequence for the genome feature
        Return:
            (seq object, [int]): tuple containing the spliced sequence [0] and an int array of 
            the true positions for the nucrec [1]
        '''
        sequence_string = ""
        true_positions = []
        for part in self.curr_feat.location.parts:
            #print(part)
            true_positions.extend(list(range(int(part.start), int(part.end)+1)))
            sequence_string += str(nucrec[part.start:part.end].seq)
        #print(sequence_string)
        #print(true_positions,'\n')
        return (Seq(sequence_string), true_positions)

    def get_true_pos(self, spliced_sequence_pos):
        '''
        Returns the true nucrec position given a spliced sequence position
        Parameters:
            spliced_sequence_pos (int): position in the spliced DNA sequence
        Return:
            true_pos (int): nucrec position
        '''
        return self.splice_seq_true_positions[spliced_sequence_pos]



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

    def get_original_seq(self):
        '''
        Adds space between every 3 base pairs in original sequence
        Return:
            (string): string for the original sequence with spacing between codons
        '''
        seq_string = str(self.genome_feature.spliced_sequence)
        return ' '.join(seq_string[i:i+3] for i in range(0,len(seq_string),3))

    def get_frameshifted_seq(self):
        '''
        Returns a string of the frameshifted sequence with a space between 3 base pairs
        Return:
            (string): string for the frameshifted sequence with spacing between codons
        '''
        frameshifted_seq = ''
        if self.case == 'Downstream':
            cur_pos = self.start_pos
            while cur_pos <= self.stop_pos + 3:
                if cur_pos == self.heptamer_location:
                    frameshifted_seq += str(self.genome_feature.spliced_sequence[cur_pos:cur_pos+3]) + ' '
                    frameshifted_seq += str(self.genome_feature.spliced_sequence[cur_pos+3:cur_pos+4]) + ' '
                    cur_pos += 4
                else:
                    frameshifted_seq += str(self.genome_feature.spliced_sequence[cur_pos:cur_pos+3]) + ' '
                    cur_pos += 3
            return frameshifted_seq



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

def find_heptamer(feature, signals, start_pos, stop_pos, fs_stop_pos):
    '''
    Search for heptamer in geneome feature within the given range
    Parameters:
        feature (GenomeFeature obj): object with the sequence to scan
        signals [("String", Int)]: array of tuples ("Heptamer", Score)
        start_pos (int): start position of scan
        stop_pos (int): stop position of scan
        fs_stop (int): the stop position of the frameshift, will be equal to stop_pos in most cases unless 
        we are near the end of the sequence.
    Return:
        None, prints frameshifted sequence + creates and appends Frameshift object detected_frameshifts array
    '''
    if args.verbose:print('Searching from', start_pos, 'to', stop_pos)
    current_pos = start_pos
    for signal in signals:
        while current_pos < stop_pos-6:
            if str(feature.spliced_sequence[current_pos:current_pos+7]) in signal[0]:
                print('Found Downstream Heptamer')
                frameshift = Frameshift(feature, signal[0], 'Downstream', current_pos, 0, fs_stop_pos)
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
    if args.verbose:
        print('\n********** Searching for downstream frameshift **********')
        print('\nSource Frame: ', end='')
        for x in range(0, len(feature.spliced_sequence)//3):
            if x == feature.spliced_start_pos:
                print(colored(feature.spliced_sequence[x*3:x*3+3], 'green'), end=' ')
            elif x == feature.spliced_stop_pos//3:
                print(colored(feature.spliced_sequence[x*3:x*3+3], 'red'), end=' ')
            else:
                print(feature.spliced_sequence[x*3:x*3+3], end=' ')
        print()
    ustream_count = 0
    # Frameshift + 1 to destination and go downstream from the annotated stop until the first frameshifted stop is reached
    destination_stop_codon_pos = []
    current_pos = feature.spliced_stop_pos + shift
    while feature.spliced_sequence[current_pos:current_pos+3] not in stop_codons and current_pos < len(feature.spliced_sequence):
        current_pos += 3

    downstream_stop_pos = current_pos
    if args.verbose: print('Found stop', feature.spliced_sequence[current_pos:current_pos+3])
    current_pos -=3
    # From the first downstream stop codon, go back upstream and find the positions of all stop codons in the destination
    while ustream_count < ustream_limit and current_pos >= 0:
        if feature.spliced_sequence[current_pos:current_pos+3] in stop_codons:
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

    # Search for the heptamers in the ROIS
    if args.verbose: print('\n\n...Searching for heptamers...')
    for roi_count in range(0,len(destination_stop_codon_pos)):
        if roi_count == 0: 
            # source start to destination first stop
            find_heptamer(feature, signals, feature.spliced_start_pos, destination_stop_codon_pos[roi_count]-1, destination_stop_codon_pos[roi_count]-1)
            if len(destination_stop_codon_pos) == 1:
                # destination stop to source stop
                find_heptamer(feature, signals, destination_stop_codon_pos[roi_count]-1, feature.spliced_stop_pos, downstream_stop_pos)
                break
            else:
                # destination first stop to next stop
                find_heptamer(feature, signals, destination_stop_codon_pos[roi_count]-1, destination_stop_codon_pos[roi_count+1]-1, destination_stop_codon_pos[roi_count+1]-1)
        elif roi_count == len(destination_stop_codon_pos)-1: 
            # destination stop to source stop
            find_heptamer(feature, signals, destination_stop_codon_pos[roi_count]-1, feature.spliced_stop_pos, downstream_stop_pos)
        else:
            # destination stop to next destination stop
            find_heptamer(feature, signals, destination_stop_codon_pos[roi_count]-1, destination_stop_codon_pos[roi_count+1]-1, destination_stop_codon_pos[roi_count+1]-1)

def write_to_txt(output_filename):
    '''
    Create and write to text file with frameshift information
    '''
    with open(output_filename+'.txt','w') as outfile:
        for fs in detected_frameshifts:
            outfile.write('\n\nAccession: ' + fs.genome_feature.accession)
            outfile.write('\nDescription: ' + fs.genome_feature.description)
            outfile.write('\nLocus Tag: ' + fs.genome_feature.curr_feat.qualifiers['locus_tag'][0])
            outfile.write('\nProtein ID: ' + fs.genome_feature.curr_feat.qualifiers['protein_id'][0])
            outfile.write('\nProduct: ' + fs.genome_feature.curr_feat.qualifiers['product'][0])
            outfile.write(fs.genome_feature.get_next_feat_metadata())
            outfile.write('\nStrand: ' + str(fs.genome_feature.curr_feat.strand))
            outfile.write('\nCase: ' + str(fs.case))
            outfile.write('\nSignal Found: ' + fs.signal_found)
            outfile.write('\nOriginal Location: [' + str(fs.genome_feature.curr_feat.location))
            outfile.write('\nOriginal Sequence:\n' + fs.get_original_seq())
            outfile.write('\nFrameshift Location: [' + str(fs.genome_feature.get_true_pos(fs.start_pos)) + 
            ':' + str(fs.genome_feature.get_true_pos(fs.stop_pos)) + ']')
            outfile.write('\nFrameshifted Sequence:\n' + fs.get_frameshifted_seq())


        
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
            if (feat.type == 'CDS') and (feat.strand == 1):
                genome_features.append(GenomeFeature(nuc_acc=nucrec.name, nuc_desc=nucrec.description, curr_feat=feat))
    
        # For each GenomeFeature in genome_features, search for upstream and downstream frameshifts
        for i in range(0, len(genome_features)):
            if i == len(genome_features) - 1:
                genome_features[i].set_next_feat(None)
            else:
                genome_features[i].set_next_feat(genome_features[i+1]) # Set next feature to search in downstream
            find_downstream_frameshift(genome_features[i], params['frame'], params['ustream_limit'], params['stop_codons'], params['signals'])

    write_to_txt(params['outfile_name'])