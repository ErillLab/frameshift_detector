'''
Author: Sara Tahir
Date: 04/15/2022
Purpose: To identify new +1 programmed ribosomal frameshifts through the detection of specific heptameric sequences
Example usage: python3 frameshift_detector.py input.json --verbose 
'''

from calendar import c
from termcolor import colored
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from lib2to3.pgen2.tokenize import generate_tokens
from posixpath import split
import ncbi.datasets
from ete3 import NCBITaxa
import argparse
import json
import csv
import time
import os

detected_frameshifts = [] # Holds all Frameshift objects created after successful heptamer detection
genome_features = [] # Holds all GenomeFeature objects created from features read in from genbank files

# start an api_instance 
api_instance = ncbi.datasets.GenomeApi(ncbi.datasets.ApiClient())
ncbi = NCBITaxa()
#ncbi.update_taxonomy_database()

class GenomeFeature:
    '''
    Contains information about a specific feature pulled from a Genbank file.
    '''
    def __init__(self, nucrec, nuc_acc, nuc_desc, feature, strand):
        # Basic feature metadata
        self.nucrec = nucrec
        self.accession = nuc_acc
        self.description = nuc_desc
        self.feature = feature
        self.type = feature.type
        self.strand = strand
        self.location = feature.location
        self.location_str = self.location
        if strand == '-1':
            self.location_str = ''
            if len(feature.location.parts) > 1:
                     self.location_str += 'join'
            for part in feature.location.parts:
                self.location_str += '[' + str(len(self.nucrec) - part.start) + ',' + str(len(self.nucrec) - part.end) + '](+)'
        self.locus_tag = feature.qualifiers['locus_tag'][0]
        try:
            feature.qualifiers['protein_id'][0]
            self.protein_id = feature.qualifiers['protein_id'][0]
        except:
            self.protein_id = ''

        try:
            feature.qualifiers['product'][0]
            self.product = feature.qualifiers['product'][0]
        except:
            self.product = ''
        self.annotated = False
        # If CDS has multiple exons
        if self.feature.location_operator == 'join':
            if(self.feature.location.parts[1].start - self.feature.location.parts[0].end == 1):
                self.spliced_seq = self.nucrec[self.location.start:self.location.end].seq
                self.spliced_seq_true_pos = list(range(int(self.location.start), int(self.location.end)))
                self.annotated = True
            else:
                self.get_spliced_DNA()
        # If CDS has one exon
        else:
            self.spliced_seq = self.nucrec[self.location.start:self.location.end].seq
            self.spliced_seq_true_pos = list(range(int(self.location.start), int(self.location.end)))
        # Set downstream region (x nucleotides after annotated stop)
        self.get_downstream_region()
        # Set upstream region (x nucleotides before annotated start)
        self.get_upstream_region()

    def splice(self):
        '''
        Combines all exons into one sequence for the genome feature
        Return:
            (seq object, [int]): tuple containing the spliced sequence [0] and an int array of 
            the true positions for the nucrec
        '''
        sequence_string = ""
        for part in self.location.parts:
            #print(part)
            sequence_string += str(self.nucrec[part.start:part.end].seq)
        #print(true_positions,'\n')
        return Seq(sequence_string)

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
            true_positions.extend(list(range(int(part.start), int(part.end))))
            sequence_string += str(self.nucrec[part.start:part.end].seq)
        #print(true_positions,'\n')
        self.spliced_seq = Seq(sequence_string)
        self.spliced_seq_true_pos = true_positions
    
    def get_downstream_region(self):
        '''
        Adds x downstream nucleotides to an array 
        Return:
            (seq object, [int]): tuple containing the downstream sequence [0] and an int array of 
            the true positions for the nucrec [1]
        '''
        last_exon_end = self.location.parts[-1].end # end of last exon
        self.downstream_region_seq = self.nucrec[last_exon_end:last_exon_end+10000].seq
        self.downstream_region_true_pos = list(range(int(last_exon_end), int(last_exon_end)+10000))

    def get_upstream_region(self):
        '''
        Adds x upstream nucleotides to an array 
        Return:
            (seq object, [int]): tuple containing the upstream sequence [0] and an int array of 
            the true positions for the nucrec [1]
        '''
        first_exon_start = self.location.parts[0].start # start of first exon
        if first_exon_start < 10000:
            upstream_seq = self.nucrec[1:first_exon_start].seq
            true_positions = list(range(1, int(first_exon_start)+1))
        else:
            upstream_seq = self.nucrec[first_exon_start-10000:first_exon_start].seq
            true_positions = list(range(first_exon_start-10000, int(first_exon_start)+1))
        self.upstream_region_seq = upstream_seq
        self.upstream_region_true_pos = true_positions

    def get_true_pos_downstream(self, spliced_sequence_pos):
        '''
        Returns nucrec position given a spliced sequence position for downstream case
        Parameters:
            spliced_sequence_pos (int): position in the spliced DNA sequence
        Return:
            true_pos (int): nucrec position
        '''
        # combine spliced sequence and downstream region position arrays, lookup pos from combined seq
        combined_sequence_pos = self.spliced_seq_true_pos + self.downstream_region_true_pos
        if self.strand == '+1':
            return combined_sequence_pos[spliced_sequence_pos]
        else:
            return len(self.nucrec) - combined_sequence_pos[spliced_sequence_pos]
   
    def get_true_pos_upstream(self, spliced_sequence_pos):
        '''
        Returns nucrec position given a spliced sequence position for upstream case
        Parameters:
            spliced_sequence_pos (int): position in the spliced DNA sequence
        Return:
            true_pos (int): nucrec position
        '''
        # combine spliced sequence and downstream region position arrays, lookup pos from combined seq
        combined_sequence_pos = self.upstream_region_true_pos + self.spliced_seq_true_pos
        if self.strand == '+1':
            return combined_sequence_pos[spliced_sequence_pos]
        else:
            return len(self.nucrec) - combined_sequence_pos[spliced_sequence_pos]


class Frameshift:
    '''
    Holds information for a detected frameshift
    '''
    def __init__(self, genome_feature, signal_found, signal_score, fs_case, heptamer_location, start_pos, stop_pos):
        self.genome_feature = genome_feature
        self.signal_found = signal_found
        self.signal_score = signal_score
        self.case = fs_case
        self.heptamer_location = heptamer_location
        self.start_pos = start_pos 
        self.stop_pos = stop_pos
        self.seq_end = 0 
        self.stop_codon = 'None' # Default value of None in-case the frameshift does not have a real stop before the limit
        self.get_frameshift_seq()

    def get_original_seq(self):
        '''
        Adds space after every 3 nucleotides in original sequence
        Return:
            (string): string for the original sequence with spacing between codons
        '''
        seq_string = str(self.genome_feature.spliced_seq)
        return seq_string
        #return ' '.join(seq_string[i:i+3] for i in range(0,len(seq_string),3))

    def get_frameshift_seq(self):
        '''
        Returns a string of the frameshift sequence with a space after every 3 nucleotides also handling frameshift spacing
        Return:
            (string): string for the frameshift sequence with spacing between codons
        '''
        frameshift_seq = ''
        frameshift_translation = ''
        
        if self.case == 'Downstream':
            extended_seq = str(self.genome_feature.spliced_seq) + str(self.genome_feature.downstream_region_seq)
            cur_pos = self.start_pos
        elif self.case == 'Upstream':
            extended_seq = str(feature.upstream_region_seq) + str(feature.spliced_seq)
            cur_pos = self.start_pos
        # Traverse extended seq until hit a stop codon or the downstream limit
        while extended_seq[cur_pos:cur_pos+3] not in params['stop_codons'] and (cur_pos-self.start_pos < params['dstream_limit']):
            # At frameshift postion
            if cur_pos == self.heptamer_location:
                frameshift_seq += str(extended_seq[cur_pos:cur_pos+3])
                frameshift_translation += str(Seq(extended_seq[cur_pos:cur_pos+3]).translate())
                #frameshift_seq += str(extended_seq[cur_pos+3:cur_pos+4])
                #print(str(extended_seq[cur_pos:cur_pos+3]), end=' ')
                #print(str(extended_seq[cur_pos+3:cur_pos+4]), end=' ')
                cur_pos += 4

            # At non-frameshift position
            else:
                frameshift_seq += str(extended_seq[cur_pos:cur_pos+3])
                #print(str(extended_seq[cur_pos:cur_pos+3]), end=' ')
                frameshift_translation += str(Seq(extended_seq[cur_pos:cur_pos+3]).translate())
                #print(extended_seq[cur_pos:cur_pos+3], end='')
                cur_pos += 3
        
        # Append stop codon to frameshift sequence
        frameshift_seq += str(extended_seq[cur_pos:cur_pos+3])

        if str(extended_seq[cur_pos:cur_pos+3]) in params['stop_codons']:
            self.stop_codon = str(extended_seq[cur_pos:cur_pos+3])
            print(str(extended_seq[cur_pos:cur_pos+3]), end=' ')
            self.seq_end = cur_pos + 3
        #print()
        self.frameshift_seq = frameshift_seq
        self.frameshift_translation = frameshift_translation

def download_gbk_files(entrez_email, entrez_api_key, genome_path, chromids):
    '''
    Downloads genbank files to genome_path output dir using accession numbers from input json file
    Parameters:
        entrez_email (string): entrez email
        entrez_api_key (string) - entrez api key
        genome_path (string) - output directory to save genbank files to
    '''
    Entrez.email=entrez_email
    Entrez.api_key=entrez_api_key

    folder_check = os.path.isdir(genome_path)
    if folder_check == False:
        os.makedirs(genome_path)

    dir_contents = os.listdir(genome_path)
    #for each chromid in assembly, download GBK, and save using accession number
    print("Downloading GBK files and saving using accession number")
    for chromid_id in chromids:
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

def find_heptamer(sequence, signals, start_pos, stop_pos, case):
    '''
    Search for heptamer in geneome feature within the given range
    Parameters:
        feature (GenomeFeature obj): object with the sequence to scan
        signals [("String", Int)]: array of tuples ("Heptamer", Score)
        start_pos (int): start position of scan
        stop_pos (int): stop position of scan
    Return:
        None, prints frameshift sequence + creates and appends Frameshift object detected_frameshifts array
    '''
    if args.verbose:print('Searching from', start_pos, 'to', stop_pos)
    #print(sequence[start_pos:stop_pos])
    for signal in signals:
        current_pos = start_pos
        while current_pos < stop_pos-6:
            if str(sequence[current_pos:current_pos+7]) in signal[0]:
                print('Found Heptamer', feature.strand, feature.accession, feature.locus_tag, feature.protein_id, feature.product)
                if case == 'Downstream':
                    frameshift = Frameshift(feature, signal[0], signal[1], case, current_pos, 0, stop_pos)
                elif case == 'Upstream':
                    frameshift = Frameshift(feature, signal[0], signal[1], case, current_pos, start_pos, feature.location.end)
                detected_frameshifts.append(frameshift)
                print_pos = start_pos
                while print_pos <= stop_pos + 1:
                    if print_pos == current_pos:
                        print(colored(sequence[print_pos:print_pos+3], 'green'), end=' ')
                        print(colored(sequence[print_pos+3:print_pos+4], 'green'), end=' ')
                        print(colored(sequence[print_pos+4:print_pos+7], 'green'), end=' ')
                        print_pos += 7
                    else:
                        print(sequence[print_pos:print_pos+3], end=' ')
                        print_pos += 3
                print()
                break   
            #print(sequence[current_pos:current_pos+7]) 
            current_pos += 3

def find_start_codon(sequence, start_pos, stop_pos):
    current_pos = start_pos
    while current_pos < stop_pos:
        if sequence[current_pos:current_pos+3] == 'ATG':
            return current_pos
        current_pos += 3
    return -1

def find_upstream_frameshift(feature, shift, ustream_limit, stop_codons, signals):
    '''
    Search for +1 upstream frameshift in given feature - this is the case where we frameshift
    into the annotated gene
    Parameters:
        feature - search for frameshift for this feature
    '''
    if args.verbose: print('\n********** Searching for uptream frameshift in ', feature.protein_id, feature.locus_tag, '**********')
    
    extended_seq = str(feature.upstream_region_seq) + str(feature.spliced_seq)
    destination_start_codon_pos = len(extended_seq) - len(feature.spliced_seq)
    destination_stop_codon_pos = len(extended_seq) - 3
    #print(extended_seq[destination_start_codon_pos:destination_stop_codon_pos+3])
    #print(feature.nucrec[feature.get_true_pos_upstream(destination_start_codon_pos):feature.get_true_pos_upstream(destination_stop_codon_pos+2)+1].seq)
    #print(feature.nucrec[feature.get_true_pos_upstream(destination_start_codon_pos):feature.get_true_pos_upstream(destination_start_codon_pos)+3].seq)
    #print(feature.nucrec[feature.get_true_pos_upstream(destination_stop_codon_pos):feature.get_true_pos_upstream(destination_stop_codon_pos)+3].seq)
    roi_right = destination_stop_codon_pos
    if args.verbose:
        # Print destination frame with spacing 
        if args.verbose:
            print('\nDestination Frame: ', end='')
            for x in range(0, len(feature.spliced_seq)//3):
                if x == 0:
                    print(colored(feature.spliced_seq[x*3:x*3+3], 'green'), end=' ')
                elif x == len(feature.spliced_seq)//3 - 1:
                    print(colored(feature.spliced_seq[x*3:x*3+3], 'red'))
                else:
                    print(feature.spliced_seq[x*3:x*3+3], end=' ')

    # Find from upstream stop in destination frame
    ustream_count = 0
    current_pos = destination_start_codon_pos
    roi_left = -1
    while ustream_count < ustream_limit and current_pos >= 0:
        if extended_seq[current_pos:current_pos+3] in stop_codons:
            roi_left = current_pos
            break
        current_pos -= 3
        ustream_count += 3
    if roi_left == -1:
        return
    # Print destination after finding the first upstream stop
    if args.verbose:
        print('\nSource Frame Before -1: ', end='')
        print_pos = current_pos
        while print_pos <= (len(extended_seq)):
            if extended_seq[print_pos:print_pos+3] == 'ATG':
                print(colored(extended_seq[print_pos:print_pos+3], 'green'), end=' ')
            elif extended_seq[print_pos:print_pos+3] in stop_codons:
                print(colored(extended_seq[print_pos:print_pos+3], 'red'), end=' ')
            else:
                print(extended_seq[print_pos:print_pos+3], end=' ')
            print_pos += 3
        print()

    # Shift to source frame to -1 from first upstream stop in destination
    if feature.annotated == False:
        current_pos = current_pos - shift
    # Find first ATG upstream in source frame
    roi_left = -1
    while ustream_count < ustream_limit and current_pos >= 0:
        if extended_seq[current_pos:current_pos+3] == 'ATG':
            roi_left = current_pos
            break
        current_pos -= 3
        ustream_count += 3
    if roi_left == -1:
        return
    
    # Print -1 shifted source frame
    if args.verbose:
        print('\nSource Frame After -1: ', end='')
        print_pos = current_pos
        while print_pos <= (len(extended_seq)):
            if extended_seq[print_pos:print_pos+3] == 'ATG':
                print(colored(extended_seq[print_pos:print_pos+3], 'green'), end=' ')
            elif extended_seq[print_pos:print_pos+3] in stop_codons:
                print(colored(extended_seq[print_pos:print_pos+3], 'red'), end=' ')
            else:
                print(extended_seq[print_pos:print_pos+3], end=' ')
            print_pos += 3
        print('\n')

    # Search for stop codons in -1 frame
    source_stop_codon_pos = []
    while current_pos < destination_stop_codon_pos:
        if extended_seq[current_pos:current_pos+3] in stop_codons:
            source_stop_codon_pos.append(current_pos)
            #print(current_pos, extended_seq[current_pos:current_pos+3])
        current_pos += 3

    for roi in range(0,len(source_stop_codon_pos)):
         # if first roi, look for start codon from roi left to source first stop, if found search for heptamer
         # from start codon to end 
        if roi == 0:
            #start_codon_pos = find_start_codon(extended_seq, roi_left-1, source_stop_codon_pos[roi]-1)
            #print(start_codon_pos, roi_left-1, source_stop_codon_pos[roi]-1)
            find_heptamer(extended_seq, signals, roi_left, source_stop_codon_pos[roi], 'Upstream')
        # Else, look for heptamer from destination stop to next destination stop
        if roi == len(source_stop_codon_pos)-1:
            start_codon_pos = find_start_codon(extended_seq, source_stop_codon_pos[roi], destination_stop_codon_pos)
            #print(start_codon_pos, source_stop_codon_pos[roi]-1, destination_stop_codon_pos)
            if start_codon_pos != -1:find_heptamer(extended_seq, signals, start_codon_pos, destination_stop_codon_pos, 'Upstream')
        else:
            start_codon_pos = find_start_codon(extended_seq, source_stop_codon_pos[roi], source_stop_codon_pos[roi+1])
            #print(start_codon_pos, source_stop_codon_pos[roi]-1, source_stop_codon_pos[roi+1]-1)
            if start_codon_pos != -1:find_heptamer(extended_seq, signals, start_codon_pos, source_stop_codon_pos[roi+1], 'Upstream')
            
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
    if args.verbose: print('\n********** Searching for downstream frameshift in ', feature.protein_id, feature.locus_tag, '**********')
    source_start_codon_pos = 0
    source_stop_codon_pos = len(feature.spliced_seq) - 3 # the stop codon in the source frame

    # Print source frame with spacing 
    if args.verbose:
        print('\nSource Frame: ', end='')
        for x in range(0, len(feature.spliced_seq)//3):
            if x == 0:
                print(colored(feature.spliced_seq[x*3:x*3+3], 'green'), end=' ')
            elif x == len(feature.spliced_seq)//3 - 1:
                print(colored(feature.spliced_seq[x*3:x*3+3], 'red'))
            else:
                print(feature.spliced_seq[x*3:x*3+3], end=' ')

    ustream_count = 0
    # Frameshift + 1 to destination, start at stop + 1, go upstream until we find the first stop codon
    destination_stop_codon_pos = []
    if feature.annotated:
        current_pos = source_stop_codon_pos
    else:
        current_pos = source_stop_codon_pos + shift

    while ustream_count < ustream_limit and current_pos >= 0:
        if feature.spliced_seq[current_pos:current_pos+3] in stop_codons:
            roi_left = current_pos # found first +1 destination frame stop codon
            destination_stop_codon_pos.append(current_pos)
        current_pos -= 3
        ustream_count += 3
    
    # Print destination frame with spacing 
    if args.verbose:
        print('\nDestination Frame: ', end='')
        for x in range(0, len(feature.spliced_seq)//3):
            if x*3+1 in destination_stop_codon_pos:
                print(colored(feature.spliced_seq[x*3+1:x*3+1+3], 'red'), end=' ')
            else:
                print(feature.spliced_seq[x*3+1:x*3+1+3], end=' ')        
        print()

    destination_stop_codon_pos = sorted(destination_stop_codon_pos)

    if args.verbose:
        print('\nDestination stop positions: ', end='')
        for stop_pos in destination_stop_codon_pos:
            print(stop_pos, feature.spliced_seq[stop_pos:stop_pos+3], end=', ')
    
    if args.verbose: print('\n\n...Searching for heptamers...')
    # search for heptamer in roi's
    for roi in range(0,len(destination_stop_codon_pos)):
        # if first roi, look for heptamer from source start to destination first stop
        if roi == 0: 
            find_heptamer(feature.spliced_seq, signals, source_start_codon_pos, destination_stop_codon_pos[roi]-1, 'Downstream')
        # If last roi, look for heptamer from destination stop to source stop
        if roi == len(destination_stop_codon_pos)-1: 
            find_heptamer(feature.spliced_seq, signals, destination_stop_codon_pos[roi]-1, source_stop_codon_pos, 'Downstream')
        # Else, look for heptamer from destination stop to next destination stop
        else:
            find_heptamer(feature.spliced_seq, signals, destination_stop_codon_pos[roi]-1, destination_stop_codon_pos[roi+1]-1, 'Downstream')

def write_to_txt(output_filename):
    '''
    Create and write frameshift information to txt file
    '''
    with open(output_filename+'.txt','w') as outfile:
        for fs in detected_frameshifts:
            if fs.stop_codon != 'None':
                outfile.write('\n\nAccession: ' + fs.genome_feature.accession)
                outfile.write('\nDescription: ' + fs.genome_feature.description)
                outfile.write('\nLocus Tag: ' + fs.genome_feature.locus_tag)
                outfile.write('\nProtein ID: ' + fs.genome_feature.protein_id)
                outfile.write('\nKnown: ' + 'Yes' if fs.genome_feature.annotated else 'No')
                outfile.write('\nProduct: ' + fs.genome_feature.product)
                outfile.write('\nStrand: ' + str(fs.genome_feature.strand))
                outfile.write('\nCase: ' + str(fs.case))
                outfile.write('\nSignal Found: ' + fs.signal_found)
                outfile.write('\nSignal Score: ' + str(fs.signal_score))
                outfile.write('\nFrameshift Stop Codon: ' + fs.stop_codon)
                outfile.write('\nAnnotated Gene Location: [' + str(fs.genome_feature.location_str))
                outfile.write('\nAnnotated Gene Sequence:\n' + fs.get_original_seq())
                outfile.write('\nAnnotated Gene Product Length: ' + str(len(fs.genome_feature.splice().translate())))
                outfile.write('\nAnnotated Gene Product:\n' + str(fs.genome_feature.splice().translate()))
                if fs.case == 'Downstream':
                    outfile.write('\nFrameshift Location: [' + str(fs.genome_feature.get_true_pos_downstream(fs.start_pos)).replace('\'','') + 
                ':' + str(fs.genome_feature.get_true_pos_downstream(fs.seq_end)).replace('\'','') + '](+)')
                elif fs.case == 'Upstream':
                    outfile.write('\nFrameshift Location: [' + str(fs.genome_feature.get_true_pos_upstream(fs.start_pos)).replace('\'','') + 
                ':' + str(fs.genome_feature.get_true_pos_upstream(fs.seq_end)+1).replace('\'','') + '](+)')
                outfile.write('\nSpliced Frameshift Sequence:\n' + str(fs.frameshift_seq))
                outfile.write('\nFrameshift Product Length: ' + str(len(fs.frameshift_translation)))
                outfile.write('\nFrameshift Product:\n' + str(fs.frameshift_translation))

def write_to_csv(output_filename):
    '''
    Create and write frameshift information to csv file
    '''
    fields = ['Accession', 'Description', 'Locus Tag', 'Protein ID', 'Known', 'Product', 'Strand', 'Case', 'Signal', 'Signal Score', 'Frameshift Stop Codon', 
    'Annotated Gene Location', 'Frameshift Location', 'Annotated Gene Product Length', 'Frameshift Product Length', 'Annotated Gene Product', 
    'Frameshift Product', 'Spliced Annotated Gene Sequence', 'Spliced Frameshift Sequence']
    with open(output_filename + '.csv', 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(fields)
        for fs in detected_frameshifts:
            if fs.stop_codon != 'None':
                output = []
                output.append(fs.genome_feature.accession)
                output.append(fs.genome_feature.description)
                output.append(fs.genome_feature.locus_tag)
                output.append(fs.genome_feature.protein_id)
                output.append('Yes' if fs.genome_feature.annotated else 'No')
                output.append(fs.genome_feature.product)
                output.append(str(fs.genome_feature.strand))
                output.append(str(fs.case))
                output.append(fs.signal_found)
                output.append(str(fs.signal_score))
                output.append(fs.stop_codon)
                output.append(str(fs.genome_feature.location_str))
                if fs.case == 'Downstream':
                    output.append([str(fs.genome_feature.get_true_pos_downstream(fs.start_pos)).replace('\'','') + ':' + 
                    str(fs.genome_feature.get_true_pos_downstream(fs.seq_end)).replace('\'','')])
                    print(fs.genome_feature.get_true_pos_downstream(fs.start_pos), fs.genome_feature.get_true_pos_downstream(fs.seq_end))
                    print(fs.genome_feature.nucrec[fs.genome_feature.get_true_pos_downstream(fs.start_pos):fs.genome_feature.get_true_pos_downstream(fs.seq_end)].seq)
                elif fs.case == 'Upstream':
                    output.append([str(fs.genome_feature.get_true_pos_upstream(fs.start_pos)+1).replace('\'','') + ':' + 
                    str(fs.genome_feature.get_true_pos_upstream(fs.seq_end)+1).replace('\'','')])
                    print(fs.genome_feature.get_true_pos_upstream(fs.start_pos), fs.genome_feature.get_true_pos_upstream(fs.seq_end)+1)
                    print(fs.genome_feature.nucrec[fs.genome_feature.get_true_pos_upstream(fs.start_pos):fs.genome_feature.get_true_pos_upstream(fs.seq_end)+1].seq)
                output.append(str(len(fs.genome_feature.splice().translate())))
                output.append(str(len(fs.frameshift_translation)))
                output.append(str(fs.genome_feature.splice().translate()))
                output.append(str(fs.frameshift_translation))
                output.append(fs.get_original_seq())
                output.append(str(fs.frameshift_seq))
                
                csvwriter.writerow(output)

def species_to_taxid(species):
    return ncbi.get_name_translator([species])['Saccharomyces cerevisiae S288C'][0]

def read_input_file(input_csv):
    species_heptamers_dict = {}
    with open(input_csv, newline='') as csvfile:
        file_reader = csv.reader(csvfile, delimiter='\n', quotechar='|')
        next(file_reader) # Skip header row
        for row in file_reader:
            split_row = str(row[0]).split(',', 1)
            species = split_row[0]
            heptamers = split_row[1].strip('"').strip('[').strip(']').split('], ')
            heptamers_and_scores = []
            for heptamer in heptamers:
                heptamer_score_pair = []
                heptamer_score = heptamer.strip('[').split(',')
                heptamer_score_pair.append(heptamer_score[0])
                heptamer_score_pair.append(int(heptamer_score[1]))
                heptamers_and_scores.append(heptamer_score_pair)
            species_heptamers_dict[species] = heptamers_and_scores

    return species_heptamers_dict

def generate_jsons(input_csv, path):
    species_heptamers_dict = read_input_file(input_csv)
    for species in species_heptamers_dict.keys():
        species_dict = {}
        species_dict["species_name"] = species
        species_dict["outfile_name"] = species.replace(' ', '_')
        species_dict["signals"] = species_heptamers_dict[species]
        print(species, '\n', species_dict["signals"])
        taxid = species_to_taxid(species)
        genome_summary = api_instance.assembly_descriptors_by_taxon(
            taxon=str(taxid),
            page_size=100,
            filters_assembly_source='refseq')

        print(f"Number of assemblies: {genome_summary.total_count}")
        for assembly in map(lambda d: d.assembly, genome_summary.assemblies):
            if not assembly.annotation_metadata:
                continue
            n_chr = len(assembly.chromosomes) if assembly.assembly_level == 'Chromosome' else None
            
            print('name:', assembly.org.title, '\n',
                'tax id:', assembly.org.tax_id, '\n',
                'assm_level:', assembly.assembly_level, '\n',
                'num_chromosomes:', n_chr, '\n',
                'accession:', assembly.assembly_accession,'\n'
                'chromosomes:')
            if (assembly.assembly_level == 'Complete Genome'):
                chromosomes = []
                for chromosome in assembly.chromosomes:
                    #print('chromosome ', chromosome.name, ': ', chromosome.accession_version)
                    chromosomes.append(chromosome.accession_version)
                species_dict["assembly_chromids"] = chromosomes

        folder_check = os.path.isdir(path)
        if folder_check == False:
            os.makedirs(path)
    
        with open(path + '/' + species_dict["outfile_name"] + '_input.json', 'w', encoding='utf-8') as f:
            json.dump(species_dict, f, ensure_ascii=False, indent=1)

def create_fasta_file(fasta_file_name):
    print('creating fasta file')
    fasta_file = open(fasta_file_name, 'w+')
    for file in os.listdir(params["results_dir"]):
        if(file[-3:] == 'csv'):
            print(file)
            with open(params["results_dir"] + '/' + file, newline='') as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    fasta_row = '>' + file.replace('_', ' ')[:-4] + ' | ' + row['Locus Tag'] + ' | ' + row['Protein ID'] + ' | ' + row['Product'] + '\n' + row['Frameshift Product'] + '\n'
                    fasta_file.write(fasta_row)
    fasta_file.close()

def makeblastdb(fasta_file_name):
    output_file = 'test_blast.fasta'
    cmd = 'makeblastdb -in {input_file} -out {output_file} -dbtype prot'
    os.system(cmd.format(input_file=fasta_file_name, output_file=output_file))

if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(description='Find some frameshifts')
    # Require input json file
    parser.add_argument('input', nargs='+', help='an input json file with frameshift search parameters')
    # Optional --verbose flag for printing debug statements
    parser.add_argument('--verbose', action='store_true', help='print helpful things')
    # Optional --blast_only flag to skip frameshift detection and only find frameshift conservation
    parser.add_argument('--find_conservation', action='store_true', help='skip frameshift detection and only find frameshift conservation')
    parser.set_defaults(verbose=False, find_conservation=False)

    args = parser.parse_args()

    print('Opening input File: ', args.input[0])
    global params
    with open(args.input[0]) as json_file:  
        params = json.load(json_file)

    if not args.find_conservation:
        generate_jsons(params['csv_input_file'], params['fs_inputs_path'])

        for input_file in os.listdir(params['fs_inputs_path']):
            print(input_file)
            with open(params['fs_inputs_path'] + '/' +input_file) as json_file:  
                species_params = json.load(json_file)
            download_gbk_files(params['entrez_email'], params['entrez_apikey'], params['genome_path'] + '/' + species_params['outfile_name'], species_params['assembly_chromids'])
            #download_gbk_files(params['entrez_email'], params['entrez_apikey'], params['genome_path'], species_params['assembly_chromids'])
            # For each nucleotide record associated with assembly
            for chromid_id in species_params['assembly_chromids']:
                nfile = open(params['genome_path'] + '/' + species_params['outfile_name'] +'/'+chromid_id + '.gb', "r")
                #nfile = open(params['genome_path'] + '/' + chromid_id + '.gb', "r")
                nucrec = SeqIO.read(nfile, "genbank")
                nfile.close()
                
                nucrec_id = nucrec.id
                nucrec_desc = nucrec.description
                print("Processing: " + nucrec.id)
                # For each feature in the genebank file, create a GenomeFeature object and store in global genome_features array
                for i in range(len(nucrec.features)):
                    feat = nucrec.features[i]
                    if (feat.type == 'CDS'):
                        if (feat.strand == 1):
                            genome_features.append(GenomeFeature(nucrec=nucrec, nuc_acc=nucrec_id, nuc_desc=nucrec_desc, feature=feat, strand='+1'))
                
                # Same loop for the reverse strand
                nucrec = nucrec.reverse_complement()
                for i in range(len(nucrec.features)):
                    feat = nucrec.features[i]
                    if (feat.type == 'CDS'):
                        if (feat.strand == 1):
                            genome_features.append(GenomeFeature(nucrec=nucrec, nuc_acc=nucrec_id, nuc_desc=nucrec_desc, feature=feat, strand='-1'))

            # For each GenomeFeature in genome_features, search for upstream and downstream frameshifts
            for feature in genome_features:
                find_upstream_frameshift(feature, params['frame'], params['ustream_limit'], params['stop_codons'], species_params['signals'])
                find_downstream_frameshift(feature, params['frame'], params['ustream_limit'], params['stop_codons'], species_params['signals'])
            print('Writing results')

            folder_check = os.path.isdir(params["results_dir"] )
            if folder_check == False:
                os.makedirs(params["results_dir"] )

            write_to_txt(params["results_dir"] + '/' + species_params['outfile_name'])
            write_to_csv(params["results_dir"] + '/' + species_params['outfile_name'])
    
    create_fasta_file('frameshifts.fasta') 
    makeblastdb('frameshifts.fasta')