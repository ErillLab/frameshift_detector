'''
Author: Sara Tahir
Date: 04/15/2022
Purpose: To identify +1 programmed ribosomal frameshifts through the detection of specific 
heptameric sequences, and show frameshift conservation
Example usage: python3 frameshift_detector.py input.json --verbose 
'''

from calendar import c
from ctypes import alignment
from termcolor import colored
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from lib2to3.pgen2.tokenize import generate_tokens
from posixpath import split
import ncbi.datasets
from ete3 import NCBITaxa
import argparse
import json
import csv
import time
import os
import shutil
import zipfile

detected_frameshifts = {} # Holds all Frameshift objects created after successful heptamer detection
genome_features = [] # Holds all GenomeFeature objects created from features read in from genbank files

# start an api_instance for ncbi datasets
api_instance = ncbi.datasets.GenomeApi(ncbi.datasets.ApiClient())
ncbi = NCBITaxa()
#ncbi.update_taxonomy_database()

class GenomeFeature:
    '''
    Contains information about a specific feature pulled from a Genbank file.
    '''
    def __init__(self, nucrec, nuc_acc, nuc_desc, feature, strand, species):
        # Basic feature metadata
        self.species = species
        self.nucrec = nucrec
        self.accession = nuc_acc
        self.description = nuc_desc
        self.feature = feature
        self.type = feature.type
        self.strand = strand
        self.location = feature.location
        self.location_str = str(self.location).replace("(+)", "")
        # format -1 strand location to match +1 strand location format when there is a spliced sequence
        if strand == '-1':
            self.location_str = ''
            if len(feature.location.parts) > 1:
                     self.location_str += 'join{'
            for part in feature.location.parts:
                self.location_str += '[' + str(len(self.nucrec) - part.start) + ',' + str(len(self.nucrec) - part.end) + ']'
            self.location_str += '}'
        try:
            feature.qualifiers['locus_tag'][0]
            self.locus_tag = feature.qualifiers['locus_tag'][0]
        except:
            self.locus_tag = ''

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
            # If there is only one nucleotide in between two exons, mark this as an annotated frameshift
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
            None, sets spliced_seq and spliced_seq_true_pos member vars
        '''
        sequence_string = ""
        true_positions = []
        for part in self.location.parts:
            true_positions.extend(list(range(int(part.start), int(part.end))))
            sequence_string += str(self.nucrec[part.start:part.end].seq)
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
            upstream_seq = self.nucrec[0:first_exon_start].seq
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
            extended_seq = str(self.genome_feature.upstream_region_seq) + str(self.genome_feature.spliced_seq)
            cur_pos = self.start_pos
            print(self.start_pos)
        # Traverse extended seq until hit a stop codon or the downstream limit
        while extended_seq[cur_pos:cur_pos+3] not in ["TAG","TAA","TGA"] and (cur_pos-self.start_pos < 10000):
            # At frameshift postion
            if cur_pos == self.heptamer_location:
                frameshift_seq += str(extended_seq[cur_pos:cur_pos+3])
                #frameshift_translation += str(Seq(extended_seq[cur_pos:cur_pos+3]).translate())
                #frameshift_seq += str(extended_seq[cur_pos+4:cur_pos+7])
                #frameshift_translation += str(Seq(extended_seq[cur_pos+4:cur_pos+7]).translate())
                #print(str(extended_seq[cur_pos:cur_pos+3]), end=' ')
                #print(str(extended_seq[cur_pos+3:cur_pos+4]), end=' ')
                cur_pos += 4

            # At non-frameshift position
            else:
                frameshift_seq += str(extended_seq[cur_pos:cur_pos+3])
                #print(str(extended_seq[cur_pos:cur_pos+3]), end=' ')
                #frameshift_translation += str(Seq(extended_seq[cur_pos:cur_pos+3]).translate())
                #print(extended_seq[cur_pos:cur_pos+3], end='')
                cur_pos += 3
        
        # Append stop codon to frameshift sequence
        frameshift_seq += str(extended_seq[cur_pos:cur_pos+3])
        #frameshift_translation += str(Seq(extended_seq[cur_pos:cur_pos+3]).translate())
        if str(extended_seq[cur_pos:cur_pos+3]) in ["TAG","TAA","TGA"]:
            self.stop_codon = str(extended_seq[cur_pos:cur_pos+3])
            self.seq_end = cur_pos + 3

        self.frameshift_seq = frameshift_seq
        self.frameshift_translation = str(Seq(frameshift_seq).translate())

def download_gbk_files(data_path, assembly_accessions):
    '''
    Downloads genbank files to genome_path output dir using accession numbers from input json file
    Parameters:
        entrez_email (string): entrez email
        entrez_api_key (string) - entrez api key
        genome_path (string) - output directory to save genbank files to
    '''
    chromosomes = ['']
    exclude_sequence = False
    include_annotation_type = ['GENOME_GB']

    api_response = api_instance.download_assembly_package(
        assembly_accessions,
        chromosomes=chromosomes,
        exclude_sequence=exclude_sequence,
        include_annotation_type=include_annotation_type,
        # Because we are streaming back the results to disk, 
        # we should defer reading/decoding the response
        _preload_content=False
    )

    with open('genome_data.zip', 'wb') as f:
        f.write(api_response.data)

    with zipfile.ZipFile('genome_data.zip', 'r') as zip_ref:
        zip_ref.extractall(data_path)

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
        print(chromid_id)
        if chromid_id is not None:
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
    '''
def find_heptamer(feature, sequence, signals, start_pos, stop_pos, case, args):
    '''
    Search for heptamer in geneome feature within the given range
    Parameters:
        feature (GenomeFeature obj): object with the sequence to scan
        signals [("String", Int)]: array of tuples ("Heptamer", Score)
        start_pos (int): start position of scan
        stop_pos (int): stop position of scan
    Return:
        frameshifts [Frameshift]: array of Frameshift objects
    '''
    if args.verbose:print('Searching from', start_pos, 'to', stop_pos)
    frameshifts = []
    #print(sequence[start_pos:stop_pos])
    for signal in signals:
        current_pos = start_pos
        while current_pos < stop_pos and (stop_pos - start_pos >=7):
            #print(sequence[current_pos:current_pos+3], end=" ")
            if str(sequence[current_pos:current_pos+7]) in signal[0]:
                print('Found Heptamer', case, feature.strand, feature.accession, feature.locus_tag, feature.protein_id, feature.product)
                if case == 'Downstream':
                    #__init__(self, genome_feature, signal_found, signal_score, fs_case, heptamer_location, start_pos, stop_pos):
                    frameshift = Frameshift(feature, signal[0], signal[1], case, current_pos, 0, stop_pos)
                elif case == 'Upstream':
                    if (sequence[start_pos:start_pos+3]) != 'ATG':
                        start_pos = find_codon('ATG', sequence, start_pos, 0)
                    frameshift = Frameshift(feature, signal[0], signal[1], case, current_pos, start_pos, feature.location.end)
                #detected_frameshifts[feature.species].append(frameshift)
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
                frameshifts.append(frameshift)
                break   
            #print(sequence[current_pos:current_pos+7]) 
            current_pos += 3
    return frameshifts

def find_start_codon(sequence, start_pos, stop_pos):
    current_pos = start_pos
    while current_pos < stop_pos:
        if sequence[current_pos:current_pos+3] == 'ATG':
            return current_pos
        current_pos += 3
    return -1

def find_codon(codon, sequence, start_pos, stop_pos):
    direction = 0
    if (stop_pos - start_pos) > 0:
        direction = 1
    else:
        direction = -1
    for i in range(start_pos, stop_pos, direction * 3):
        #print(sequence[i:i+3], end=" ")
        if sequence[i:i+3] == codon:
            return i
    return -1


def find_upstream_frameshift(feature, shift, ustream_limit, stop_codons, signals, args):
    '''
    Search for +1 upstream frameshift in given feature - this is the case where we frameshift
    into the annotated gene
    Parameters:
        feature - search for frameshift for this feature
    '''
    if args.verbose: print('\n********** Searching for uptream frameshift in ', feature.protein_id, feature.locus_tag, '**********')
    
    frameshifts = []

    extended_seq = str(feature.upstream_region_seq) + str(feature.spliced_seq)
    destination_start_codon_pos = len(extended_seq) - len(feature.spliced_seq)
    destination_stop_codon_pos = len(extended_seq) - 4
    #print(extended_seq[destination_start_codon_pos:destination_stop_codon_pos+3])
    #print(feature.nucrec[feature.get_true_pos_upstream(destination_start_codon_pos):feature.get_true_pos_upstream(destination_stop_codon_pos+2)+1].seq)
    #print(feature.nucrec[feature.get_true_pos_upstream(destination_start_codon_pos):feature.get_true_pos_upstream(destination_start_codon_pos)+3].seq)
    #print(feature.nucrec[feature.get_true_pos_upstream(destination_stop_codon_pos):feature.get_true_pos_upstream(destination_stop_codon_pos)+3].seq)
    roi_right = destination_stop_codon_pos
    if args.verbose:
        # Print destination frame with spacing 
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
    while ustream_count < ustream_limit and current_pos >= -1:
        if extended_seq[current_pos:current_pos+3] in stop_codons:
            roi_left = current_pos
            break
        if current_pos == -1:
            roi_left = 0
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
    #if feature.annotated == False:
    current_pos = current_pos - shift
    first_ustream_stop = current_pos
    # Find first ATG upstream in source frame
    first_ustream_start = -1
    while ustream_count < ustream_limit and current_pos >= 0:
        if extended_seq[current_pos:current_pos+3] == 'ATG':
            first_ustream_start = current_pos
            break
        current_pos -= 3
        ustream_count += 3
    if first_ustream_start == -1:
        return frameshifts
    
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
        #print(extended_seq[current_pos:current_pos+3])
        if extended_seq[current_pos:current_pos+3] in stop_codons:
            source_stop_codon_pos.append(current_pos)
            #print(current_pos, extended_seq[current_pos:current_pos+3])
        current_pos += 3

    for roi in range(0,len(source_stop_codon_pos)):
         # if first roi, look for start codon from roi left to source first stop, if found search for heptamer
         # from start codon to end 
        if roi == 0:
            #start_codon_pos = find_start_codon(extended_seq, roi_left, source_stop_codon_pos[roi])
            #print(start_codon_pos, roi_left-1, source_stop_codon_pos[roi]-1)
            fs = find_heptamer(feature, extended_seq, signals, first_ustream_stop, source_stop_codon_pos[roi], 'Upstream', args)
            frameshifts = frameshifts + fs
        # Else, look for heptamer from destination stop to next destination stop
        if roi == len(source_stop_codon_pos)-1:
            start_codon_pos = find_start_codon(extended_seq, source_stop_codon_pos[roi], destination_stop_codon_pos)
            if start_codon_pos < destination_start_codon_pos:
                #print(start_codon_pos, source_stop_codon_pos[roi]-1, destination_stop_codon_pos)
                if start_codon_pos != -1:
                    fs = find_heptamer(feature, extended_seq, signals, start_codon_pos, destination_stop_codon_pos, 'Upstream', args)
                    frameshifts = frameshifts + fs
        else:
            start_codon_pos = find_start_codon(extended_seq, source_stop_codon_pos[roi], source_stop_codon_pos[roi+1])
            #print(start_codon_pos, source_stop_codon_pos[roi]-1, source_stop_codon_pos[roi+1]-1)
            if start_codon_pos < destination_start_codon_pos:
                if start_codon_pos != -1:
                    fs = find_heptamer(feature, extended_seq, signals, start_codon_pos, source_stop_codon_pos[roi+1], 'Upstream', args)
                    frameshifts = frameshifts + fs
    
    return frameshifts

def find_stop_codon(feature, stop_codons):
    '''
    Search for stop codon in given feature
    Parameters:
        feature (GenomeFeature): search for stop codon for this feature
    Return:
        stop_codon_pos: stop codon position
    '''
    for x in range(0, len(feature.spliced_seq)//3):
        if feature.spliced_seq[x*3:x*3+3] in stop_codons:
            return x*3
    return -1
            
def find_downstream_frameshift(feature, shift, ustream_limit, stop_codons, signals, args):
    '''
    Search for +1 frameshifts downstream in given feature - this is the case where we frameshift
    out of the annotated gene.
    Parameters:
        feature (GenomeFeature): search for frameshift for this feature
        shift (int): integer value that determines search shift i.e. + 1
        ustream_limit (int):  how far up we should look
        stop_codons ["",""]: array of strings of stop codons
        signals [("String", Int)]: array of tuples ("Heptamer", Score)
    Return:
        frameshifts [Frameshift]: array of Frameshift objects
    '''
    frameshifts = []
    if args.verbose: print('\n********** Searching for downstream frameshift in ', feature.protein_id, feature.locus_tag, '**********')
    source_start_codon_pos = 0
    #source_stop_codon_pos = len(feature.spliced_seq) - 3 # the stop codon in the source frame # Do not assume annotated stop
    source_stop_codon_pos = find_stop_codon(feature, stop_codons)
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
        #print(feature.spliced_seq[current_pos:current_pos+3])
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
            fs =  find_heptamer(feature, feature.spliced_seq, signals, source_start_codon_pos, destination_stop_codon_pos[roi]-1, 'Downstream', args)
            frameshifts = frameshifts + fs
        # If last roi, look for heptamer from destination stop to source stop
        if roi == len(destination_stop_codon_pos)-1: 
            fs = find_heptamer(feature, feature.spliced_seq, signals, destination_stop_codon_pos[roi]-1, source_stop_codon_pos, 'Downstream', args)
            frameshifts = frameshifts + fs
        # Else, look for heptamer from destination stop to next destination stop
        else:
            fs = find_heptamer(feature, feature.spliced_seq, signals, destination_stop_codon_pos[roi]-1, destination_stop_codon_pos[roi+1]-1, 'Downstream', args)
            frameshifts = frameshifts + fs
    
    return frameshifts

def write_to_txt(output_filename, species):
    '''
    Create and write frameshift information to txt file
    '''
    with open(output_filename+'.txt','w') as outfile:
        for fs in detected_frameshifts[species]:
            if fs.stop_codon != 'None':
                outfile.write('\n\nSpecies: ' + fs.genome_feature.species)
                outfile.write('\nAccession: ' + fs.genome_feature.accession)
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
                ':' + str(fs.genome_feature.get_true_pos_downstream(fs.heptamer_location-1)).replace('\'','') + ',[' + str(fs.genome_feature.get_true_pos_downstream(fs.heptamer_location)).replace('\'','') + ':' + str(fs.genome_feature.get_true_pos_downstream(fs.seq_end)).replace('\'','') + '](+)')
                elif fs.case == 'Upstream':
                    outfile.write('\nFrameshift Location: [' + str(fs.genome_feature.get_true_pos_upstream(fs.start_pos)).replace('\'','') + 
                ':' + str(fs.genome_feature.get_true_pos_upstream(fs.seq_end)+1).replace('\'','') + '](+)')
                outfile.write('\nSpliced Frameshift Sequence:\n' + str(fs.frameshift_seq))
                outfile.write('\nFrameshift Product Length: ' + str(len(fs.frameshift_translation)))
                outfile.write('\nFrameshift Product:\n' + str(fs.frameshift_translation))

def write_to_csv(output_filename, species):
    '''
    Create and write frameshift information to csv file
    '''
    print('writing')
    with open(output_filename + '.csv', 'a', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        for fs in detected_frameshifts[species]:
            if fs.stop_codon != 'None':
                output = []
                output.append(fs.genome_feature.species)
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
                    output.append('join{['+str(fs.genome_feature.get_true_pos_downstream(fs.start_pos)).replace('\'','') + ':' + str(fs.genome_feature.get_true_pos_downstream(fs.heptamer_location+3)).replace('\'','') + '],[' +  str(fs.genome_feature.get_true_pos_downstream(fs.heptamer_location+4)).replace('\'','')  + ':' + str(fs.genome_feature.get_true_pos_downstream(fs.seq_end)).replace('\'','') + ']}')
                    #print(fs.genome_feature.get_true_pos_downstream(fs.start_pos), fs.genome_feature.get_true_pos_downstream(fs.seq_end))
                    #print(fs.genome_feature.nucrec[fs.genome_feature.get_true_pos_downstream(fs.start_pos):fs.genome_feature.get_true_pos_downstream(fs.seq_end)].seq)
                elif fs.case == 'Upstream':
                    output.append('join{['+str(fs.genome_feature.get_true_pos_upstream(fs.start_pos)).replace('\'','') + ':' + str(fs.genome_feature.get_true_pos_upstream(fs.heptamer_location+3)).replace('\'','') + '],[' +  str(fs.genome_feature.get_true_pos_upstream(fs.heptamer_location+4)).replace('\'','')  + ':' + str(fs.genome_feature.get_true_pos_upstream(fs.seq_end)).replace('\'','') + ']}')
                    #print(fs.genome_feature.get_true_pos_upstream(fs.start_pos), fs.genome_feature.get_true_pos_upstream(fs.seq_end)+1)
                    #print(fs.genome_feature.nucrec[fs.genome_feature.get_true_pos_upstream(fs.start_pos):fs.genome_feature.get_true_pos_upstream(fs.seq_end)+1].seq)
                output.append(str(len(fs.genome_feature.splice().translate())))
                output.append(str(len(fs.frameshift_translation)))
                output.append(str(len(fs.frameshift_translation) - len(fs.genome_feature.splice().translate())))
                output.append(str(fs.genome_feature.splice().translate()))
                output.append(str(fs.frameshift_translation))
                output.append(fs.get_original_seq())
                output.append(str(fs.frameshift_seq))
                output.append("None")
                
                csvwriter.writerow(output)

def species_to_taxid(species):
    return ncbi.get_name_translator([species])[species][0]

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
    folder_check = os.path.isdir(path)
    if folder_check == True:
        shutil.rmtree(path)
    os.makedirs(path)

    species_heptamers_dict = read_input_file(input_csv)
    for species in species_heptamers_dict.keys():
        species_dict = {}
        species_dict["species_name"] = species
        species_dict["outfile_name"] = species.replace(' ', '_')
        species_dict["signals"] = species_heptamers_dict[species]
        taxid = species_to_taxid(species)
        genome_summary = api_instance.assembly_descriptors_by_taxon(
            taxon=str(taxid),
            page_size=100,
            filters_assembly_source='refseq')
        #print(genome_summary)
        print('Species:', species)
        #print(f"Number of assemblies: {genome_summary.total_count}")
        
        if genome_summary.total_count != None:
            for assembly in map(lambda d: d.assembly, genome_summary.assemblies):
                if not assembly.annotation_metadata:
                    continue
                # TO-DO Apply reference genome filter
                #if (assembly.assembly_level == 'Complete Genome'):
                n_chr = len(assembly.chromosomes) 
                if assembly.assembly_level == 'Chromosome' or assembly.assembly_level == 'Complete Genome':
                    print('assembly: ',assembly.assembly_accession)
                    species_dict["assembly_accession"] = assembly.assembly_accession
                    species_dict["assembly_level"] = assembly.assembly_level
                    chromosomes = []
                    for chromosome in assembly.chromosomes:
                        print('chromosome ', chromosome.name, ': ', chromosome.accession_version)
                        chromosomes.append(chromosome.accession_version)
                    species_dict["assembly_chromids"] = chromosomes
                
                    with open(path + '/' + species_dict["outfile_name"] + '_input.json', 'w', encoding='utf-8') as f:
                        json.dump(species_dict, f, ensure_ascii=False, indent=1)
                else:
                    None

def create_fasta_file_for_blast_db(fasta_file_name):
    if (os.path.exists(fasta_file_name)):
        os.remove(fasta_file_name)

    print('Creating fasta file for blast db')
    fasta_file = open(fasta_file_name, 'w+')

    # Open frameshift csv file(s) and add to fasta file for blast db
    for file in os.listdir(params["results_dir"]):
        if(file[-3:] == 'csv'):
            with open(params["results_dir"] + '/' + file, newline='') as csvfile:
                frameshift_csv = csv.DictReader(csvfile)
                for row in frameshift_csv:
                    fasta_row = '>' + row['Species'] + ' | ' + row['Locus Tag'] + ' | ' + row['Protein ID'] + ' | ' + row['Product'] + ' | ' + row['Frameshift Product Length'] + '\n' + row['Frameshift Product'] + '\n'
                    fasta_file.write(fasta_row)
    fasta_file.close()

def create_blast_query_file(csv_row, species):
    fasta_file_name = 'query_files/' + csv_row['Protein ID'] + '.fasta'
    fasta_file = open( fasta_file_name, 'w+')
    fasta_row = '>' + species + ' | ' + csv_row['Locus Tag'] + ' | ' + csv_row['Protein ID'] + ' | ' + csv_row['Product'] + ' | ' + csv_row['Frameshift Product Length'] + '\n' + csv_row['Frameshift Product'] + '\n'
    fasta_file.write(fasta_row)
    fasta_file.close()
    return fasta_file_name

def makeblastdb(fasta_file_name):
    blast_db = './db/blast_db'
    cmd = 'makeblastdb -in {input} -out {out} -dbtype prot'
    os.system(cmd.format(input=fasta_file_name, out=blast_db))

def append_ortho_group(protein_id, group_num):
    global frameshift_list
    for row in frameshift_list:
        if row['Protein ID'] == protein_id:
            row['Orthology Group'] = group_num

def blast_search():
    print('Beginning blast search')
    cmd = 'blastp -query {query} -db {db} -evalue {e} -out {out} -outfmt 5'
    blast_db = './db/blast_db'
    
    folder_check = os.path.isdir('query_files')
    if folder_check == True:
        shutil.rmtree('query_files')
    os.makedirs('query_files')

    folder_check = os.path.isdir('blast_output')
    if folder_check == True:
        shutil.rmtree('blast_output')
    os.makedirs('blast_output')
    
    e_val = 10E-10

    ortho_group = 1
    for file in os.listdir(params["results_dir"]):
        if(file == 'frameshifts.csv'): # only open master frameshift file
            print(file)
            with open(params["results_dir"] + '/' + file, newline='') as fs_csv_file:
                frameshift_csv = csv.DictReader(fs_csv_file)
                global frameshift_list
                frameshift_list = list(frameshift_csv)
                for row in frameshift_list:
                    if row['Orthology Group'] == "None":
                        append_ortho_group(row['Protein ID'], ortho_group)
                        output_file = 'blast_output/' + row['Protein ID'] + '.xml'
                        query_fasta = create_blast_query_file(row, file.replace('_', ' ')[:-4])
                        os.system(cmd.format(query=query_fasta, db=blast_db, e=e_val ,out=output_file))
                        parse_blast_output(output_file, ortho_group)
                        ortho_group += 1
                
                # Write results
                print('Updating fs conservation csv')
                with open(params["results_dir"] + '/fs_conservation.csv', 'w', newline='') as conservation_csv_file:
                    fields = ['Species','Accession', 'Description', 'Locus Tag', 'Protein ID', 'Known', 'Product', 'Strand', 'Case', 'Signal', 'Signal Score', 'Frameshift Stop Codon', 
        'Annotated Gene Location', 'Frameshift Location', 'Annotated Gene Product Length', 'Frameshift Product Length', 'Product Length Diff', 'Annotated Gene Product', 
        'Frameshift Product', 'Spliced Annotated Gene Sequence', 'Spliced Frameshift Sequence', 'Orthology Group']
                    writer = csv.DictWriter(conservation_csv_file, fieldnames=fields)
                    writer.writeheader()
                    for row in frameshift_list:
                        writer.writerow(row)
                        

def parse_blast_output(output_file, ortho_group):
    with open(output_file) as blast_output:
        blast_records = list(NCBIXML.parse(blast_output))
        for blast_record in blast_records:
            print('\n*************************************************')     
            print(blast_record.query)
            input_seq_len = blast_record.query.split(" | ")[4]
            print("Input sequence length: ", input_seq_len)
            print("Alignments:")
            for alignment in blast_record.alignments:
                coverage = float(alignment.hsps[0].query_end - alignment.hsps[0].query_start + 1) / float(input_seq_len)
                print('coverage: ', coverage)
                for hsp in alignment.hsps:
                    if hsp.expect < params['blast_e_val_threshold'] and coverage > params['blast_coverage_threshold']:
                        protein_id = alignment.title.split(' | ')[2]
                        print(protein_id, alignment.title.split(' | ')[3])
                        append_ortho_group(protein_id, ortho_group)
            print('*************************************************\n')

if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(description='Find some frameshifts')
    # Require input json file
    parser.add_argument('input', nargs='+', help='an input json file with frameshift search parameters')
    # Optional --verbose flag for printing debug statements
    parser.add_argument('--verbose', action='store_true', help='print helpful things')
    # Optional --blast_only flag to skip frameshift detection and only find frameshift conservation
    parser.add_argument('--find_conservation', action='store_true', help='skip frameshift detection and only find frameshift conservation')
    # Optional --protein_id flag to detect frameshift for only a specified protein
    parser.add_argument('--protein_id', help='detect frameshift for only a specified protein')

    parser.set_defaults(verbose=False, find_conservation=False)

    args = parser.parse_args()

    print('Opening input File: ', args.input[0])
    print(args.protein_id)
    global params
    with open(args.input[0]) as json_file:  
        params = json.load(json_file)

    if not args.find_conservation:
        if args.protein_id == None:
            folder_check = os.path.isdir(params["results_dir"] )
            if folder_check == True:
                shutil.rmtree(params["results_dir"])
            os.makedirs(params["results_dir"])

            fields = ['Species','Accession', 'Description', 'Locus Tag', 'Protein ID', 'Known', 'Product', 'Strand', 'Case', 'Signal', 'Signal Score', 'Frameshift Stop Codon', 
        'Annotated Gene Location', 'Frameshift Location', 'Annotated Gene Product Length', 'Frameshift Product Length', 'Product Length Diff', 'Annotated Gene Product', 
        'Frameshift Product', 'Spliced Annotated Gene Sequence', 'Spliced Frameshift Sequence', 'Orthology Group']
            with open(params["results_dir"] + '/frameshifts.csv', 'a', newline='') as csvfile:
                csvwriter = csv.writer(csvfile)
                csvwriter.writerow(fields)
            
            generate_jsons(params['csv_input_file'], params['fs_inputs_path'])
            assembly_accessions = []
            for input_file in os.listdir(params['fs_inputs_path']):
                print(input_file)
                with open(params['fs_inputs_path'] + '/' +input_file) as json_file:  
                    species_params = json.load(json_file)
                    assembly_accessions.append(species_params['assembly_accession'])
            download_gbk_files(params['genome_path'] + '/', assembly_accessions)
            
        data_path = params['genome_path'] + '/ncbi_dataset/data/'
        for input_file in os.listdir(params['fs_inputs_path']):
            print(input_file)
            with open(params['fs_inputs_path'] + '/' +input_file) as json_file:  
                species_params = json.load(json_file)
                detected_frameshifts[species_params['species_name']] = []
                nucrecs = SeqIO.parse(data_path + species_params['assembly_accession'] + '/genomic.gbff', "genbank")
                for nucrec in nucrecs:
                    nucrec_id = nucrec.id
                    print('Processing: ', nucrec_id)
                    nucrec_desc = nucrec.description
                    #print("Processing: " + nucrec.id)
                    # For each feature in the genebank file, create a GenomeFeature object and store in global genome_features array
                    for i in range(len(nucrec.features)):
                        feat = nucrec.features[i]
                        if (feat.type == 'CDS'):
                            if (feat.strand == 1):
                                if args.protein_id == None:
                                    #genome_features.append(GenomeFeature(nucrec=nucrec, nuc_acc=nucrec_id, nuc_desc=nucrec_desc, feature=feat, strand='+1', species=species_params['species_name']))
                                    feature = GenomeFeature(nucrec=nucrec, nuc_acc=nucrec_id, nuc_desc=nucrec_desc, feature=feat, strand='+1', species=species_params['species_name'])
                                    us_frameshifts = find_upstream_frameshift(feature, params['frame'], params['ustream_limit'], params['stop_codons'], species_params['signals'], args)
                                    if len(us_frameshifts) > 0:
                                        for fs in us_frameshifts:
                                            detected_frameshifts[feature.species].append(fs)
                                    ds_frameshifts = find_downstream_frameshift(feature, params['frame'], params['ustream_limit'], params['stop_codons'], species_params['signals'], args)
                                    if len(ds_frameshifts) > 0:
                                        for fs in ds_frameshifts:
                                            detected_frameshifts[feature.species].append(fs)
                                else:
                                    try:
                                        feat.qualifiers['protein_id'][0]
                                        if (feat.qualifiers['protein_id'][0]) == args.protein_id:
                                            print('found ',args.protein_id)
                                            feature = GenomeFeature(nucrec=nucrec, nuc_acc=nucrec_id, nuc_desc=nucrec_desc, feature=feat, strand='+1', species=species_params['species_name'])
                                            us_frameshifts = find_upstream_frameshift(feature, params['frame'], params['ustream_limit'], params['stop_codons'], species_params['signals'], args)
                                            if len(us_frameshifts) > 0:
                                                for fs in us_frameshifts:
                                                    detected_frameshifts[feature.species].append(fs)
                                            ds_frameshifts = find_downstream_frameshift(feature, params['frame'], params['ustream_limit'], params['stop_codons'], species_params['signals'], args)
                                            if len(ds_frameshifts) > 0:
                                                for fs in ds_frameshifts:
                                                    detected_frameshifts[feature.species].append(fs)
                                    except:
                                        nucrec.features[i]
                    # Same loop for the reverse strand
                    nucrec = nucrec.reverse_complement()
                    for i in range(len(nucrec.features)):
                        feat = nucrec.features[i]
                        if (feat.type == 'CDS'):
                            if (feat.strand == 1):
                                if args.protein_id == None:
                                    feature = GenomeFeature(nucrec=nucrec, nuc_acc=nucrec_id, nuc_desc=nucrec_desc, feature=feat, strand='-1', species=species_params['species_name'])
                                    us_frameshifts = find_upstream_frameshift(feature, params['frame'], params['ustream_limit'], params['stop_codons'], species_params['signals'], args)
                                    if len(us_frameshifts) > 0:
                                        for fs in us_frameshifts:
                                            detected_frameshifts[feature.species].append(fs)
                                    ds_frameshifts = find_downstream_frameshift(feature, params['frame'], params['ustream_limit'], params['stop_codons'], species_params['signals'], args)
                                    if len(ds_frameshifts) > 0:
                                        for fs in ds_frameshifts:
                                            detected_frameshifts[feature.species].append(fs)
                                else:
                                    try:
                                        feat.qualifiers['protein_id'][0]
                                        if (feat.qualifiers['protein_id'][0]) == args.protein_id:
                                            print('found ',args.protein_id)
                                            feature = GenomeFeature(nucrec=nucrec, nuc_acc=nucrec_id, nuc_desc=nucrec_desc, feature=feat, strand='-1', species=species_params['species_name'])
                                            us_frameshifts = find_upstream_frameshift(feature, params['frame'], params['ustream_limit'], params['stop_codons'], species_params['signals'], args)
                                            if len(us_frameshifts) > 0:
                                                for fs in us_frameshifts:
                                                    detected_frameshifts[feature.species].append(fs)
                                            ds_frameshifts = find_downstream_frameshift(feature, params['frame'], params['ustream_limit'], params['stop_codons'], species_params['signals'], args)
                                            if len(ds_frameshifts) > 0:
                                                for fs in ds_frameshifts:
                                                    detected_frameshifts[feature.species].append(fs)
                                    except:
                                        nucrec.features[i]
        # For each GenomeFeature in genome_features, search for upstream and downstream frameshifts
        for feature in genome_features:
            us_frameshifts = find_upstream_frameshift(feature, params['frame'], params['ustream_limit'], params['stop_codons'], species_params['signals'], args)
            if len(us_frameshifts) > 0:
                for fs in us_frameshifts:
                    detected_frameshifts[feature.species].append(fs)
            ds_frameshifts = find_downstream_frameshift(feature, params['frame'], params['ustream_limit'], params['stop_codons'], species_params['signals'], args)
            if len(ds_frameshifts) > 0:
                for fs in ds_frameshifts:
                    detected_frameshifts[feature.species].append(fs)
        
        if args.protein_id == None:
            for species_name in detected_frameshifts.keys():
                print('Writing results to', species_name)
                write_to_txt(params["results_dir"] + '/' + species_name , species_name )
                write_to_csv(params["results_dir"] + '/frameshifts', species_name)
    if args.protein_id == None:  
        print()      
        create_fasta_file_for_blast_db('frameshifts.fasta') 
        makeblastdb('frameshifts.fasta')
        blast_search()