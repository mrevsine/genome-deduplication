### New attempt to deduplicate a fasta

# sample_regions are inclusive at the start, exclusive at the end
# masked_starts are the start coordinates of the masked kmers
# skipped_regions are inclusive at the start, exclusive at the end

###=============================================================================
### Imports

import argparse
from Bio import SeqIO
import gzip
import json
import numpy as np
import os
import pickle
import random as rng
import re
import sys
import tracemalloc


###=============================================================================
### Helper functions

## Component Functions ================

def encode_kmer(kmer):
    char_map = {'A':0, 'C':1, 'G':2, 'T':3}
    kmer_num = 0
    for c in kmer:
        kmer_num = (kmer_num << 2) | char_map[c]
    return kmer_num


def decode_kmer(kmer_num, k=32):
    char_map = {0:'A', 1:'C', 2:'G', 3:'T'}
    kmer = []
    for i in range(k):
        nucleotide_code = kmer_num & 3 # 0x11
        kmer.append(char_map[nucleotide_code])
        kmer_num >>= 2
    return ''.join(reversed(kmer))


def get_fasta_basename(fasta):
    n_suffixes = 2 if fasta.endswith(".gz") else 1
    fasta_basename = '.'.join(os.path.basename(fasta).split('.')[:-(n_suffixes)])
    return fasta_basename


# Get total number of duplicate bases from seq_info array, considering value of k
# e.g.   sequence =  T G C A A G A A A C C A A A T G  and seen_kmers = {AAA, AAG, AGA}
#        seq_info = [1,1,1,4,4,1,4,1,1,1,1,4,1,1,1,1] and k = 3
#  duplicate_idxs = [0,0,0,1,1,1,1,1,1,0,0,1,1,1,0,0]
# duplicate count = 9
def get_duplicate_base_count(seq_info, k, duplicate_codes=[4,5]):
    duplicate_idxs = np.zeros(seq_info.shape[0], dtype=bool)
    for i,n in enumerate(seq_info):
        if n in duplicate_codes:
            duplicate_idxs[i:i+k] = True
    return duplicate_idxs.sum()


def get_cigar(seq_info):

    # 0 - nothing, 1 - unique (sample), 2 - omitted (ignored), 3 - ambiguous, 4 - local repeat, 5 - global repeat
    region_type_char_codes = ["X", "U", "O", "A", "L", "G"] 

    # Compute cigar string
    cigar = ""
    curr_type = None
    curr_start_idx = 0
    for i,n in enumerate(seq_info):
        if n != curr_type:
            if curr_type is not None:
                cigar += f"{i - curr_start_idx}{region_type_char_codes[curr_type]}"
            curr_type = n
            curr_start_idx = i
    cigar += f"{seq_info.shape[0] - curr_start_idx}{region_type_char_codes[curr_type]}"

    return cigar


# Group masked indices into contiguous bed regions
# e.g. [2,3,4,7,8,20] -> [(2,5), (7,9), (20,21)]
def condense_masked_regions(masked):
    masked_regions = []
    if len(masked) == 0:
        return masked_regions
    region_start = masked[0]
    for i in range(len(masked)-1):
        if masked[i] + 1 != masked[i+1]:
            masked_regions.append((region_start, masked[i] + 1))
            region_start = masked[i+1]
    masked_regions.append((region_start, masked[-1] + 1))
    return masked_regions


# Get all kmers from a given sample set
# Used when resuming deduplication from previous run to not have to re-duplicate completed inputs
def compute_seen_kmers_from_bed_and_fasta(samples_bed, fasta, k):
    seen_kmers = set()
    # Read in fasta sequences
    seq_dict = {}
    with open_maybe_gzip(fasta) as f:
        for record in SeqIO.parse(f, "fasta"):
            # Read sequence
            sequence = str(record.seq)
            # Clean sequence
            ambiguous_chars_regex = re.compile(r'[^ACGTN]')
            clean_sequence = ambiguous_chars_regex.sub('N', str(sequence.upper()))
            # Add to seq_dict
            seqname = record.id
            seq_dict[seqname] = clean_sequence
    # Read in bed file and extract kmers from sampled regions
    with open_maybe_gzip(samples_bed) as bedfile:
        for line in bedfile:
            fields = line.strip().split('\t')
            seqname = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            sequence = seq_dict.get(seqname, None)
            if sequence is None:
                raise ValueError(f"Error: sequence {seqname} found in bed file but not in fasta file.")
            for i in range(start, end - k + 1):
                kmer = sequence[i:i+k]
                if 'N' not in kmer:
                    kmer_num = encode_kmer(kmer)
                    seen_kmers.add(kmer_num)
    return seen_kmers


## I/O Functions ================

def type_check(file):
    """
    Check if file is gzipped or not
    """
    allowed_suffixes = [".gz", ".fasta", ".fa", ".fna", ".txt", ".list"]
    good_suffix = False
    for suffix in allowed_suffixes:
        if file.endswith(suffix):
            good_suffix = True
            break
    if not good_suffix:
        raise ValueError(f"Error: could not determine file type. Supported types are {', '.join(allowed_suffixes)}")


# Flexible function to open either a regular or gzipped file
def open_maybe_gzip(fname):
    with open(fname, "rb") as raw:
        signature = raw.read(2)
    if signature == b"\x1f\x8b":
        return gzip.open(fname, "rt")
    return open(fname, "rt")


# Extra term is a value that will be written to every row of the 4th column
def writeout_bed(name, regions, bedfile, extra_col_val=None):
    """
    Create bed files for regions
    """
    if len(regions) < 1:
        return
    extra_term = "" if extra_col_val is None else f"\t{extra_col_val}"
    for start,end in regions:
        bedfile.write(f"{name}\t{start}\t{end}{extra_term}\n")


def writeout_kmers(seen_kmer_set, file_basename):
    """
    Dump seen kmers into a pickle file that can be used as input for next dedup file
    """
    outfile = file_basename + ".kmers.pkl"
    with open(outfile, 'wb') as f:
        pickle.dump(seen_kmer_set, f)


def output_dump(local_dict, file_basename, k):
    """
    Output all bed files
    """
    with open(file_basename + ".samples.bed", 'w') as sample_regions:
        with open(file_basename + ".masks.bed", 'w') as masked_regions:
            with open(file_basename + ".ignored.bed", 'w') as skipped_regions:
                with open(file_basename + ".ambiguous.bed", 'w') as ambiguous_regions:
                    for seqname in local_dict.keys():
                        writeout_bed(seqname, local_dict[seqname]["sample_regions"], sample_regions)
                        writeout_bed(seqname, local_dict[seqname]["masked_regions"], masked_regions)
                        writeout_bed(seqname, local_dict[seqname]["skipped_regions"], skipped_regions)
                        writeout_bed(seqname, local_dict[seqname]["ambiguous_regions"], ambiguous_regions)


## Testing Functions ================

def test_with_fasta(fasta, k, sample_len, seen_kmers):
    local_dict = {}
    file_basename = fasta().split('.')[0]
    with open(fasta, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            sequence = str(record.seq)
            seqname = record.id
            sample_regions, masked_starts, skipped_regions, seen_kmers = deduplicate(sequence, k, sample_len, seen_kmers=seen_kmers)

            print(f"Sample regions: {sample_regions}")
            print(f"Masked starts: {masked_starts}")
            print(f"Skipped regions: {skipped_regions}")

            # Keep dict associating the regions with the sequence name
            local_dict[seqname] = {
                "sample_regions": sample_regions,
                "masked_starts": masked_starts,
                "skipped_regions": skipped_regions,
            }
            print(f"Seen kmers: {seen_kmers}")

        output_dump(local_dict, file_basename)

    writeout_kmers(seen_kmers, file_basename + f".pickle")


def test_with_sequence(args):
    with open(args.input[0], 'r') as f:
        seq = f.read().strip()
    sample_regions, masked_starts, skipped_regions, seen_kmers = deduplicate_seq(seq, args.seen_kmers, args)
    print("Sample regions: ", sample_regions)
    print("Masked starts: ", masked_starts)
    print("Skipped regions: ", skipped_regions)
    print("Seen kmers: ", seen_kmers)
    return sample_regions, masked_starts, skipped_regions, seen_kmers

## Core Deduplication Functions ================

# Assess a potential sequence as a whole
# All samples are analyzed in their entirety, but many are rejected
# We want to retain information from one sequence to the next for speed
# Therefore, we input and output information from the previous run
#
# Step 1: go back to the parts of seq_info before seq_info_offset that might have changed by nature of moving forward in the seq
#   - this is just local repeats, which may no longer be repeats if their original match has been advanced past
#   - all other types should be unchanged
#   - Step 1a:  create a first draft of sample_seen_kmers from all the positions in seq_info previously annotated as unique
#   - Step 1b: re-check all positions of seq_info previously annotated as local/internal repeats vs. sample_seen_kmers
#
# Step 2: from seq_info_offset to the end of the sample, label all new Ns as ambiguous
#
# Step 3: from seq_info_offset to the end of the sample, find all new contigs based on new positions of ambiguous bases
# 
# Step 4: iterate over all new contigs. In each contig, iterate over all kmers
# 
def check_sample_overall(seq, seq_info, seq_info_offset, seen_kmers, k, min_sample_len, no_overlap_samples, dedup_args):

    # if seq_info is None:
    #     seq_info = np.ones(len(seq))
    #     seq_info_offset = 0
    #     seq_info[-k+1:] = region_type_codes["unannotated"]

    ###=========================================================================
    ### Global data

    # Dictionary of region types to their int encodings, i.e. enums
    region_type_codes = {
        "unannotated": 0,
        "unique": 1,
        "ignored": 2,
        "ambiguous": 3, 
        "internal repeat": 4,
        "global repeat": 5,
    }

    # Relevant global variables
    final_start_idx = len(seq) - k
    sample_end_coord = len(seq)
    duplicate_start_idx = -1
    ambiguous_idx = -1
    ignored_region = None
    next_start_offset = len(seq) if no_overlap_samples else len(seq) - k + 1
    
    # First draft of sample_seen_kmers based on all positions in seq_info labelled as unique
    sample_seen_kmers = {}
    unique_idxs = [i for i,n in enumerate(seq_info[:min(seq_info_offset, final_start_idx+1)]) if n == region_type_codes["unique"]]
    for i in unique_idxs:
        kmer_num = encode_kmer(seq[i:i+k])
        if kmer_num not in sample_seen_kmers:
            sample_seen_kmers[kmer_num] = i

    ###=========================================================================
    ### Re-check all local repeats before seq_info_offset 
    ### This is the only info retained from the previous iteration that could be wrong now
    
    # First, account for all positions previously classified as internal repeats, which may no longer be accurate
    local_repeat_idxs = [i for i,n in enumerate(seq_info) if n == region_type_codes["internal repeat"]]
    for kmer_start_idx in local_repeat_idxs:

        # Get kmer
        kmer = seq[kmer_start_idx:kmer_start_idx+k]
        kmer_num = encode_kmer(kmer)

        # Add kmer to sample seen kmers if not seen before
        # We know that this position was not a global repeat in the prior round
        # The fact that we're seeing this position again means the last sample was not accepted, so the global kmers are the same
        # Actually, this might not be true ^ if we allow for some crazy overlap between samples
        # 
        # Also, we know that one of these kmers being included in sample_seen_kmers will never match at a later start position, 
        # or else that position would've been marked as a local repeat in the previous iteration. At this point, 
        # sample_seen_kmers is only comprised of kmers that were previously unique.
        # So any match here can only be w.r.t. a kmer that comes earlier in this sequence
        if kmer_num not in sample_seen_kmers:
            # sample_seen_kmers.add(kmer_num)
            sample_seen_kmers[kmer_num] = kmer_start_idx
            seq_info[kmer_start_idx] = region_type_codes["unique"] # Reset the classification at these positions
        
        # Will this case ever occur? Maybe a speedup by eliminating
        elif kmer_num in seen_kmers:
            seq_info[kmer_start_idx] = region_type_codes["global repeat"]

    ###=========================================================================
    ### Find new ambiguous bases after seq_info_offset

    # Look for new ambiguous bases
    # new_ambiguous_bases = [seq_info_offset + i for i,c in enumerate(seq[seq_info_offset:(final_start_idx+1)]) if c == 'N']
    new_ambiguous_bases = [seq_info_offset + i for i,c in enumerate(seq[seq_info_offset:]) if c == 'N']
    seq_info[new_ambiguous_bases] = region_type_codes["ambiguous"]

    ###=========================================================================
    ### Get new contigs after seq_info_offset

    # Get contigs based on ambiguous chars
    # Make sure no contig eclipses the seq len - k + 1 mark
    contig_boundaries = []
    contig_start = seq_info_offset
    for N_idx in new_ambiguous_bases:
        if N_idx > contig_start:
            contig_boundaries.append((contig_start, N_idx))
        contig_start = N_idx + 1
    if contig_start < len(seq):
        contig_boundaries.append((contig_start, len(seq)))

    ###=========================================================================
    ### Check all new contigs

    # Next, loop through all contigs
    for contig_start,contig_end in contig_boundaries:

        # The final k-1 bases are ignored
        # If the whole contig is shorter than k bases, the whole thing is ignored
        ignored_region_end = contig_end
        ignored_region_start = max(contig_start, contig_end - k + 1)
        seq_info[ignored_region_start:ignored_region_end] = region_type_codes["ignored"]

        # Find all kmers in this contig
        for kmer_start_idx in range(contig_start, contig_end - k + 1):

            # Get kmer
            kmer = seq[kmer_start_idx:kmer_start_idx+k]
            if 'N' in kmer:
                print(f"ERROR: N in kmer at {kmer_start_idx}, kmer = {kmer}")
                print(seq_info)
                print(contig_boundaries)
                print(contig_start, contig_end)
                print(seq)
                exit()
            kmer_num = encode_kmer(kmer)

            # Add kmer to sample seen kmers if not seen before 
            # And mark position as unique
            if kmer_num not in sample_seen_kmers and kmer_num not in seen_kmers:
                sample_seen_kmers[kmer_num] = kmer_start_idx
                seq_info[kmer_start_idx] = region_type_codes["unique"]

            # If kmer has been seen, add its position info to a list of duplicate regions
            else:
                if kmer_num in sample_seen_kmers:
                    seq_info[kmer_start_idx] = region_type_codes["internal repeat"]
                else:
                    seq_info[kmer_start_idx] = region_type_codes["global repeat"]

    ###=========================================================================
    ### Compute sequence-level metrics on ambiguous and duplicate content

    # Get total number of ambiguous bases
    # Account for the final k-1 bases that we skip up until now
    n_ambiguous_bases = sum(seq_info == region_type_codes["ambiguous"])

    # Compute total number of duplicates
    n_duplicate_kmers = sum((seq_info == region_type_codes["global repeat"]) | (seq_info == region_type_codes["internal repeat"]))
    n_duplicate_bases = get_duplicate_base_count(seq_info, k)

    ###=========================================================================
    ### Decide whether to keep this sample

    evaluation_method = dedup_args["evaluation_method"]
    # selection_strategy = "agnostic_retain_pct" # Option 1
    # selection_strategy = "content_threshold" # Option 2
    # selection_strategy = "content_motivated_retain_pct" # Option 3

    valid_sequence = True

    ## Option 1: keep or discard at some fixed % regardless of sequence content if sequence has duplication/ambiguity
    if evaluation_method == "per_sample_agnostic":

        # Input arg(s)
        agnostic_retain_pct = dedup_args["agnostic_retain_pct"]

        # Only flip coin if the sequence has duplication; all good sequences are accepted anyway
        if n_ambiguous_bases > 0 or (n_duplicate_bases > 0 and rng.random() >= agnostic_retain_pct):

            valid_sequence = False

    ## Option 2: keep all samples with deduplication and ambiguous char rates below thresholds
    elif evaluation_method == "per_sample_threshold":

        # Input arg(s)
        n_allowed_ambiguous_bases = dedup_args["ambiguous_base_threshold"] * len(seq)
        n_allowed_duplicate_bases = dedup_args["duplicate_base_threshold"] * len(seq)

        ## Could march backwards from seq len to min seq len, attempting to take the biggest possible sample

        if n_ambiguous_bases > n_allowed_ambiguous_bases or n_duplicate_bases > n_allowed_duplicate_bases:

            valid_sequence = False

    ## Option 3: chance of discarding is a function of the amount of deduplication or ambiguous chars
    ## Not currently a viable selection in the code, but may be implemented in the future
    elif evaluation_method == "content_motivated_retain_pct":

        # Input arg(s)
        def content_motivated_retain_pct(n_duplicate, n_ambiguous):
            return (n_duplicate + n_ambiguous) / len(seq)
        
        if (n_ambiguous_bases > 0 or n_duplicate_bases > 0) \
            and rng.random() >= content_motivated_retain_pct(n_duplicate_bases, n_ambiguous_bases):
            
            valid_sequence = False
        
    ###=========================================================================
    ### Find the first problematic region for determining the cause of the sample being rejected

    if not valid_sequence:

        # # Always advance by 1 base, since it could chance the overall metrics of the sequence
        # # Record initial base
        # if seq_info[0] == region_type_codes["ambiguous"]:
        #     ambiguous_idx = 0
        # elif seq_info[0] == region_type_codes["internal repeat"]:
        #     ignored_region = (0,1)
        # elif seq_info[0] == region_type_codes["global repeat"]:
        #     duplicate_start_idx = 0
        # next_start_offset = 1

        first_issue_base = None
        for i,n in enumerate(seq_info):
            if n in [region_type_codes["ambiguous"], region_type_codes["internal repeat"], region_type_codes["global repeat"]]:
                first_issue_base = (i,n)
                break

        # If an ambiguous character comes first:
        if first_issue_base[1] == region_type_codes["ambiguous"]:
            ambiguous_idx = first_issue_base[0]
            if ambiguous_idx > 0:
                ignored_region = (0, ambiguous_idx)
            next_start_offset = ambiguous_idx + 1
        
        # If an internal repeat comes first
        elif first_issue_base[1] == region_type_codes["internal repeat"]:

            # Use the original kmer as the skipping point, not its local repeat
            repeat_kmer = seq[first_issue_base[0]:(first_issue_base[0]+k)]
            repeat_kmer_num = encode_kmer(repeat_kmer)
            repeat_first_idx = sample_seen_kmers[repeat_kmer_num] # The position of the first instance of this repeated kmer

            # Record the entire region up to and including the original kmer as ignored
            # Importantly, we do not denote the original kmer as a duplicate, since it is 
            # both not in the global set and is not added to a sample      
            ignored_region = (0, repeat_first_idx + 1)
            # The next sample should start at the original kmer's position + 1
            next_start_offset = repeat_first_idx + 1

        # If a global repeat comes first
        else:
            duplicate_start_idx = first_issue_base[0]
            if duplicate_start_idx > 0:
                ignored_region = (0, duplicate_start_idx)
            next_start_offset = duplicate_start_idx + 1

        sample_end_coord = -1
        sample_seen_kmers = None

    ###=========================================================================
    ### Return needed info

    return sample_end_coord, duplicate_start_idx, ambiguous_idx, ignored_region, next_start_offset, sample_seen_kmers, seq_info


# Function to check a sample for duplicates, allowing some rate of duplicates
def check_sample_per_kmer(seq, seen_kmers, k, allowed_duplicate_rate, min_sample_len, no_overlap_samples, N_policy="allow_none"):

    ## Initialize all data structures returned by this function

    # Data structure for new kmers in this sample
    sample_seen_kmers = {}
    
    # Sequence attributes
    sample_end_coord = len(seq)
    duplicate_start_idx = -1
    ambiguous_idx = -1
    ignored_region = None
    next_start_offset = len(seq) if no_overlap_samples else len(seq) - k + 1
    
    # # Check for Ns
    if N_policy == "allow_none":
        N_index = seq.find('N')
        # If we find an N
        if N_index > -1:
            # If the first N occurs before the minimum sample length, there's no valid sample here; handle like a failed sample
            if N_index < min_sample_len:
                # Record position of the N in the list of ambiguous characters
                ambiguous_idx = N_index
                # Record entire region before the N as skipped
                if N_index > 0:
                    ignored_region = (0, N_index)
                next_start_offset = N_index + 1
                sample_end_coord = -1
            else:
                # Set N index to be the end of this possible sample
                sample_end_coord = N_index

    # Loop through all kmers in this possible sample
    if sample_end_coord >= min_sample_len:
        for kmer_start_idx in range(sample_end_coord-k+1):

            # Get kmer and its encoding for storage
            kmer = seq[kmer_start_idx:kmer_start_idx+k]
            kmer_num = encode_kmer(kmer)

            # If we haven't seen this kmer before, record it in the sample kmers
            if kmer_num not in seen_kmers and kmer_num not in sample_seen_kmers:
                sample_seen_kmers[kmer_num] = kmer_start_idx

            # If we've seen this kmer before, stop analyzing this sample
            else:

                # Opportunity to process this duplicate with some function to determine whether to ignore it
                # Simple chance, function based on previous occurrence count, etc.
                if allowed_duplicate_rate > 0 and rng.random() < allowed_duplicate_rate:
                    continue

                # At this point, we have found a definite duplicate

                # If internal repeat (offender is a repeat of original):
                    # If offender comes before min_sample_len:
                        # Although we don't save sample seen kmers, we know that starting again 
                        # before the original kmer will run into the same repeat
                        # Therefore, next sample needs to start after the original
                        # -> Invalid sample
                        # -> Ignored region from start to original + 1
                        # -> No duplicate found
                        # -> Next start = original + 1
                    # else:
                        # We will save this sample and add its kmers to the global set, so any new sample 
                        # before the offender will be invalidated by the original
                        # Therefore, next sample needs to start after the offender
                        # -> valid sample up to offender position
                        # -> Next start = offender + 1
                # else (repeat matches the global, not the local, set)
                    # If offender comes before min_sample_len:
                        # -> Invalid sample
                        # -> Next start = offender + 1
                    # else:
                        # -> valid sample up to offender position
                        # -> Next start = offender + 1

                # If the offending kmer is within min_sample_len
                if kmer_start_idx < min_sample_len:

                    # This is an invalid sample
                    sample_end_coord = -1

                    # Check if this is an internal duplicate
                    repeat_match_idx = sample_seen_kmers.get(kmer_num, -1)

                    # If an internal match
                    if repeat_match_idx >= 0:

                        # The next sample should start at the original kmer's position + 1
                        next_start_offset = repeat_match_idx + 1

                        # Record the entire region up to and including the original kmer as ignored
                        # Importantly, we do not denote the original kmer as a duplicate, since it is
                        # both not in the global set and is not added to a sample
                        ignored_region = (0, next_start_offset)

                    # Else if a match to a global kmer
                    else:

                        # The next sample should start at the current position + 1
                        next_start_offset = kmer_start_idx + 1

                        # Record the current kmer as a duplicate
                        duplicate_start_idx = kmer_start_idx

                        # Record the entire region up to but not including the offending kmer as ignored
                        if kmer_start_idx > 0:
                            ignored_region = (0, kmer_start_idx)

                    # Since this sample is invalid, we can discard its seen kmers
                    sample_seen_kmers = None

                # If the offending kmer is past min_sample_len
                # Whether the repeat is internal or global, the current kmer is the offending one
                else:

                    # The current kmer is the end boundary of the sample
                    sample_end_coord = kmer_start_idx

                    # The current kmer is a duplicate
                    duplicate_start_idx = kmer_start_idx

                    # The next sample should start at 1 + the current position
                    next_start_offset = kmer_start_idx + 1
                
                # We have finished processing this sample
                break

    return sample_end_coord, duplicate_start_idx, ambiguous_idx, ignored_region, next_start_offset, sample_seen_kmers


###=============================================================================
### Main functions

## Main Deduplication ================
def deduplicate_seq(seq, seen_kmers, args):

    ###=========================================================================
    ### Set up needed data

    # Collect needed args
    k = args.kmer
    sample_len = args.sample_len
    min_sample_len = args.min_sample_len
    no_overlap_samples = args.no_overlap
    evaluation_method = args.evaluation_method
    seed = args.seed

    # Instantiate seen_kmers if None
    if seen_kmers is None:
        seen_kmers = set()

    # Set up storage for samples, masked regions, and ignored regions
    sample_regions = []
    masked_starts = []
    ignored_regions = []
    ambiguous_positions = []

    # Evaluation method-specific needed data
    if evaluation_method == "per_kmer":
        allowed_duplicate_rate = args.per_kmer_retain_pct
    else:
        sample_seq_info = np.array([])
        seq_info_offset = 0
        dedup_args = {
            "evaluation_method": evaluation_method,
            "ambiguous_base_threshold": args.ambiguous_base_threshold,
            "duplicate_base_threshold": args.duplicate_base_threshold,
            "agnostic_retain_pct": args.agnostic_retain_pct,
        }

    ###=========================================================================
    ### Early exit conditions
    ### If input sequence is too short, return null values

    # Check validity of inputs
    if min_sample_len is not None and len(seq) < min_sample_len:
        print("Warning: the sequence length is less than min_sample_len. Skipping this sequence.")
        return [], [], [(0,len(seq))], [], seen_kmers 
    if min_sample_len is None and len(seq) < sample_len:
        print("Warning: the sequence length is less than sample_len. Skipping this sequence.")
        return [], [], [(0,len(seq))], [], seen_kmers 

    ###=========================================================================
    ### Main loop logic; process each potential sample in the sequence

    # Needed global vars
    max_start_idx = len(seq) - min_sample_len # final possible start index
    sample_start = 0 # start index of the current sample

    # Investigate every possible sample
    while sample_start <= max_start_idx:

        # Get boundary for this possible sample
        sample_end = min(len(seq), sample_start + sample_len) # Checking against seq len is only necessary for final sample

        ###=====================================================================
        ### Evaluate this sample

        # In any sample mode, instantiate seq info
        if evaluation_method in ["per_sample_agnostic", "per_sample_threshold"]:
            seqlen = sample_end - sample_start
            seq_info = np.zeros(shape=seqlen, dtype=int) # 0 means unannotated
            seq_info[:sample_seq_info.shape[0]] = sample_seq_info

        ###=====================================================================
        ### Evaluate this sample
        
        if evaluation_method == "per_kmer":
            checked_sample_len, duplicate_start_idx, ambiguous_idx, ignored_region, next_start_offset, sample_seen_kmers = \
                check_sample_per_kmer(seq[sample_start:sample_end], seen_kmers, k, allowed_duplicate_rate, min_sample_len, no_overlap_samples)
        else:
            checked_sample_len, duplicate_start_idx, ambiguous_idx, ignored_region, next_start_offset, sample_seen_kmers, sample_seq_info = \
                check_sample_overall(seq[sample_start:sample_end], seq_info, seq_info_offset, seen_kmers, k, min_sample_len, no_overlap_samples, dedup_args)

        ###=====================================================================
        ### Update regions with result of evaluation method call

        # Add sample to samples if valid
        if checked_sample_len > -1:
            sample_regions.append((sample_start, sample_start + checked_sample_len))
        
        # Record the duplicate if any was found
        if duplicate_start_idx > -1:
            masked_starts.append(sample_start + duplicate_start_idx)

        # Record the ambiguous index if any was found
        if ambiguous_idx > -1:
            ambiguous_positions.append(sample_start + ambiguous_idx)

        # Add ignored region if any was found
        if ignored_region is not None:
            ignored_regions.append((sample_start+ignored_region[0], sample_start+ignored_region[1]))

        # Update seen kmers with any sample-specific kmers
        if sample_seen_kmers is not None:
            seen_kmers.update(sample_seen_kmers)

        # Set start index of next iteration of the loop
        sample_start = sample_start + next_start_offset

        ###=====================================================================
        ### Update regions with result of evaluation method call

        # Account for seq_info array in overall mode
        if evaluation_method in ["per_sample_agnostic", "per_sample_threshold"]:

            # next_start_offset is the first position where we are clear of the last issue, where we will start the next sample
            # next_iteration_amount is the number of bases to check in the next iteration
            # The final next_iteration_amount bases sample_seq_info will start as unannotated
            # Must always be at least k, so that we look at at least 1 kmer next time. Otherwise there's no point
            # Must always be less than or equal to sample length, or else we will look at too much sequence
            next_iteration_amount = min(sample_len, next_start_offset + k - 1)
            seq_info_offset = sample_seq_info.shape[0] - next_iteration_amount
            sample_seq_info = sample_seq_info[next_start_offset:]

    ###=========================================================================
    ### Final housekeeping 

    # Record possible ignored sequence at the end
    if sample_start < len(seq):
        ignored_regions.append((sample_start, len(seq)))

    # Convert masked starting indices and ambiguous positions to regions for more condensed bed files
    masked_regions = condense_masked_regions(masked_starts)
    ambiguous_regions = condense_masked_regions(ambiguous_positions)

    return sample_regions, masked_regions, ignored_regions, ambiguous_regions, seen_kmers


def deduplicate_genome(fasta, seen_kmers, save_kmers_to_file, args):

    # Keep dict associating the regions with the sequence name
    local_dict = {} 

    # Set output prefix
    fasta_basename = get_fasta_basename(fasta)
    out_prefix = os.path.join(args.output_dir, fasta_basename)

    # Read in fasta
    with open_maybe_gzip(fasta) as f:
        for record in SeqIO.parse(f, "fasta"):
            sequence = str(record.seq)
            seqname = record.id

            # Convert all chars to uppercase and replace ambiguous chars with N
            ambiguous_chars_regex = re.compile(r'[^ACGTN]')
            clean_sequence = ambiguous_chars_regex.sub('N', str(sequence.upper()))

            # Deduplicate this sequence
            sample_regions, masked_regions, skipped_regions, ambiguous_regions, seen_kmers = deduplicate_seq(clean_sequence, seen_kmers, args)

            print(f"Seen kmer size: {sys.getsizeof(seen_kmers)}")

            # Associate deduplication data with the sequence name
            local_dict[seqname] = {
                "sample_regions": sample_regions,
                "masked_regions": masked_regions,
                "skipped_regions": skipped_regions,
                "ambiguous_regions": ambiguous_regions
            }

        output_dump(local_dict, out_prefix, args.kmer)

    if save_kmers_to_file:
        writeout_kmers(seen_kmers, out_prefix)

    #print(f"Done with {fasta}")
    return seen_kmers


def deduplicate(args):
    
    # Create output directory if it doesn't already exist
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)
    else:
        delete_check=input("Output directory already exists. Continue and potentially overwrite files? (y/n): ")
        if delete_check.lower() != 'y' and delete_check.lower() != 'yes':
           print("Exiting...")
           sys.exit(1)

    # Write (processed) input args to file for reproducibility
    with open(os.path.join(args.output_dir, "config.json"), 'w') as f:
        json.dump(vars(args), f, indent=4)

    # Read input file of genome locations
    # Single fasta input or series of individual fasta entries
    if len(args.input) > 1 or args.input[0].endswith(('.fa', '.fasta', '.fasta.gz', '.fna', '.fna.gz')):
        fastas = args.input
    # File with a list of fasta files as input
    elif len(args.input) == 1 and args.input[0].endswith(('.txt', '.list')):
        with open(args.input[0], 'r') as f:
            fastas = [line.rstrip('\n') for line in f.readlines()]

    # Filter down to only the fastas that we can find and issue warning about
    # any we can't find
    valid_fastas = [fasta for fasta in fastas if os.path.isfile(fasta)]
    invalid_fastas = list(set(fastas).difference(set(valid_fastas)))
    if len(invalid_fastas) > 0:
        print("Warning: could not find the following fastas:")
        print(invalid_fastas)

    # Write basename to file map for all valid fastas
    basename_fasta_file = os.path.join(args.output_dir, "basename_fasta_match.txt")
    with open(basename_fasta_file, 'w') as f:
        for fasta in valid_fastas:
            fasta_basename = get_fasta_basename(fasta)
            f.write(f"{fasta_basename}\t{fasta}\n")

    # Initialize seen kmers here
    if args.seen_kmers is None:
        seen_kmers = set()
    else:
        print("Loading seen kmers")
        seen_kmers = pickle.load(open(args.seen_kmers, "rb"))

    # Iterate over fastas and deduplicate each one
    for i,fasta in enumerate(valid_fastas):

        # Deduplicate this fasta
        # print(f"Deduplicating {fasta}")
        save_kmers_to_file = args.save_every > 0 and (i+1) % args.save_every == 0

        # If we already completed this file, compute its kmer set and advance
        fasta_basename = get_fasta_basename(fasta)
        samples_bed = os.path.join(args.output_dir, f"{fasta_basename}.samples.bed")
        if os.path.exists(samples_bed):
            print(f"Found existing output for {fasta}, skipping deduplication and computing seen kmers from existing samples")
            fasta_seen_kmers = compute_seen_kmers_from_bed_and_fasta(samples_bed, fasta, args.kmer)
            seen_kmers.update(fasta_seen_kmers)
            continue

        # Otherwise, deduplicate as normal
        print(f"Deduplicating {fasta}")
        seen_kmers = deduplicate_genome(fasta, seen_kmers, save_kmers_to_file, args)

    # Optionally save seen kmers at the end
    if not args.no_save_kmers_at_end:
        final_seen_kmer_file_basename = "final"
        idx_suffix = 1
        while os.path.exists(os.path.join(args.output_dir, f"{final_seen_kmer_file_basename}.kmers.pkl")):
            final_seen_kmer_file_basename = f"final_{idx_suffix}"
            idx_suffix += 1
        writeout_kmers(seen_kmers, os.path.join(args.output_dir, final_seen_kmer_file_basename))


###=============================================================================
### Call to main

def __main__():

    tracemalloc.start()

    ## Collect input args
    parser = argparse.ArgumentParser()
    parser.add_argument("input", nargs="+", help="Input list of FASTA files or a txt file with one FASTA file per line")
    parser.add_argument("-a", "--agnostic_retain_pct", type=float, default=0.0, 
                        help="Likelihood a sample with any duplication will be retained in per_sample_agnostic evaluation mode")
    parser.add_argument("-b", "--ambiguous_base_threshold", type=float, default=0.0, 
                        help="pct. of bases allowed to be ambiguous in a sample for it to be retained in per_sample_threshold evaluation mode")
    parser.add_argument("-d", "--duplicate_base_threshold", type=float, default=0.0, 
                        help="pct. of bases allowed to be duplicate in a sample for it to be retained in per_sample_threshold evaluation mode")
    parser.add_argument("-e", "--evaluation_method", type=str, default="per_kmer", choices=["per_kmer", "per_sample_agnostic", "per_sample_threshold"], 
                        help="Method for deduplication evaluation (default: per_kmer)")
    parser.add_argument("-k", "--kmer", type=int, default=32, help="Kmer size (default: 32)")
    parser.add_argument("-l", "--sample_len", type=int, default=1000, help="Sample length (default: 1000)")
    parser.add_argument("-m", "--min_sample_len", type=int, default=None, help="Minimum sample length (default: 50)")
    parser.add_argument("-o", "--output_dir", default="dedup_out", help="Output directory (default: dedup_out/)")
    parser.add_argument("-p", "--seen_kmers", default=None, help="Pickle file containing seen kmers (default: None)")
    parser.add_argument("-r", "--per_kmer_retain_pct", type=float, default=0.0, 
                        help="Likelihood a duplicate kmer will be allowed through in per_kmer evaluation mode")
    parser.add_argument("-s", "--save_every", type=int, default=0, help="Save seen kmer set every n samples (default: 0, don't save any before the end)")
    parser.add_argument("--no-overlap", action="store_true", help="Keep neighboring samples discrete rather than overlapping by k-1 bases")
    parser.add_argument("--no-save-kmers-at-end", action="store_true", help="Opt out of saving the seen kmer set at the end of program execution")
    parser.add_argument("-seed", "--seed", type=int, default=123, help="Random seed for reproducibility")
    # Hidden test function -- pass in a file (.fa or .txt) here with expected results to verify correctness
    parser.add_argument("-T", "--test", action="store_true", help=argparse.SUPPRESS)
    # Hidden test kmer set input for testing. If not included, will start with empty kmer set
    parser.add_argument("-I", "--test-input-kmers", type=str, help=argparse.SUPPRESS)
    parser.add_argument("-O", "--test-output-kmers", type=str, help=argparse.SUPPRESS)
    args = parser.parse_args()

    # Set random seed for reproducibility
    rng.seed(args.seed)

    ## Process input args, checking for validity

    # Input fasta must be a valid file
    for f in args.input:
        if not os.path.isfile(f):
            raise("Error: could not find supplied list of fasta files")
        type_check(f)

    # Kmer must be a positive number
    if args.kmer < 1:
        raise("Error: kmer size must be a positive integer")
    
    # If kmer is very small or large, issue a warning
    if args.kmer < 16:
        print("Warning: small kmer sizes will result in very strict deduplication. Consider increasing the kmer size in the range 16-32 (default: 32)")
    if args.kmer > 32:
        print("Warning: kmer sizes above 32 may not work or may result in very slow evaluation. Consider decreasing the kmer size in the range 16-32 (default: 32)")

    # Check parameters related to evaluation method
    if args.evaluation_method not in ["per_kmer", "per_sample_agnostic", "per_sample_threshold"]:
        raise("Error: evaluation method must be one of per_kmer, per_sample_agnostic, or per_sample_threshold")
    if args.evaluation_method == "per_kmer" and (args.per_kmer_retain_pct > 1.0 or args.per_kmer_retain_pct < 0.0):
        raise("Error: per kmer retain pct rate must be between 0.0 and 1.0")
    if args.evaluation_method == "per_sample_agnostic" and (args.agnostic_retain_pct < 0.0 or args.agnostic_retain_pct > 1.0):
        raise("Error: agnostic keep pct must be between 0.0 and 1.0")
    if args.evaluation_method == "per_sample_threshold" and (args.ambiguous_base_threshold < 0.0 or args.ambiguous_base_threshold > 1.0):
        raise("Error: ambiguous base threshold must be between 0.0 and 1.0")
    if args.evaluation_method == "per_sample_threshold" and (args.duplicate_base_threshold < 0.0 or args.duplicate_base_threshold > 1.0):
        raise("Error: duplicate base threshold must be between 0.0 and 1.0")

    # Set min kmer length to sample length if None
    if args.min_sample_len is None:
        args.min_sample_len = args.sample_len

    # Seen kmers must be either none or a valid pickle file
    if args.seen_kmers is not None:
        if not os.path.isfile(args.seen_kmers):
            raise("Error: could not find supplied seen kmers file")

    # If sample len is smaller than k, issue warning that no deduplication will occur. Not an error bc this can be a way to create the raw sample sets
    if args.sample_len < args.kmer:
        print("Warning: sample len is smaller than k, meaning no deduplication will occur")

    # If min sample len is higher than sample len, set min to equal the default sample len
    if args.min_sample_len > args.sample_len:
        print("Warning: min sample length cannot be bigger than the default sample length; defaulting to equal the standard sample length")
        args.min_sample_len = args.sample_len

    ## Run deduplication
    print(args)
    if not args.test:
        print("Running deduplication")
        deduplicate(args)
    else:
        input=args.input[0]
        if input.endswith(".txt"):
            with open(input) as f:
                if f.readline().endswith(".fa\n") or f.readline().endswith(".fasta\n") or f.readline().endswith(".fna\n") or f.readline().endswith(".fa.gz\n") or f.readline().endswith(".fasta.gz\n") or f.readline().endswith(".fna.gz\n"):
                    print("Testing with fasta")
                    test_with_fasta(input, args.kmer, args.sample_len, input_kmers=set())
            print("Testing with sequence")
            test_with_sequence(args)
        else:
            print("Testing with fasta")
            file_basename = input().split('.')[0]
            file_outputname=file_basename + ".pickle"
            input_kmers = pickle.load(open(args.test_input_kmers, "rb")) if args.test_input_kmers is not None else set()
            assert args.ouptput_kmers is not None, "Error: must provide output kmer file when testing"
            output_compare = pickle.load(open(args.output_kmers, "rb"))
            test_with_fasta(input, args.kmer, args.sample_len, input_kmers)
            assert output_compare == pickle.load(file_outputname), "Error: output kmers do not match expected output"


if __name__ == "__main__":
    __main__()

