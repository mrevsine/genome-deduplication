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


###=============================================================================
### Helper functions

## Classes ================

class SeqInfo:

	def __init__(self, seqlen, k, encoding_dict=None):

		if encoding_dict is None:
			encoding_dict = {
				"unannotated": 0,
				"unique": 1, #
				"ignored": 2,
				"ambiguous": 3, #
				"internal repeat": 4,
				"global repeat": 5 #
			}
		self.encoding_dict = encoding_dict
		self.arr = np.zeros(seqlen, dtype=np.uint8) 
		self.n = seqlen
		self.current_sample_length = seqlen
		self.k = k
		self.current_global_idx = 0
		self.sample_seen_kmers = {}

		# The first position in seq_info's array that is new to this sample
		# i.e. analyzing sample chr1:1000-2000, but we already have data on 
		# chr1:1000-1500 from the previous sample, then seq_info_offset will be 500, 
		# since the first 500 bases of this sample are already annotated in seq_info
		self.offset = 0

		# Internally, self.arr is treated as a circular array
		# self.array_start is the position in the array that we treat as index 0
		self.array_start = 0

	def record_kmer(self, kmer_num, relative_idx):
		self.sample_seen_kmers[kmer_num] = self.current_global_idx + relative_idx

	def get_kmer_idx(self, kmer_num):
		global_kmer_idx = self.sample_seen_kmers.get(kmer_num, -1)
		if global_kmer_idx == -1:
			return -1
		else:
			return global_kmer_idx - self.current_global_idx

	def reset_annotations(self, n_bases):
		self[:n_bases] = self.encoding_dict["unannotated"]

	def purge_kmers(self, min_idx=None):
		if min_idx is None:
			self.sample_seen_kmers = {}    
		else:
			if not isinstance(min_idx, int):
				raise ValueError(f"Error: min_idx must be an integer, but got {type(min_idx)}")
			for kmer, idx in list(self.sample_seen_kmers.items()):
				if idx - self.current_global_idx < min_idx:
					del self.sample_seen_kmers[kmer]

	def move_internal_pointer(self, n_bases, last_sample_accepted=False):
		self.reset_annotations(n_bases)
		if last_sample_accepted:
			self.purge_kmers()
			# self.offset = 0
		else:
			self.purge_kmers(n_bases)
			# self.offset = max(0, self.current_sample_length - (self.k-1) - n_bases)
		self.offset = 0
		self.array_start = self._resolve_index(n_bases)

	def _resolve_index(self, idx):
		if idx < 0 or idx >= self.n:
			raise IndexError(f"index {idx} is out of bounds for size {self.n}")
		adj_idx = self.array_start + idx
		if adj_idx >= self.n:
			adj_idx -= self.n
		return adj_idx

	def _resolve_slice(self, s):
		indices = range(*s.indices(self.n))
		return np.array([self._resolve_index(i) for i in indices])

	def __getitem__(self, key):
		if isinstance(key, slice):
			return self.arr[self._resolve_slice(key)]
		return self.arr[self._resolve_index(key)]

	def __setitem__(self, key, value):
		if isinstance(key, slice):
			self.arr[self._resolve_slice(key)] = value
		else:
			self.arr[self._resolve_index(key)] = value


	def get_cigar(self, n_bases=None):

		region_type_char_codes = ["X", "U", "O", "A", "L", "G"] 
		cigar_arr = self[:self.current_sample_length] if n_bases is None else self[:(min(n_bases, self.current_sample_length))]
		cigar = ""
		curr_type = None
		curr_start_idx = 0
		for i,n in enumerate(cigar_arr):
			if n != curr_type:
				if curr_type is not None:
					cigar += f"{i - curr_start_idx}{region_type_char_codes[curr_type]}"
				curr_type = n
				curr_start_idx = i
		cigar += f"{self.n - curr_start_idx}{region_type_char_codes[curr_type]}"

		return cigar

	
	def __len__(self):
		return self.n


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
	for _ in range(k):
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

# Get contigs of unambiguous sequence from a chromosome, which may be used for downstream deduplication
# Return these contigs as a list of (start, end, [internal N region tpls]) coordinates
# Also return N regions for write to ambiguous.bed

def get_contigs_from_chromosome(seq, N_consecutive_allowed_Ns):
		
	# Get indices of all Ns in this sequence
	N_idxs = [i for i,c in enumerate(seq) if c == 'N']

	# Condense consecutive Ns into regions
	N_regions = condense_masked_regions(N_idxs)

	# Get valid contigs based on N regions and the allowed number of consecutive Ns
	contigs = []
	start_idx = 0
	internal_N_regions = []
	for N_start, N_end in N_regions:
		N_region_length = N_end - N_start
		if N_region_length > N_consecutive_allowed_Ns:
			if start_idx < N_start:
				contigs.append((start_idx, N_start, internal_N_regions))
			start_idx = N_end 
			internal_N_regions = []
		else:
			internal_N_regions.append((N_start, N_end))
	if start_idx < len(seq):
		contigs.append((start_idx, len(seq), internal_N_regions))

	# Return valid contigs and N regions
	return contigs, N_regions


# Function to check a sample for duplicates, allowing some rate of duplicates
def check_sample(seq, internal_N_regions, seen_kmers, k, dedup_parameter, min_sample_len, no_overlap_samples, evaluation_mode):

	###=========================================================================
	### Global data

	# Data structure for new kmers in this sample
	sample_seen_kmers = {}
	
	# Sequence attributes that will be returned
	sample_end_coord = len(seq)
	duplicate_start_idx = -1
	ignored_regions = []
	skipped_region = None
	next_start_offset = len(seq) if no_overlap_samples else len(seq) - k + 1

	###=========================================================================
	### Using the internal N regions, get all ranges of start indices of valid kmers 
	### Also record ignored regions in the k-1 bases before ambiguous chars

	# Get contigs based on ambiguous chars
	valid_kmer_start_ranges = []
	valid_region_start = 0
	ignored_regions_before_ambiguous_bases = []
	for N_start, N_end in internal_N_regions:
		final_kmer_start = N_start - k

		# Record final k-1 bases before N as ignored
		ignored_region_start = max(valid_region_start, final_kmer_start + 1)
		if ignored_region_start < N_start:
			ignored_regions_before_ambiguous_bases.append((ignored_region_start, N_start))

		# Record region before this as valid
		if final_kmer_start >= valid_region_start:
			valid_kmer_start_ranges.append((valid_region_start, final_kmer_start + 1))

		valid_region_start = N_end

	# No need to record anything for final k-1 bases, which should remain unannotated
	# Just record final valid region
	final_kmer_start = len(seq) - k
	if final_kmer_start >= valid_region_start:
		valid_kmer_start_ranges.append((valid_region_start, final_kmer_start + 1))

	###=========================================================================
	### Accumulate all needed variables regardless of evaluation mode

	# Used by all methods
	done_evaluating = False
	found_first_duplicate = False
	first_duplicate_idx = -1 # Index of first instance of duplicate kmer
	first_duplicate_match_idx = -1 # Index of first instance of local repeat

	# Only for per-sample agnostic mode
	already_decided_agnostic = False

	# Only for per-sample threshold mode
	n_duplicate_kmers = 0 # Currently unused, but could be used in a future version
	n_duplicate_bases = 0
	last_duplicate_end_idx = -1
	n_allowed_full_length_duplicate_bases = int(np.floor(dedup_parameter * len(seq)))
	below_allowed_duplication_threshold = True
	exceed_duplication_threshold_idx = -1

	###=========================================================================
	### Check all contigs

	for valid_start, valid_end in valid_kmer_start_ranges:

		if done_evaluating:
			break

		for kmer_start_idx in range(valid_start, valid_end):

			# Get kmer and its encoding for storage
			kmer = seq[kmer_start_idx:kmer_start_idx+k]
			kmer_num = encode_kmer(kmer)

			# If we haven't seen this kmer before, record it in the sample kmers
			duplicate_match_idx = sample_seen_kmers.get(kmer_num, -1)
			if duplicate_match_idx == -1 and kmer_num not in seen_kmers:
				sample_seen_kmers[kmer_num] = kmer_start_idx

			# If we've seen this kmer before, potentially stop analyzing this sample
			elif not already_decided_agnostic:

				# In per-kmer mode, decide whether to reject this same based on a coin flip with respect to this kmer
				if evaluation_mode == "per_kmer":

					# Continue analyzing kmers, "pretend" this wasn't a duplicate
					if dedup_parameter > 0 and rng.random() < dedup_parameter:
						continue 

					# Done with sample; decide whether to accept or reject based on whether we have reached min_sample_len
					else:
						done_evaluating = True 
						if not found_first_duplicate:
							found_first_duplicate = True
							first_duplicate_idx = kmer_start_idx
							first_duplicate_match_idx = duplicate_match_idx

				# In per-sample agnostic mode, decide on the whole sample based on a coin flip when we reach the first offending kmer
				elif evaluation_mode == "per_sample_agnostic":

					if dedup_parameter > 0 and rng.random() < dedup_parameter:
						already_decided_agnostic = True # Accept sample immediately, but keep analyzing kmers just to record them in sample_seen_kmers
						continue
					else:
						done_evaluating = True # Don't analyze any more kmers, accept or reject now based on whether we have reached min_samaple_len
						if not found_first_duplicate:
							found_first_duplicate = True
							first_duplicate_idx = kmer_start_idx
							first_duplicate_match_idx = duplicate_match_idx

				elif evaluation_mode == "per_sample_threshold":

					# Record possible first duplicate
					if not found_first_duplicate:
						found_first_duplicate = True
						first_duplicate_idx = kmer_start_idx
						first_duplicate_match_idx = duplicate_match_idx

					# Update duplication counts
					n_duplicate_kmers += 1
					if kmer_start_idx >= last_duplicate_end_idx:
						n_duplicate_bases += k
					else:
						n_duplicate_bases += k - (last_duplicate_end_idx - kmer_start_idx)
					last_duplicate_end_idx = kmer_start_idx + k

					# Make a note if we have yet to exceed the allowed rate of duplication
					# This can be used to salvage a sample < len but >= min_sample_len
					if (n_duplicate_bases / (kmer_start_idx+k)) <= dedup_parameter:
						below_allowed_duplication_threshold = True
						exceed_duplication_threshold_idx = -1
					elif below_allowed_duplication_threshold == True:
						below_allowed_duplication_threshold = False
						exceed_duplication_threshold_idx = kmer_start_idx

					# Early exit if we have already exceeded the max allowed duplication
					if n_duplicate_bases > n_allowed_full_length_duplicate_bases:
						done_evaluating = True

				# If the previous step causes us to reject the sample already, record needed data and prepare to return
				if done_evaluating:

					# We have finished processing this sample
					break

	###=========================================================================
	### Decide what to do with this sample

	# The position of the kmer that caused us to cut the sample off early
	truncating_duplicate_idx = exceed_duplication_threshold_idx \
								if evaluation_mode == "per_sample_threshold" \
								else first_duplicate_idx
	invalid_sample = truncating_duplicate_idx > -1 and truncating_duplicate_idx < min_sample_len

	# If the offending kmer is within min_sample_len, this sample is invalid
	if invalid_sample:

		# Since this sample is invalid, we can denote its length as -1 and discard its seen kmers
		sample_end_coord = -1
		sample_seen_kmers = None

		# If a local duplicate
		is_local_duplicate = first_duplicate_match_idx > -1
		if is_local_duplicate:

			# The next sample should start at the original kmer's position + 1
			next_start_offset = first_duplicate_match_idx + 1

			# Record the entire region up to and including the original kmer as ignored
			# Importantly, we do not denote the original kmer as a duplicate, since it is
			# both not in the global set and is not added to a sample
			skipped_region = (0, next_start_offset)

		# Else if a match to a global kmer
		else:

			# The next sample should start at the current position + 1
			next_start_offset = first_duplicate_idx + 1

			# Record the current kmer as a duplicate
			duplicate_start_idx = first_duplicate_idx

			# Record the entire region up to but not including the offending kmer as ignored
			if first_duplicate_idx > 0:
				skipped_region = (0, first_duplicate_idx)

	# If truncating_duplicate_idx is > -1, this is a shortened sample
	# Since the sequence before it is accepted, truncating duplicate_idx is now a global repeat
	# If the offending kmer is past min_sample_len, there is still a valid sample with what we have
	# Whether the repeat is internal or global, the current kmer is the offending one
	# ^ This is because the sample is now accepted, meaning the repeat is now global
	elif truncating_duplicate_idx > -1:

		# The final valid kmer is truncating_duplicate_idx - 1, so the sample end is truncating_duplicate_idx - 1 + k
		sample_end_coord = truncating_duplicate_idx - 1 + k 

		# The current kmer is a duplicate
		duplicate_start_idx = truncating_duplicate_idx

		# The next sample should start at 1 + the current position
		next_start_offset = sample_end_coord if no_overlap_samples else truncating_duplicate_idx + 1
		# next_start_offset = truncating_duplicate_idx + 1

	###=========================================================================
	### Get all ignored regions before next start offset

	# If skipped region is defined, we are rejecting the sample
	# In this case, the entirety of what we skip over is ignored, so 
	# skipped region is all-encompassing
	# We will address any other repeats in a future iteration
	if skipped_region is not None:

		ignored_regions = [skipped_region]

	# If ignored region is not defined, either the sample was rejected at idx=0
	# or it was accepted. In these cases, it is possible that some of the 
	# ignored regions before Ns could be before next start offset
	else:
		
		for ignored_start, ignored_end in ignored_regions_before_ambiguous_bases:
			if ignored_start >= next_start_offset:
				break
			ignored_regions.append((ignored_start, min(ignored_end, next_start_offset)))

	###=========================================================================
	### Return all needed data

	return sample_end_coord, duplicate_start_idx, ignored_regions, next_start_offset, sample_seen_kmers


# Function to check a sample for duplicates, allowing some rate of duplicates
def check_sample_retain_info(seq, internal_N_regions, global_seen_kmers, k, dedup_parameter, min_sample_len, no_overlap_samples, evaluation_mode, seq_info):

	###=========================================================================
	### Global data

	# Relevant global variables
	sample_end_coord = len(seq)
	duplicate_start_idx = -1
	ignored_regions = []
	skipped_region = None
	next_start_offset = len(seq) if no_overlap_samples else len(seq) - k + 1

	###=========================================================================
	### Record positions of ambiguous bases

	for N_start, N_end in internal_N_regions:
		seq_info[N_start:N_end] = seq_info.encoding_dict["ambiguous"]

	###=========================================================================
	### Get new contigs after seq_info_offset
	### Also record ignored regions before ambiguous chars

	# Get contigs based on ambiguous chars
	valid_kmer_start_ranges = []
	ignored_regions_before_ambiguous_bases = []
	valid_region_start = seq_info.offset
	for N_start, N_end in internal_N_regions:
		final_kmer_start = N_start - k

		# Record final k-1 bases before N as ignored
		ignored_region_start = max(valid_region_start, final_kmer_start + 1)
		if ignored_region_start < N_start:
			ignored_regions_before_ambiguous_bases.append((ignored_region_start, N_start))

		# Record region before this as valid
		if final_kmer_start >= valid_region_start:
			valid_kmer_start_ranges.append((valid_region_start, final_kmer_start + 1))

		valid_region_start = N_end

	# No need to record anything for final k-1 bases, which should remain unannotated
	# Just record final valid region
	final_kmer_start = len(seq) - k
	if final_kmer_start >= valid_region_start:
		valid_kmer_start_ranges.append((valid_region_start, final_kmer_start + 1))
	
	# Denote all ignored regions
	for ignored_region_start, ignored_region_end in ignored_regions_before_ambiguous_bases:
		seq_info[ignored_region_start:ignored_region_end] = seq_info.encoding_dict["ignored"]

	print(valid_kmer_start_ranges)

	###=========================================================================
	### Accumulate all needed variables regardless of evaluation mode

	# Used by all methods
	done_evaluating = False
	found_first_duplicate = False
	first_duplicate_idx = -1 # Index of first instance of duplicate kmer
	first_duplicate_match_idx = -1 # Index of first instance of local repeat

	# Only for per-sample agnostic mode
	already_decided_agnostic = False

	# Only for per-sample threshold mode
	n_duplicate_kmers = 0 # Currently unused, but could be used in a future version
	n_duplicate_bases = 0
	last_duplicate_end_idx = -1
	n_allowed_full_length_duplicate_bases = int(np.floor(dedup_parameter * len(seq)))
	below_allowed_duplication_threshold = True
	exceed_duplication_threshold_idx = -1

	###=========================================================================
	### Check all contigs

	## Loop through all kmers in this possible sample
	for valid_start, valid_end in valid_kmer_start_ranges:

		if done_evaluating:
			break

		for kmer_start_idx in range(valid_start, valid_end):

			repeat_kmer = False
			duplicate_match_idx = -1

			idx_info = seq_info[kmer_start_idx]
			if idx_info in [
					seq_info.encoding_dict["ambiguous"], # Ns are Ns
					seq_info.encoding_dict["unique"], # already declared distinct from all kmers before it
					seq_info.encoding_dict["ignored"] # see below
				]:
				continue

			elif idx_info == seq_info.encoding_dict["global repeat"]:
				repeat_kmer = True

			# Ignored can never be un-ignored. 3 types of ignored:
			# 	1. K-1 bases before Ns; always consistent
			# 	2. Min_sample_len-1 bases before end of contig; always consistent and not even set in this function
			# 	3. Spans from start of proposed sample to first problem when rejected; 
			# 		always in the past, i.e. never looked at again after being labelled thusly, 
			# 		because next start offset is always the end of the ignored region + 1

			# If we haven't already recorded this as a global kmer, check for it being a repeat
			if not repeat_kmer:

				# Get kmer and its encoding for storage
				kmer = seq[kmer_start_idx:kmer_start_idx+k]
				kmer_num = encode_kmer(kmer)

				# If local repeat
				duplicate_match_idx = seq_info.get_kmer_idx(kmer_num)
				if duplicate_match_idx >= 0: # aka if this kmer is in the sample_seen_kmers
					seq_info[kmer_start_idx] = seq_info.encoding_dict["internal repeat"]
					repeat_kmer = True

				# If global repeat 
				# No need to check global kmer set if this was previously annotated as local repeat 
				# This is because 
				elif idx_info == seq_info.encoding_dict["unannotated"] and kmer_num in global_seen_kmers:
					seq_info[kmer_start_idx] = seq_info.encoding_dict["global repeat"]
					repeat_kmer = True
				
				# If not a repeat
				else:
					seq_info[kmer_start_idx] = seq_info.encoding_dict["unique"]
					seq_info.record_kmer(kmer_num, kmer_start_idx)

			# If this kmer was a repeat, handle it according to the evaluation mode
			if repeat_kmer and not already_decided_agnostic:

				# In per-kmer mode, decide whether to reject this same based on a coin flip with respect to this kmer
				if evaluation_mode == "per_kmer":

					# Continue analyzing kmers, "pretend" this wasn't a duplicate
					if dedup_parameter > 0 and rng.random() < dedup_parameter:
						continue 

					# Done with sample; decide whether to accept or reject based on whether we have reached min_sample_len
					else:
						done_evaluating = True 
						if not found_first_duplicate:
							found_first_duplicate = True
							first_duplicate_idx = kmer_start_idx
							first_duplicate_match_idx = duplicate_match_idx

				# In per-sample agnostic mode, decide on the whole sample based on a coin flip when we reach the first offending kmer
				elif evaluation_mode == "per_sample_agnostic":

					if dedup_parameter > 0 and rng.random() < dedup_parameter:
						already_decided_agnostic = True # Accept sample immediately, but keep analyzing kmers just to record them in sample_seen_kmers
						continue
					else:
						done_evaluating = True # Don't analyze any more kmers, accept or reject now based on whether we have reached min_samaple_len
						if not found_first_duplicate:
							found_first_duplicate = True
							first_duplicate_idx = kmer_start_idx
							first_duplicate_match_idx = duplicate_match_idx

				elif evaluation_mode == "per_sample_threshold":

					# Record possible first duplicate
					if not found_first_duplicate:
						found_first_duplicate = True
						first_duplicate_idx = kmer_start_idx
						first_duplicate_match_idx = duplicate_match_idx

					# Update duplication counts
					n_duplicate_kmers += 1
					if kmer_start_idx >= last_duplicate_end_idx:
						n_duplicate_bases += k
					else:
						n_duplicate_bases += k - (last_duplicate_end_idx - kmer_start_idx)
					last_duplicate_end_idx = kmer_start_idx + k

					# Make a note if we have yet to exceed the allowed rate of duplication
					# This can be used to salvage a sample < len but >= min_sample_len
					if (n_duplicate_bases / (kmer_start_idx+k)) <= dedup_parameter:
						below_allowed_duplication_threshold = True
						exceed_duplication_threshold_idx = -1
					elif below_allowed_duplication_threshold == True:
						below_allowed_duplication_threshold = False
						exceed_duplication_threshold_idx = kmer_start_idx

					# Early exit if we have already exceeded the max allowed duplication
					if n_duplicate_bases > n_allowed_full_length_duplicate_bases:
						done_evaluating = True

			if done_evaluating:
				break

	###=========================================================================
	### Decide what to do with this sample

	# The position of the kmer that caused us to cut the sample off early
	truncating_duplicate_idx = exceed_duplication_threshold_idx \
								if evaluation_mode == "per_sample_threshold" \
								else first_duplicate_idx
	invalid_sample = truncating_duplicate_idx > -1 and truncating_duplicate_idx < min_sample_len

	# If the offending kmer is within min_sample_len, this sample is invalid
	if invalid_sample:

		# Since this sample is invalid, we can denote its length as -1 and discard its seen kmers
		sample_end_coord = -1

		# If a local duplicate
		is_local_duplicate = first_duplicate_match_idx > -1
		if is_local_duplicate:

			# The next sample should start at the original kmer's position + 1
			next_start_offset = first_duplicate_match_idx + 1

			# Record the entire region up to and including the original kmer as ignored
			# Importantly, we do not denote the original kmer as a duplicate, since it is
			# both not in the global set and is not added to a sample
			skipped_region = (0, next_start_offset)

		# Else if a match to a global kmer
		else:

			# The next sample should start at the current position + 1
			next_start_offset = first_duplicate_idx + 1

			# Record the current kmer as a duplicate
			duplicate_start_idx = first_duplicate_idx

			# Record the entire region up to but not including the offending kmer as ignored
			if first_duplicate_idx > 0:
				skipped_region = (0, first_duplicate_idx)

	# If truncating_duplicate_idx is > -1, this is a shortened sample
	# Since the sequence before it is accepted, truncating duplicate_idx is now a global repeat
	# If the offending kmer is past min_sample_len, there is still a valid sample with what we have
	# Whether the repeat is internal or global, the current kmer is the offending one
	# ^ This is because the sample is now accepted, meaning the repeat is now global
	elif truncating_duplicate_idx > -1:

		# The final valid kmer is truncating_duplicate_idx - 1, so the sample end is truncating_duplicate_idx - 1 + k
		sample_end_coord = truncating_duplicate_idx - 1 + k 

		# The current kmer is a duplicate
		duplicate_start_idx = truncating_duplicate_idx

		# The next sample should start at 1 + the current position
		next_start_offset = sample_end_coord if no_overlap_samples else truncating_duplicate_idx + 1

	###=========================================================================
	### Get all ignored regions before next start offset

	# If skipped region is defined, we are rejecting the sample
	# In this case, the entirety of what we skip over is ignored, so 
	# skipped region is all-encompassing
	# We will address any other repeats in a future iteration
	if skipped_region is not None:

		ignored_regions = [skipped_region]

	# If ignored region is not defined, either the sample was rejected at idx=0
	# or it was accepted. In these cases, it is possible that some of the 
	# ignored regions before Ns could be before next start offset
	else:
		
		for ignored_start, ignored_end in ignored_regions_before_ambiguous_bases:
			if ignored_start >= next_start_offset:
				break
			ignored_regions.append((ignored_start, min(ignored_end, next_start_offset)))

	###=========================================================================
	### Return all needed data

	return sample_end_coord, duplicate_start_idx, ignored_regions, next_start_offset


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
	dedup_param = args.dedup_param

	# Instantiate seen_kmers if None
	if seen_kmers is None:
		seen_kmers = set()

	# Set up storage for samples, masked regions, and ignored regions
	sample_regions = []
	masked_starts = []
	ignored_regions = []

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

	# Get positions of all valid kmers and all ambiguous characters
	N_consecutive_allowed_Ns = args.allowed_consecutive_ambiguous_chars 
	valid_contigs, ambiguous_regions = get_contigs_from_chromosome(seq, N_consecutive_allowed_Ns)

	# Analyze all valid contigs, looking at all possible samples within each contig
	for contig_start, contig_end, contig_N_regions in valid_contigs:
   
		# Needed global vars
		sample_start = contig_start
		max_start_idx = contig_end - min_sample_len # final possible start index

		# Investigate every possible sample
		while sample_start <= max_start_idx:

			# Get boundary for this possible sample
			sample_end = min(len(seq), sample_start + sample_len) # Checking against seq len is only necessary for final sample

			###=====================================================================
			### Get internal N regions for this sample
			
			# Create sample_N_regions containing the relative coordinates of Ns within this sample, ranging from 0 to sample_len
			first_overlapping_region_idx = 0
			found_first_overlapping_region = False
			sample_N_regions = []
			for i, (N_start, N_end) in enumerate(contig_N_regions):
				if N_start >= sample_end:
					break
				if N_end > sample_start:
					if not found_first_overlapping_region:
						first_overlapping_region_idx = i
						found_first_overlapping_region = True  
					sample_N_regions.append((max(N_start, sample_start) - sample_start, min(N_end, sample_end) - sample_start))
				
			# To avoid redundant computations, drop any N regions that we have passed
			contig_N_regions = contig_N_regions[first_overlapping_region_idx:]

			###=====================================================================
			### Evaluate this sample

			checked_sample_len, duplicate_start_idx, sample_ignored_regions, next_start_offset, sample_seen_kmers = \
				check_sample(seq[sample_start:sample_end], 
				 				sample_N_regions, 
								seen_kmers, 
								k, 
								dedup_param, 
								min_sample_len, 
								no_overlap_samples, 
								evaluation_method)
			
			print(f"Checked sample from {sample_start} to {sample_end}: checked_sample_len: {checked_sample_len}, duplicate_start_idx: {duplicate_start_idx}, ignored_regions: {sample_ignored_regions}, next_start_offset: {next_start_offset}")
			
			###=====================================================================
			### Update regions with result of evaluation method call

			# Add sample to samples if valid
			if checked_sample_len > -1:
				sample_regions.append((sample_start, sample_start + checked_sample_len))
			
			# Record the duplicate if any was found
			if duplicate_start_idx > -1:
				masked_starts.append(sample_start + duplicate_start_idx)

			# Add ignored regions if any were found
			for ignored_region in sample_ignored_regions:
				ignored_regions.append((sample_start+ignored_region[0], sample_start+ignored_region[1]))

			# Update seen kmers with any sample-specific kmers
			if sample_seen_kmers is not None:
				seen_kmers.update(sample_seen_kmers)

			# Set start index of next iteration of the loop
			sample_start = sample_start + next_start_offset

		###=========================================================================
		### Final housekeeping 

		# Record possible ignored sequence at the end
		if sample_start < contig_end:
			ignored_regions.append((sample_start, contig_end))

	# Convert masked starting indices to regions for more condensed bed files
	masked_regions = condense_masked_regions(masked_starts)

	return sample_regions, masked_regions, ignored_regions, ambiguous_regions, seen_kmers


def deduplicate_seq_retain_info(seq, seen_kmers, args):

	###=========================================================================
	### Set up needed data

	# Collect needed args
	k = args.kmer
	sample_len = args.sample_len
	min_sample_len = args.min_sample_len
	no_overlap_samples = args.no_overlap
	evaluation_method = args.evaluation_method
	dedup_param = args.dedup_param

	# Instantiate seen_kmers if None
	if seen_kmers is None:
		seen_kmers = set()

	# Set up storage for samples, masked regions, and ignored regions
	sample_regions = []
	masked_starts = []
	ignored_regions = []

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

	# Get positions of all valid kmers and all ambiguous characters
	N_consecutive_allowed_Ns = args.allowed_consecutive_ambiguous_chars 
	valid_contigs, ambiguous_regions = get_contigs_from_chromosome(seq, N_consecutive_allowed_Ns)
	# print(valid_contigs)
	# print(ambiguous_regions)

	# Analyze all valid contigs, looking at all possible samples within each contig
	for i, (contig_start, contig_end, contig_N_regions) in enumerate(valid_contigs):

		# Instantiate object to store sample info between iterations for this contig
		seq_info = SeqInfo(sample_len, k)
   
		# Needed global vars
		sample_start = contig_start
		max_start_idx = contig_end - min_sample_len # final possible start index

		# Investigate every possible sample
		while sample_start <= max_start_idx:

			# Get boundary for this possible sample
			sample_end = min(len(seq), sample_start + sample_len) # Checking against seq len is only necessary for final sample

			###=====================================================================
			### Get internal N regions for this sample
			
			# Create sample_N_regions containing the relative coordinates of Ns within this sample, ranging from 0 to sample_len
			first_overlapping_region_idx = 0
			found_first_overlapping_region = False
			sample_N_regions = []
			for i, (N_start, N_end) in enumerate(contig_N_regions):
				if N_start >= sample_end:
					break
				if N_end > sample_start + seq_info.offset:
					if not found_first_overlapping_region:
						first_overlapping_region_idx = i
						found_first_overlapping_region = True  
					sample_N_regions.append((max(N_start, sample_start + seq_info.offset) - sample_start, min(N_end, sample_end) - sample_start))
				
			# To avoid redundant computations, drop any N regions that we have passed
			contig_N_regions = contig_N_regions[first_overlapping_region_idx:]

			###=====================================================================
			### Evaluate this sample

			seq_info.current_global_idx = sample_start
			seq_info.current_sample_length = sample_end - sample_start

			# print(f"BEFORE: CIGAR={seq_info.get_cigar()}; SEEN KMERS={len(seq_info.sample_seen_kmers)}")

			checked_sample_len, duplicate_start_idx, sample_ignored_regions, next_start_offset = \
				check_sample_retain_info(seq[sample_start:sample_end], 
													sample_N_regions, 
													seen_kmers, 
													k, 
													dedup_param, 
													min_sample_len,
													no_overlap_samples, 
													evaluation_method,
													seq_info)

			# print(f"Checked sample from {sample_start} to {sample_end}: checked_sample_len: {checked_sample_len}, duplicate_start_idx: {duplicate_start_idx}, ignored_regions: {sample_ignored_regions}, next_start_offset: {next_start_offset}")
			# print(f"AFTER: CIGAR={seq_info.get_cigar()}; SEEN KMERS={len(seq_info.sample_seen_kmers)}")
			# print("------------------")

			###=====================================================================
			### Update regions with result of evaluation method call

			# Add sample to samples if valid
			if checked_sample_len > -1:
				sample_regions.append((sample_start, sample_start + checked_sample_len))
			
			# Record the duplicate if any was found
			if duplicate_start_idx > -1:
				masked_starts.append(sample_start + duplicate_start_idx)

			# Add ignored region if any was found
			for ignored_region in sample_ignored_regions:
				ignored_regions.append((sample_start+ignored_region[0], sample_start+ignored_region[1]))

			# Update seen kmers with any sample-specific kmers
			if checked_sample_len > -1:
				# print(f"Adding {len(seq_info.sample_seen_kmers)} sample seen kmers to global seen kmers")
				seen_kmers.update(seq_info.sample_seen_kmers)
			# else:
			#     print("Adding 0 sample seen kmers to global seen kmers")

			# Set start index of next iteration of the loop
			sample_start = sample_start + next_start_offset

			# Update seq_info with the number of cleared bases
			seq_info.move_internal_pointer(next_start_offset, last_sample_accepted=(checked_sample_len > -1))

		###=========================================================================
		### Final housekeeping 

		# Record possible ignored sequence at the end
		if sample_start < contig_end:
			ignored_regions.append((sample_start, contig_end))

	# Convert masked starting indices to regions for more condensed bed files
	masked_regions = condense_masked_regions(masked_starts)

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
			if args.retain_info:
				sample_regions, masked_regions, skipped_regions, ambiguous_regions, seen_kmers = deduplicate_seq_retain_info(clean_sequence, seen_kmers, args)
			else:
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

	## Collect input args
	parser = argparse.ArgumentParser()
	parser.add_argument("input", nargs="+", help="Input list of FASTA files or a txt file with one FASTA file per line")
	parser.add_argument("-d", "--dedup_param", type=float, default=0.0, 
						help="Parameter controlling the amount of allowed duplication. Set to 0 for strict deduplication in any mode." + 
						"Per-kmer mode: per-kmer retention rate. Per-sample-agnostic mode: per-sample retention rate. Per-sample-threshold mode: duplicate base % threshold")
	parser.add_argument("-e", "--evaluation_method", type=str, default="per_kmer", choices=["per_kmer", "per_sample_agnostic", "per_sample_threshold"], 
						help="Method for deduplication evaluation (default: per_kmer)")
	parser.add_argument("-k", "--kmer", type=int, default=32, help="Kmer size (default: 32)")
	parser.add_argument("-l", "--sample_len", type=int, default=1000, help="Sample length (default: 1000)")
	parser.add_argument("-m", "--min_sample_len", type=int, default=None, help="Minimum sample length (default: 50)")
	parser.add_argument("-n", "--allowed_consecutive_ambiguous_chars", type=int, default=0, help="Number of allowed consecutive ambiguous characters in a valid sample (default: 0)")
	parser.add_argument("-o", "--output_dir", default="dedup_out", help="Output directory (default: dedup_out/)")
	parser.add_argument("-p", "--seen_kmers", default=None, help="Pickle file containing seen kmers (default: None)")
	parser.add_argument("-r", "--retain_info", action="store_true", help="Whether to retain information as scanning each contig. Recommended for large sample lengths.")
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
	if args.dedup_param < 0.0 or args.dedup_param > 1.0:
		raise("Error: deduplication parameter must be between 0.0 and 1.0")

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
	# print(args)
	if not args.test:
		# print("Running deduplication")
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

