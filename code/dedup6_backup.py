### Deduplication of genomic sequences 

###=============================================================================
### Imports

import argparse
from Bio import SeqIO
import gzip
import json
import numpy as np
import os
import random as rng
import re
import struct
import sys


###=============================================================================
### Helper functions

def decode_kmer(kmer_num, k=32):
	char_map = {0:'A', 1:'C', 2:'G', 3:'T'}
	kmer = []
	for _ in range(k):
		nucleotide_code = kmer_num & 3 # 0x11
		kmer.append(char_map[nucleotide_code])
		kmer_num >>= 2
	return ''.join(reversed(kmer))

def get_repeat_bases_from_cigar(cigar, n_bases):
	# 0 - nothing, 1 - unique (sample), 2 - omitted (ignored), 3 - ambiguous, 4 - local repeat, 5 - global repeat
	repeat_bases = np.zeros(n_bases)
	i = 0
	for count, region_type in re.findall(r'(\d+)([XUOALG])', cigar):
		if region_type in ['L', 'G']:
			repeat_bases[i:i+int(count)+31] = 1
		i += int(count)
		if i >= n_bases - 31:
			break
	return int(repeat_bases.sum())

## Classes ================

class SeqInfo:

	def __init__(self, seqlen, k, current_global_idx=0, encoding_dict=None):

		if encoding_dict is None:
			encoding_dict = {
				"unannotated": 0,
				"unique": 1, 
				"ignored": 2,
				"ambiguous": 3, 
				"internal repeat": 4,
				"global repeat": 5 
			}
		self.encoding_dict = encoding_dict
		self.arr = np.zeros(seqlen, dtype=np.uint8) 
		self.n = seqlen
		self.current_sample_length = seqlen
		self.k = k
		self.current_global_idx = current_global_idx
		self.sample_seen_kmers = {}
		self.n_repeat_bases = 0
		self.last_repeat_end_idx = 0
		self.first_repeat_idx = None

		# Internally, self.arr is treated as a circular array
		# self.array_start is the position in the array that we treat as index 0
		self.array_start = 0


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

	
	def __len__(self):
		return self.n


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


	def get_kmer_idx(self, kmer_num):
		kmer_positions = self.sample_seen_kmers.get(kmer_num, None)
		local_kmer_idx = -1 if kmer_positions is None else kmer_positions[0]
		if local_kmer_idx == -1:
			return -1
		else:
			return local_kmer_idx - self.current_global_idx
		
	
	def get_closest_repeat_idx(self, offset=0, max_distance=None, forward=True, ignore_offset=False):

		# Process input args
		max_position = self.n - self.k if forward else 0
		if max_distance is not None:
			max_position = min(max_position, offset + max_distance) if forward else max(max_position, offset - max_distance)
		increment = 1 if forward else -1
		
		# Loop for closest repeat index
		closest_repeat_idx = None
		i = offset + increment if ignore_offset else offset
		while (i <= max_position) if forward else (i >= max_position):
			if self[i] == self.encoding_dict["global repeat"] or self[i] == self.encoding_dict["internal repeat"]:
				closest_repeat_idx = i
				break
			i += increment
		return closest_repeat_idx


	# Update both the positions and total count of repeat bases
	def update_repeat_bases(self, repeat_idx):
		# if self.last_repeat_idx is None:
		# 	n_new_bases = self.k
		# else:
		# 	n_new_bases = (repeat_idx + self.k) - max(repeat_idx, (self.last_repeat_idx + self.k))
		n_new_bases = (repeat_idx + self.k) - max(repeat_idx, self.last_repeat_end_idx)
		# if self.current_global_idx == 13327501:
		# 	print(f"{repeat_idx}, last_repeat_end_idx = {self.last_repeat_end_idx} - updating repeat base count from {self.n_repeat_bases} to {self.n_repeat_bases + n_new_bases}")
		self.n_repeat_bases += n_new_bases
		# self.last_repeat_idx = repeat_idx
		self.last_repeat_end_idx = repeat_idx + self.k


	def record_local_kmer(self, relative_idx, kmer_num, is_novel):
		if is_novel:
			self.sample_seen_kmers[kmer_num] = []
		self.sample_seen_kmers[kmer_num].append(self.current_global_idx + relative_idx)


	# 1. Record this base as a global repeat
	# 2. Update repeat base count
	# 3. Possibly record this index as the first repeat for this sample
	def record_global_repeat_kmer(self, relative_idx):
		self[relative_idx] = self.encoding_dict["global repeat"]
		self.update_repeat_bases(relative_idx)
		if self.first_repeat_idx is None:
			self.first_repeat_idx = relative_idx


	# 1. Record this base as a local repeat
	# 2. Add this kmer to sample seen kmers
	# 3. Update repeat base count
	# 4. Possibly record the index of the first instance of this local repeat as the first repeat for this sample
	def record_local_repeat_kmer(self, relative_idx, kmer_num, local_match_idx):
		self[relative_idx] = self.encoding_dict["internal repeat"]
		self.record_local_kmer(relative_idx, kmer_num, is_novel=False)
		self.update_repeat_bases(relative_idx)
		if self.first_repeat_idx is None:
			self.first_repeat_idx = local_match_idx 


	# 1. Record this base as unique
	# 2. Add this kmer to sample seen kmers
	def record_unique_kmer(self, relative_idx, kmer_num):
		self[relative_idx] = self.encoding_dict["unique"]
		self.record_local_kmer(relative_idx, kmer_num, is_novel=True)


	# Remove any kmers from sample_seen_kmers that are before the current global index + min_idx.
	# Return the annotation indices that should be reset to unannotated as a result of this purge.
	# Reset indices are the next instance of each kmer that had a position scrubbed,
	# since those indices are no longer local repeats, although each subsequent instance still is.
	# e.g. GATTACA kmer is at indices 0,10,20, and we move the pointer to 5. Now, the kmer at 0 is 
	# purged, the kmer at 10 is unique, and the kmer at 20 is still a local repeat
	# TODO: if len(newly_unique_idxs) >= min_idx can we early exit due to pidgeonhole principle?
	def update_seen_kmers(self, min_idx):
		newly_unique_idxs = []
		if min_idx <= 0:
			return newly_unique_idxs
		if min_idx > self.current_sample_length - self.k:
			self.sample_seen_kmers = {}    
		else:
			# Conveniently, kmer_positions are always in ascending order
			# TODO: this is very slow bc we check all kmers. Find a way to only check the passed kmers (new data structure needed)
			# Could be another array containing the kmer encoding at every base
			kmers = list(self.sample_seen_kmers.keys())
			for kmer in kmers:
				kmer_positions = self.sample_seen_kmers[kmer]
				first_safe_kmer_idx = len(kmer_positions)
				for i, pos in enumerate(kmer_positions):
					if pos >= self.current_global_idx:
						first_safe_kmer_idx = i
						break

				# If all positions of this kmer are before the current global index + min_idx, delete it
				# There are no more indices matching this kmer in the current window
				if first_safe_kmer_idx == len(kmer_positions):
					del self.sample_seen_kmers[kmer]

				elif first_safe_kmer_idx > 0:
					newly_unique_idxs.append(kmer_positions[first_safe_kmer_idx] - self.current_global_idx)
					self.sample_seen_kmers[kmer] = kmer_positions[first_safe_kmer_idx:]
		newly_unique_idxs = sorted(newly_unique_idxs)
		return newly_unique_idxs


	def update_coordinates(self, sample_start, sample_end, seq):

		# print(f"Updating coordinates to sample_start={sample_start} and sample_end={sample_end} with current_global_idx={self.current_global_idx}")

		# Needed global variables for these steps
		n_advanced_bases = sample_start - self.current_global_idx
		passed_kmer_idxs = np.array([]) if n_advanced_bases <= 0 else np.where(
			(self[:n_advanced_bases] == self.encoding_dict["internal repeat"]) | 
			(self[:n_advanced_bases] == self.encoding_dict["global repeat"])
		)[0]

		# 1. Set current sample length
		self.current_sample_length = sample_end - sample_start

		# 2. Set new current global idx
		self.current_global_idx = sample_start

		# 3. Reset old annotations before new global idx
		if n_advanced_bases > 0:
			self[:n_advanced_bases] = self.encoding_dict["unannotated"]

		# 4. Update array start pointer to reflect new global idx
		self.array_start = self._resolve_index(n_advanced_bases)

		# 5. Update seen kmers, removing any before new global idx
		newly_unique_idxs = self.update_seen_kmers(n_advanced_bases)

		# 6. Update annotations of newly unique kmers
		for idx in newly_unique_idxs:
			self[idx] = self.encoding_dict["unique"]

		# Possible easy out for 7/8/9 if we move past all currently annotated kmers
		if n_advanced_bases > self.n - self.k:
			self.n_repeat_bases = 0
			self.first_repeat_idx = None
			self.last_repeat_end_idx = 0
			return

		# 7. Decrement n_repeat_bases based on passed_kmers and newly_unique_idxs
		n_passed_bases = 0
		old_repeat_idxs = [idx - n_advanced_bases for idx in passed_kmer_idxs] + newly_unique_idxs
		if len(old_repeat_idxs) > 0:
			for i in range(len(old_repeat_idxs)):
				curr_removed_kmer_idx = old_repeat_idxs[i]
				last_removed_kmer_idx = old_repeat_idxs[i-1] if i > 0 else None
				last_repeat_idx = None if curr_removed_kmer_idx <= 0 else self.get_closest_repeat_idx(offset=curr_removed_kmer_idx, max_distance=self.k-1, forward=False, ignore_offset=True)
				next_removed_kmer_idx = old_repeat_idxs[i+1] if i < len(old_repeat_idxs) - 1 else None
				next_repeat_idx = None if curr_removed_kmer_idx < -1 else self.get_closest_repeat_idx(offset=curr_removed_kmer_idx, max_distance=self.k-1, forward=True, ignore_offset=True) 
				closest_previous_repeat_idx = max(last_removed_kmer_idx if last_removed_kmer_idx is not None else curr_removed_kmer_idx - self.k, last_repeat_idx if last_repeat_idx is not None else curr_removed_kmer_idx - self.k)
				closest_next_repeat_idx = min(next_removed_kmer_idx if next_removed_kmer_idx is not None else curr_removed_kmer_idx + self.k, next_repeat_idx if next_repeat_idx is not None else curr_removed_kmer_idx + self.k)
				distance_to_previous_repeat = min(self.k, curr_removed_kmer_idx - closest_previous_repeat_idx)
				amount_covered_by_previous_repeat = self.k - distance_to_previous_repeat
				distance_to_next_repeat = min(self.k, closest_next_repeat_idx - curr_removed_kmer_idx)
				n_repeat_bases_unique_to_this_kmer = max(0, distance_to_next_repeat - amount_covered_by_previous_repeat)
				n_passed_bases += n_repeat_bases_unique_to_this_kmer
		self.n_repeat_bases -= n_passed_bases

		# 8. Update first repeat index
		if self.first_repeat_idx is not None:
			self.first_repeat_idx -= n_advanced_bases
			if self.first_repeat_idx < 0:
				new_first_repeat_idx = self.get_closest_repeat_idx()
				if new_first_repeat_idx is None:
					self.first_repeat_idx = None
				elif self[new_first_repeat_idx] == self.encoding_dict["global repeat"]:
					self.first_repeat_idx = new_first_repeat_idx 
				else: # Must be a local repeat, so find the position of its first match to get the correct first repeat index
					repeat_kmer = seq[new_first_repeat_idx:(new_first_repeat_idx+self.k)]
					repeat_kmer_num = encode_kmer(repeat_kmer)
					self.first_repeat_idx = self.get_kmer_idx(repeat_kmer_num) 

		# 9. Update last repeat end index
		self.last_repeat_end_idx = max(0, self.last_repeat_end_idx - n_advanced_bases) 
		if self.last_repeat_end_idx == 0 and self.first_repeat_idx is not None:
			self.last_repeat_end_idx = self.first_repeat_idx + self.k


## Component Functions ================

# Return a "clean" version of the input genomic sequence
# Convert all chars to uppercase and replace ambiguous chars with N
def get_clean_sequence(sequence):
	ambiguous_chars_regex = re.compile(r'[^ACGTN]')
	clean_sequence = ambiguous_chars_regex.sub('N', str(sequence.upper()))
	return clean_sequence


def encode_kmer(kmer):
	char_map = {'A':0, 'C':1, 'G':2, 'T':3}
	kmer_num = 0
	for c in kmer:
		kmer_num = (kmer_num << 2) | char_map[c]
	return kmer_num


def get_fasta_basename(fasta):
	n_suffixes = 2 if fasta.endswith(".gz") else 1
	fasta_basename = '.'.join(os.path.basename(fasta).split('.')[:-(n_suffixes)])
	return fasta_basename


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
def compute_seen_kmers_from_samples_bed_and_fasta(samples_bed, fasta, k):

	# Allocate set for seen kmers
	seen_kmers = set()

	# Read in fasta sequences
	seq_dict = {}
	with open_maybe_gzip(fasta) as f:
		for record in SeqIO.parse(f, "fasta"):
			# Read sequence
			sequence = str(record.seq)
			# Format it correctly
			clean_sequence = get_clean_sequence(sequence)
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


def write_seen_kmers(seen_kmer_set, file_basename):
	outfile = file_basename + ".kmers.bin"
	with open(outfile, "wb") as f:
		for n in seen_kmer_set:
			f.write(struct.pack("<Q", n))


def read_seen_kmers(seen_kmers_file):
	if not os.path.exists(seen_kmers_file):
		raise ValueError(f"Error: seen kmers file {seen_kmers_file} does not exist.")
	with open(seen_kmers_file, "rb") as f:
		data = f.read()
	n = len(data) // 8
	return set(struct.unpack(f"<{n}Q", data))


# Write bed files detailing the label of each kmer in the input fasta
def write_region_beds(kmer_region_annotations, file_prefix, 
					  write_ambiguous=False, write_ignored=False, write_masks=False):

	write_bed_types = ["samples"]
	if write_ambiguous:
		write_bed_types.append("ambiguous")
	if write_ignored:
		write_bed_types.append("ignored")
	if write_masks:
		write_bed_types.append("masks")

	for bed_type in write_bed_types:
		bed_file = file_prefix + f".{bed_type}.bed"
		with open(bed_file, 'w') as f:
			for seqname, annotations in kmer_region_annotations.items():
				regions = annotations.get(bed_type, [])
				for start, end in regions:
					f.write(f"{seqname}\t{start}\t{end}\n")


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
def check_sample(seq, sample_offset, internal_N_regions, global_seen_kmers, k, 
				 dedup_parameter, min_sample_len, overlap, evaluation_mode, seq_info):
	
	# if seq_info.current_global_idx == 182483:
	# 	print(f"first_repeat_idx: {seq_info.first_repeat_idx}")

	###=========================================================================
	### Sequence attributes that will be returned

	sample_end_coord = len(seq)
	duplicate_start_idx = -1
	ignored_regions = []
	skipped_region = None
	next_start_offset = len(seq) - overlap

	# Data structure for new kmers in this sample
	sample_seen_kmers = {} if seq_info is None else None

	# For debugging, save index of furthest progress into the sample
	furthest_progress_idx = len(seq)

	###=========================================================================
	### Record positions of ambiguous bases in retain info mode

	if seq_info is not None:
		for N_start, N_end in internal_N_regions:
			seq_info[N_start:N_end] = seq_info.encoding_dict["ambiguous"]

	###=========================================================================
	### Using the internal N regions, get all ranges of start indices of valid kmers 
	### Also record ignored regions in the k-1 bases before ambiguous chars

	# Get contigs based on ambiguous chars
	valid_kmer_start_ranges = []
	ignored_regions_before_ambiguous_bases = []
	valid_region_start = sample_offset
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
	### Record positions of ignored regions in retain info mode

	if seq_info is not None:
		for ignored_region_start, ignored_region_end in ignored_regions_before_ambiguous_bases:
			seq_info[ignored_region_start:ignored_region_end] = seq_info.encoding_dict["ignored"]
	
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
	n_duplicate_bases = 0 if seq_info is None else seq_info.n_repeat_bases
	last_duplicate_end_idx = -1
	n_allowed_full_length_duplicate_bases = int(np.floor(dedup_parameter * len(seq)))
	below_allowed_duplication_threshold = n_duplicate_bases <= n_allowed_full_length_duplicate_bases
	exceed_duplication_threshold_idx = -1 if below_allowed_duplication_threshold or seq_info is None else max(-1, seq_info.last_repeat_end_idx - k)
	# exceed_duplication_threshold_idx = -1 

	# Possible early exit if we enter the sample already above the duplicate limit
	# This can occur for example if the prior sequence was 868G100U1G31X (932); now 867G100U1G32X (931, already over if threshold is 0.9)
	# if not below_allowed_duplication_threshold:
	# 	valid_kmer_start_ranges = []
	# 	furthest_progress_idx = sample_offset - 1 # We don't actually even check this offset
	# 	if seq_info.current_global_idx >= 182725 and seq_info.current_global_idx <= 182753:
	# 		print(f"Early exiting with n_duplicate_bases = {n_duplicate_bases} and exceed_duplication_threshold_idx = {exceed_duplication_threshold_idx} ({seq_info[exceed_duplication_threshold_idx]})")

	###=========================================================================
	### Check all contigs

	## Loop through all kmers in this possible sample
	for valid_start, valid_end in valid_kmer_start_ranges:

		if done_evaluating:
			break

		for kmer_start_idx in range(valid_start, valid_end):

			is_global_repeat = False
			is_local_repeat = False
			kmer_encoding = -1
			local_match_idx = -1

			## Don't retain info mode
			if seq_info is None:

				# 1. Encode kmer
				kmer = seq[kmer_start_idx:kmer_start_idx+k]
				kmer_encoding = encode_kmer(kmer)

				# 2. Check if local repeat
				local_match_idx = sample_seen_kmers.get(kmer_encoding, -1)
				if local_match_idx >= 0:
					is_local_repeat = True

				# 3. If not local repeat, check if global repeat
				elif kmer_encoding in global_seen_kmers:
					is_global_repeat = True
			
			## Retain info mode
			else:

				# 1. Check if this kmer was previously annotated as ambiguous, unique, or ignored; if so, skip
				idx_info = seq_info[kmer_start_idx]
				# print(f"Index {kmer_start_idx} ({seq[kmer_start_idx:kmer_start_idx+6]}...): {idx_info}") # Debugging print statement
				if idx_info in [
						seq_info.encoding_dict["ambiguous"], # Ns are Ns
						seq_info.encoding_dict["unique"], # already declared distinct from all kmers before it
						seq_info.encoding_dict["ignored"] # if k-1 bases before problem, always true. If before first duplicate, will be skipped
					]:
					continue

				# 2. Check if this kmer was previously annotated as a global repeat; 
				#    if so, make a note so that we don't need to query global kmers
				if idx_info == seq_info.encoding_dict["global repeat"]:
					is_global_repeat = True

				# 3. Check if this kmer was previously annotated as an internal repeat; 
				#     if so, make a note so that we don't need to query local kmers
				elif idx_info == seq_info.encoding_dict["internal repeat"]:
					is_local_repeat = True

				# 4. If we don't already know what this kmer is, query it like normal
				else:

					# 5. Encode kmer
					kmer = seq[kmer_start_idx:kmer_start_idx+k]
					kmer_encoding = encode_kmer(kmer)

					# 6. Check if local repeat
					local_match_idx = seq_info.get_kmer_idx(kmer_encoding)
					if local_match_idx >= 0: # aka if this kmer is in the sample_seen_kmers
						seq_info.record_local_repeat_kmer(kmer_start_idx, kmer_encoding, local_match_idx)
						# seq_info[kmer_start_idx] = seq_info.encoding_dict["internal repeat"]
						# seq_info.record_kmer(kmer_encoding, kmer_start_idx) # Record this instance of the local kmer
						is_local_repeat = True

					# 7. If not local repeat, check if global repeat
					elif idx_info == seq_info.encoding_dict["unannotated"] and kmer_encoding in global_seen_kmers:
						seq_info.record_global_repeat_kmer(kmer_start_idx)
						# seq_info[kmer_start_idx] = seq_info.encoding_dict["global repeat"]
						# seq_info.update_repeat_bases(kmer_start_idx)
						is_global_repeat = True

			# If this kmer is novel, record it in the local kmers
			if not is_global_repeat and not is_local_repeat:
				if seq_info is None:
					sample_seen_kmers[kmer_encoding] = kmer_start_idx
				else:
					seq_info.record_unique_kmer(kmer_start_idx, kmer_encoding) # Record this novel kmer in the sample_seen_kmers and update annotation to unique
					# seq_info.record_kmer(kmer_encoding, kmer_start_idx) # Also update repeat base info
					# seq_info[kmer_start_idx] = seq_info.encoding_dict["unique"]

					# Needed check in per-sample threshold mode for if we are now under the duplicate %
					# Counts of duplicate bases and kmers are still unchanged since the last repeat we saw
					if evaluation_mode == "per_sample_threshold":

						# if seq_info.current_global_idx == 13327501:
						# 	print(f"Unique kmer at {kmer_start_idx}; n_duplicate_bases = {n_duplicate_bases}, below_allowed_duplication_threshold = {below_allowed_duplication_threshold}, exceed_duplication_threshold_idx = {exceed_duplication_threshold_idx}")
					
						# Because we saw a unique kmer, the duplicate % so far can only decrease, meaning we could go from over to under the threshold
						# if not below_allowed_duplication_threshold and (n_duplicate_bases / (kmer_start_idx+k)) <= dedup_parameter:
						if (n_duplicate_bases / (kmer_start_idx+k)) <= dedup_parameter:
							below_allowed_duplication_threshold = True
							exceed_duplication_threshold_idx = -1
							# if seq_info.current_global_idx == 13327501:
							# 	print(f"At {kmer_start_idx}, {n_duplicate_bases} duplicates / {(kmer_start_idx+k)} bases ({(n_duplicate_bases / (kmer_start_idx+k)):.2f}) is below allowed duplication threshold")
					

			# If this kmer was a repeat, handle it according to the evaluation mode
			# Unless we have already decided on the sample's fate in per-sample agnostic mode
			elif not already_decided_agnostic:

				# In per-kmer mode, decide whether to reject this same based on a coin flip with respect to this kmer
				if evaluation_mode == "per_kmer":

					# Continue analyzing kmers, "pretend" this wasn't a duplicate
					if dedup_parameter > 0 and rng.random() < dedup_parameter:
						continue 

					# Done with sample; decide whether to accept or reject based on whether we have reached min_sample_len
					else:
						done_evaluating = True 
						furthest_progress_idx = kmer_start_idx
						if not found_first_duplicate:
							found_first_duplicate = True
							first_duplicate_idx = kmer_start_idx
							first_duplicate_match_idx = local_match_idx

				# In per-sample agnostic mode, decide on the whole sample based on a coin flip when we reach the first offending kmer
				elif evaluation_mode == "per_sample_agnostic":

					# Accept sample immediately, but keep analyzing kmers just to record them in sample_seen_kmers
					if dedup_parameter > 0 and rng.random() < dedup_parameter:
						already_decided_agnostic = True 
						continue

					# Don't analyze any more kmers, accept or reject now based on whether we have reached min_sample_len
					else:
						done_evaluating = True 
						furthest_progress_idx = kmer_start_idx
						if not found_first_duplicate:
							found_first_duplicate = True
							first_duplicate_idx = kmer_start_idx
							first_duplicate_match_idx = local_match_idx

				elif evaluation_mode == "per_sample_threshold":

					# if seq_info.current_global_idx == 182725:
					# 	print(f"Repeat kmer at {kmer_start_idx}; n_duplicate_bases = {n_duplicate_bases}, below_allowed_duplication_threshold = {below_allowed_duplication_threshold}, exceed_duplication_threshold_idx = {exceed_duplication_threshold_idx}")

					# Record possible first duplicate
					if not found_first_duplicate:
						found_first_duplicate = True
						first_duplicate_idx = kmer_start_idx
						first_duplicate_match_idx = local_match_idx

					# Update duplication counts
					if seq_info is None:
						if kmer_start_idx >= last_duplicate_end_idx:
							n_duplicate_bases += k
						else:
							n_duplicate_bases += k - (last_duplicate_end_idx - kmer_start_idx)
						last_duplicate_end_idx = kmer_start_idx + k
					else:
						n_duplicate_bases = seq_info.n_repeat_bases

					# # Make a note if we have yet to exceed the allowed rate of duplication
					# # This can be used to salvage a sample < len but >= min_sample_len
					# if (n_duplicate_bases / (kmer_start_idx+k)) <= dedup_parameter:
					# 	below_allowed_duplication_threshold = True
					# 	exceed_duplication_threshold_idx = -1
					# 	# if seq_info.current_global_idx == 136939:
					# 	# 	print(f"At {kmer_start_idx}, {n_duplicate_bases} duplicates / {(kmer_start_idx+k)} bases ({(n_duplicate_bases / (kmer_start_idx+k)):.2f}) is below allowed duplication threshold")
					# elif below_allowed_duplication_threshold == True:
					# 	below_allowed_duplication_threshold = False
					# 	exceed_duplication_threshold_idx = kmer_start_idx
					# 	# if seq_info.current_global_idx == 136939:
					# 	# 	print(f"At {kmer_start_idx}, {n_duplicate_bases} duplicates / {(kmer_start_idx+k)} bases ({(n_duplicate_bases / (kmer_start_idx+k)):.2f}) is above allowed duplication threshold")
					
					# Because we saw a repeat, the duplicate % so far can only increase, meaning we could go from under to over the threshold
					# if below_allowed_duplication_threshold and (n_duplicate_bases / (kmer_start_idx+k)) > dedup_parameter:
					if (n_duplicate_bases / (kmer_start_idx+k)) > dedup_parameter:
						below_allowed_duplication_threshold = False
						exceed_duplication_threshold_idx = kmer_start_idx
						# if seq_info.current_global_idx == 13327501:
						# 	print(f"At {kmer_start_idx}, {n_duplicate_bases} duplicates / {(kmer_start_idx+k)} bases ({(n_duplicate_bases / (kmer_start_idx+k)):.2f}) is above allowed duplication threshold")
					
					# Early exit if we have already exceeded the max allowed duplication
					if n_duplicate_bases > n_allowed_full_length_duplicate_bases:
						# if seq_info.current_global_idx == 13327501:
						# 	print(f"Exiting early at {kmer_start_idx} due to exceeding duplication threshold ({n_duplicate_bases} > {n_allowed_full_length_duplicate_bases})")
						done_evaluating = True
						furthest_progress_idx = kmer_start_idx
					# else:
					# 	if seq_info.current_global_idx == 13327501:
					# 		print(f"Continuing at {kmer_start_idx}, n_duplicate_bases: {n_duplicate_bases} <= n_allowed_full_length_duplicate_bases: {n_allowed_full_length_duplicate_bases}")

			# If the previous step causes us to reject the sample already, record needed data and prepare to return
			if done_evaluating:
				break

	###=========================================================================
	### Decide what to do with this sample

	if seq_info is not None:
		first_duplicate_idx = -1 if seq_info.first_repeat_idx is None else seq_info.first_repeat_idx
		if first_duplicate_idx == -1:
			first_duplicate_match_idx = -1
		else:
			if seq_info[first_duplicate_idx] == seq_info.encoding_dict["unique"]:
				first_duplicate_match_idx = first_duplicate_idx
			elif seq_info[first_duplicate_idx] == seq_info.encoding_dict["global repeat"]:
				first_duplicate_match_idx = -1
			else:
				print("ERROR: first_duplicate_idx is not unique or global repeat in seq_info. This should never happen.")
				print(f"current_global_idx: {seq_info.current_global_idx}, first_duplicate_idx: {first_duplicate_idx}, seq_info[first_duplicate_idx]: {seq_info[first_duplicate_idx]}, exceed_duplication_threshold_idx: {exceed_duplication_threshold_idx}")
			# first_duplicate_match_idx = -1 if seq_info[first_duplicate_idx] != seq_info.encoding_dict["unique"] else seq_info.first_repeat_idx

	# The position of the kmer that caused us to cut the sample off early
	truncating_duplicate_idx = exceed_duplication_threshold_idx \
								if evaluation_mode == "per_sample_threshold" \
								else first_duplicate_idx
	invalid_sample = truncating_duplicate_idx > -1 and truncating_duplicate_idx < min_sample_len

	if not invalid_sample and n_duplicate_bases > dedup_parameter * len(seq):
		print(f"Warning: accepted sample has {n_duplicate_bases} duplicate bases. exceed_deduplication_threshold_idx = {exceed_duplication_threshold_idx}")

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
			skipped_region = (0, next_start_offset + k - 1)

		# Else if a match to a global kmer
		else:

			# The next sample should start at the current position + 1
			next_start_offset = first_duplicate_idx + 1
			# if seq_info.current_global_idx == 182483:
			# 	print(f"Global duplicate at {first_duplicate_idx}.  first_duplicate_match_idx = {first_duplicate_match_idx} (exceed_duplication_threshold_idx = {exceed_duplication_threshold_idx}), next_start_offset = {next_start_offset}")

			# Record the current kmer as a duplicate
			duplicate_start_idx = first_duplicate_idx

			# Record the entire region up to but not including the offending kmer as ignored
			if first_duplicate_idx > 0:
				skipped_region = (0, first_duplicate_idx + k - 1)

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
		next_start_offset = max(sample_end_coord - overlap, truncating_duplicate_idx + 1)

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

	return sample_end_coord, duplicate_start_idx, ignored_regions, next_start_offset, furthest_progress_idx, sample_seen_kmers


###=============================================================================
### Main functions

def deduplicate_seq(seq, seen_kmers, retain_info, args):

	###=========================================================================
	### Set up needed data

	# Collect needed args
	k = args.kmer
	sample_len = args.sample_len
	min_sample_len = args.min_sample_len
	overlap = args.overlap
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

	# If the seq is shorter than min_sample_len, we can still analyze it if allow-whole-contigs is set
	# This only applies if the whole contig is valid. Otherwise, the sample is invalidated by the Ns policy
	allow_whole_contigs_override = args.allow_whole_contigs and len(seq) < min_sample_len and len(valid_contigs) == 1 

	# print(f"sample_start	sample_offset	furthest_progress_idx	checked_sample_len	next_start_offset	n_repeat_bases	n_derived_repeat_bases	CIGAR")

	# Analyze all valid contigs, looking at all possible samples within each contig
	for contig_start, contig_end, contig_N_regions in valid_contigs:

		print(f"Analyzing contig from {contig_start} to {contig_end}")

		# Instantiate object to store sample info between iterations for this contig
		seq_info = SeqInfo(sample_len, k, contig_start) if retain_info else None
   
		# Needed global vars
		sample_start = contig_start
		max_start_idx = contig_end - min_sample_len # final possible start index
		sample_offset = 0 # Offset within the sample, used if overlap is > k-1

		# Allow the short sample to be analyzed if allow_whole_contigs_override is set
		if allow_whole_contigs_override:
			if max_start_idx < sample_start:
				max_start_idx = sample_start

		# Investigate every possible sample
		while sample_start <= max_start_idx:

			# if sample_start > 10010:
			# 	exit()

			# Get boundary for this possible sample
			sample_end = min(contig_end, sample_start + sample_len) # Checking against seq len is only necessary for final sample

			###=====================================================================
			### Get internal N regions for this sample
			
			# Create sample_N_regions containing the relative coordinates of Ns within this sample, ranging from 0 to sample_len
			first_overlapping_region_idx = 0
			found_first_overlapping_region = False
			sample_N_regions = []
			for i, (N_start, N_end) in enumerate(contig_N_regions):
				if N_start >= sample_end:
					break
				if N_end > sample_start + sample_offset:
					if not found_first_overlapping_region:
						first_overlapping_region_idx = i
						found_first_overlapping_region = True  
					sample_N_regions.append((
						max(N_start, sample_start + sample_offset) - sample_start, 
						min(N_end, sample_end) - sample_start))

			# To avoid redundant computations, drop any N regions that we have passed
			contig_N_regions = contig_N_regions[first_overlapping_region_idx:]

			###=====================================================================
			### Evaluate this sample

			if seq_info is not None:
				seq_info.update_coordinates(sample_start, sample_end, seq[sample_start:sample_end])
				# print(f"Just updated coordinates, n_repeat_bases={seq_info.n_repeat_bases}")

			# print(f"BEFORE: CIGAR={seq_info.get_cigar()}; SEEN KMERS={len(seq_info.sample_seen_kmers)}")
			# print(f"BEFORE: CIGAR={seq_info.get_cigar()}; SEEN KMERS={[f"{decode_kmer(kmer_num)[:6]}..." for kmer_num in seq_info.sample_seen_kmers.keys()][:10]}...") # Print first 10 seen kmers for debugging

			# if seq_info.current_global_idx >= 13327000:
			# 	print(seq[sample_start:(sample_start+100)])

			checked_sample_len, duplicate_start_idx, sample_ignored_regions, next_start_offset, furthest_progress_idx, sample_seen_kmers = \
				check_sample(
					seq[sample_start:sample_end], 
					sample_offset,
					sample_N_regions, 
					seen_kmers, 
					k, 
					dedup_param, 
					min_sample_len,
					min(overlap, (sample_end - sample_start) - 1), # Overlap cannot be longer than sample length - 1
					evaluation_method,
					seq_info
				)
			
			# if seq_info.current_global_idx >= 13_327_000:
			# 	cigar = seq_info.get_cigar(1000) if seq_info is not None else ""
			# 	print(f"{sample_start}	{sample_offset}	{furthest_progress_idx}	{checked_sample_len}	{next_start_offset}	{seq_info.n_repeat_bases}	{cigar}")
			# cigar = seq_info.get_cigar(sample_end - sample_start) if seq_info is not None else ""
			# print(f"{sample_start}	{sample_offset}	{furthest_progress_idx}	{checked_sample_len}	{next_start_offset}	{seq_info.n_repeat_bases}	{get_repeat_bases_from_cigar(cigar, sample_end - sample_start)}	{cigar}")


			# print(f"Checked sample from {sample_start} to {sample_end}: checked_sample_len: {checked_sample_len}, duplicate_start_idx: {duplicate_start_idx}, ignored_regions: {sample_ignored_regions}, next_start_offset: {next_start_offset}")
			# print(f"AFTER: CIGAR={seq_info.get_cigar()}; SEEN KMERS={len(seq_info.sample_seen_kmers)}")
			# print(f"AFTER: CIGAR={seq_info.get_cigar()}; SEEN KMERS={[f"{decode_kmer(kmer_num)[:6]}..." for kmer_num in seq_info.sample_seen_kmers.keys()][:10]}...") # Print first 10 seen kmers for debugging

			###=====================================================================
			### Update regions with result of evaluation method call

			# If we found a valid sample, record it and its seen kmers
			if checked_sample_len > -1:
				sample_regions.append((sample_start, sample_start + checked_sample_len))
				if seq_info is None:
					seen_kmers.update(sample_seen_kmers)
				else:
					seen_kmers.update(seq_info.sample_seen_kmers)
			
			# Record the duplicate if any was found
			if duplicate_start_idx > -1:
				masked_starts.append(sample_start + duplicate_start_idx)

			# Record the ignored regions if any were found
			for ignored_region in sample_ignored_regions:
				ignored_regions.append((sample_start+ignored_region[0], sample_start+ignored_region[1]))

			###=====================================================================
			### Housekeeping for next iteration

			# Set start index of next iteration of the loop
			sample_start = sample_start + next_start_offset

			# How many kmers to ignore in the next iteration
			if checked_sample_len == -1:
				if retain_info:
					sample_offset = (furthest_progress_idx + 1) - next_start_offset
				else:
					sample_offset = 0
			else:
				sample_offset = max(0, sample_end - sample_start - (k-1))

			# if seq_info.current_global_idx >= 13327000:
			# 	print("------------------")

		###=========================================================================
		### Final housekeeping at the end of a contig

		# Record possible ignored sequence at the end
		if sample_start < contig_end:
			ignored_regions.append((sample_start, contig_end))

	# Convert masked starting indices to regions for more condensed bed files
	masked_regions = condense_masked_regions(masked_starts)

	return sample_regions, masked_regions, ignored_regions, ambiguous_regions, seen_kmers


def deduplicate_genome(fasta, seen_kmers, args):

	# Keep dict associating the regions with the sequence name
	kmer_region_annotations = {} 

	# Read in fasta
	with open_maybe_gzip(fasta) as f:
		for record in SeqIO.parse(f, "fasta"):

			# Read in the sequence
			sequence = str(record.seq)
			seqname = record.id

			# Clean this sequence
			clean_sequence = get_clean_sequence(sequence)

			# Deduplicate this sequence
			sample_regions, masked_regions, skipped_regions, ambiguous_regions, seen_kmers = \
				deduplicate_seq(clean_sequence, seen_kmers, args.retain_info, args)

			# print(f"Seen kmer size: {sys.getsizeof(seen_kmers)}")

			# Associate deduplication data with the sequence name
			kmer_region_annotations[seqname] = {
				"samples": sample_regions,
				"masks": masked_regions,
				"ignored": skipped_regions,
				"ambiguous": ambiguous_regions
			}

	return seen_kmers, kmer_region_annotations


def deduplicate(args):
	
	##=========================================================================
	## Process needed files/folders

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
		print(f"Warning: could not find the following fastas: {', '.join(invalid_fastas)}")

	# Write basename to file map for all valid fastas
	basename_fasta_file = os.path.join(args.output_dir, "basename_fasta_match.txt")
	with open(basename_fasta_file, 'w') as f:
		for fasta in valid_fastas:
			fasta_basename = get_fasta_basename(fasta)
			f.write(f"{fasta_basename}\t{fasta}\n")

	##=========================================================================
	## We may be resuming from a previous run, so find the latest saved kmer checkpoint, if any

	# Needed data
	saved_kmers_idxs = []
	needed_beds = ["samples"]
	if args.write_ambiguous_beds:
		needed_beds.append("ambiguous")
	if args.write_ignored_beds:
		needed_beds.append("ignored")
	if args.write_masked_beds:
		needed_beds.append("masks")
	last_beds_idx = -1

	# Loop through the fastas looking for a checkpoint to resume from
	for i,fasta in enumerate(valid_fastas):

		# Check for a saved kmer file
		fasta_basename = get_fasta_basename(fasta)
		saved_kmers_file = os.path.join(args.output_dir, f"{fasta_basename}.kmers.bin")
		if os.path.exists(saved_kmers_file):
			saved_kmers_idxs.append(i)

		# Check for all needed beds
		all_beds_exist = True
		for bed_type in needed_beds:
			bed_file = os.path.join(args.output_dir, f"{fasta_basename}.{bed_type}.bed")
			if not os.path.exists(bed_file):
				all_beds_exist = False
				break
		if all_beds_exist:
			last_beds_idx = i

	# Take the latest saved kmer checkpoint that has the needed corresponding bed files
	last_checkpoint_idx = -1
	for saved_kmers_idx in saved_kmers_idxs:
		if saved_kmers_idx <= last_beds_idx:
			last_checkpoint_idx = saved_kmers_idx

	# Get the latest seen kmer file, if any
	last_checkpoint_kmers_file = None
	if last_checkpoint_idx > -1:
		last_checkpoint_fasta = valid_fastas[last_checkpoint_idx]
		last_checkpoint_basename = get_fasta_basename(last_checkpoint_fasta)
		last_checkpoint_kmers_file = os.path.join(args.output_dir, f"{last_checkpoint_basename}.kmers.bin")

	# Resume from the next fasta after the last saved checkpoint, or start from the beginning if no checkpoints were found
	next_fasta_idx = last_checkpoint_idx + 1
	process_fastas = valid_fastas[next_fasta_idx:]

	##=========================================================================
	## Initialize seen kmers here
	## 5 options: either or both of args.seen_kmers and last_checkpoint_kmers_file could be set;
	## If both are set, they could either be the same file or different files
	## Only difficult case is if both are set and different, 
	## in which case we will take the union of the two sets of kmers and issue a warning about it

	if args.seen_kmers is None:

		# Both are None, so start with empty set
		if last_checkpoint_kmers_file is None:
			seen_kmers = set()

		# Only checkpoint kmers file is available, so use it
		else:
			seen_kmers = read_seen_kmers(last_checkpoint_kmers_file)
			
	else:

		# Only seen kmers file is available, so use it
		if last_checkpoint_kmers_file is None:
			seen_kmers = read_seen_kmers(args.seen_kmers)

		# Both are available
		else:

			# They are the same file, so no problem
			if os.path.abspath(args.seen_kmers) == os.path.abspath(last_checkpoint_kmers_file):
				seen_kmers = read_seen_kmers(args.seen_kmers)

			# They are different files, so take the union and issue a warning
			else:
				print("Warning: a seen kmers file was supplied, but it does not match the latest checkpoint kmer set. " + 
		  				"Taking the union of the two sets of kmers and proceeding with deduplication.")
				seen_kmers_1 = read_seen_kmers(args.seen_kmers)
				seen_kmers_2 = read_seen_kmers(last_checkpoint_kmers_file)
				seen_kmers = seen_kmers_1.union(seen_kmers_2)

	##=========================================================================
	## Perform deduplication on each fasta

	# Iterate over fastas and deduplicate each one
	for i,fasta in enumerate(process_fastas):

		# Figure out whether we need to deduplicate this sample or just load its existing info
		fasta_basename = get_fasta_basename(fasta)
		all_beds_exist = True
		for bed_type in needed_beds:
			bed_file = os.path.join(args.output_dir, f"{fasta_basename}.{bed_type}.bed")
			if not os.path.exists(bed_file):
				all_beds_exist = False
				break

		out_prefix = os.path.join(args.output_dir, fasta_basename)

		# If we already completed this file, compute its kmer set from its samples
		# This is the case where we save kmers every nth fasta, and this is past that but still completed
		if all_beds_exist:
			print(f"Found existing output for {fasta}, skipping deduplication and computing seen kmers from existing samples")
			samples_bed = os.path.join(args.output_dir, f"{fasta_basename}.samples.bed")
			fasta_seen_kmers = compute_seen_kmers_from_samples_bed_and_fasta(samples_bed, fasta, args.kmer)
			seen_kmers.update(fasta_seen_kmers)

		else:

			# Otherwise, deduplicate as normal
			print(f"Deduplicating {fasta}")
			seen_kmers, kmer_region_annotations = deduplicate_genome(fasta, seen_kmers, args)

			# Write bed files for this fasta
			write_region_beds(kmer_region_annotations, out_prefix, 
					 args.write_ambiguous_beds, args.write_ignored_beds, args.write_masked_beds)

		# Optionally save seen kmers after this fasta
		save_kmers_to_file = args.save_every > 0 and (next_fasta_idx+i+1) % args.save_every == 0
		if save_kmers_to_file:
			write_seen_kmers(seen_kmers, out_prefix)

	# Optionally save seen kmers at the end
	if args.save_kmers_at_end:
		final_seen_kmer_file_basename = "final"
		idx_suffix = 1
		while os.path.exists(os.path.join(args.output_dir, f"{final_seen_kmer_file_basename}.kmers.bin")):
			final_seen_kmer_file_basename = f"final_{idx_suffix}"
			idx_suffix += 1
		write_seen_kmers(seen_kmers, os.path.join(args.output_dir, final_seen_kmer_file_basename))


###=============================================================================
### Call to main

def __main__():

	## Collect input args
	parser = argparse.ArgumentParser()
	parser.add_argument("input", nargs="+", help="Input list of FASTA files or a txt file with one FASTA file per line")
	parser.add_argument("--allow_whole_contigs", action="store_true", help="Allow whole contigs as samples, even if they are shorter than min_sample_len")
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
	parser.add_argument("-v", "--overlap", type=int, default=None, help="Overlap between samples (default: k-1)")
	parser.add_argument("--save_kmers_at_end", action="store_true", help="Save the seen kmer set at the end of program execution")
	parser.add_argument("--write_ambiguous_beds", action="store_true", help="Save a bed file of ambiguous regions for each input fasta")
	parser.add_argument("--write_ignored_beds", action="store_true", help="Save a bed file of ignored regions for each input fasta")
	parser.add_argument("--write_masked_beds", action="store_true", help="Save a bed file of masked regions for each input fasta")
	parser.add_argument("-seed", "--seed", type=int, default=123, help="Random seed for reproducibility")
	args = parser.parse_args()

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
		print("Warning: small kmer sizes will result in very strict deduplication. Consider increasing the kmer size (default: 32)")
	if args.kmer > 32:
		print("Warning: kmer sizes above 32 may not work or may result in very slow evaluation. Consider decreasing the kmer size (default: 32)")

	# Check parameters related to evaluation method
	if args.evaluation_method not in ["per_kmer", "per_sample_agnostic", "per_sample_threshold"]:
		raise("Error: evaluation method must be one of per_kmer, per_sample_agnostic, or per_sample_threshold")
	if args.dedup_param < 0.0 or args.dedup_param > 1.0:
		raise("Error: deduplication parameter must be between 0.0 and 1.0")

	# Set min kmer length to sample length if None
	if args.min_sample_len is None:
		args.min_sample_len = args.sample_len

	# Set overlap to k-1 if None
	if args.overlap is None:
		args.overlap = args.kmer - 1
	
	# Overlap must be smaller than sample length
	if args.overlap >= args.sample_len:
		raise("Error: sample overlap must be smaller than sample length")

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

	# Set random seed for reproducibility
	rng.seed(args.seed)

	## Run deduplication
	deduplicate(args)

if __name__ == "__main__":
	__main__()

