
/**
 * Faithful C++ translation of dedup5.py
 * Matches the Python script's behavior as closely as possible.
 */

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <random>
#include <tuple>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <zlib.h>

namespace fs = std::filesystem;

struct Region {
	size_t start;
	size_t end;
};

using AnnotationMap = std::unordered_map<std::string, std::vector<Region>>;
using GenomeAnnotations = std::vector<std::pair<std::string, AnnotationMap>>;

struct Args {
	std::vector<std::string> input;
	bool allow_whole_contigs = false;
	double dedup_param = 0.0;
	std::string evaluation_method = "per_kmer";
	int kmer = 32;
	int sample_len = 1000;
	int min_sample_len = -1;
	int allowed_consecutive_ambiguous_chars = 0;
	std::string output_dir = "dedup_out";
	std::optional<std::string> seen_kmers;
	bool retain_info = false;
	int save_every = 0;
	std::optional<int> overlap;
	bool save_kmers_at_end = false;
	bool write_ambiguous_beds = false;
	bool write_ignored_beds = false;
	bool write_masked_beds = false;
	int seed = 123;
};

static std::mt19937 rng_gen;

static bool ends_with(const std::string& s, const std::string& suffix) {
	return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

static bool is_gzip_file(const std::string& fname) {
	std::ifstream f(fname, std::ios::binary);
	if (!f) {
		return false;
	}
	unsigned char sig[2] = {0, 0};
	f.read(reinterpret_cast<char*>(sig), 2);
	return f.gcount() == 2 && sig[0] == 0x1f && sig[1] == 0x8b;
}

static std::vector<std::string> read_text_lines(const std::string& fname) {
	std::vector<std::string> lines;
	if (is_gzip_file(fname)) {
		gzFile gz = gzopen(fname.c_str(), "rb");
		if (!gz) {
			throw std::runtime_error("Cannot open file: " + fname);
		}
		std::string line;
		int ch;
		while ((ch = gzgetc(gz)) != -1) {
			if (ch == '\n') {
				if (!line.empty() && line.back() == '\r') {
					line.pop_back();
				}
				lines.push_back(line);
				line.clear();
			} else {
				line.push_back(static_cast<char>(ch));
			}
		}
		if (!line.empty()) {
			if (!line.empty() && line.back() == '\r') {
				line.pop_back();
			}
			lines.push_back(line);
		}
		gzclose(gz);
	} else {
		std::ifstream f(fname);
		if (!f) {
			throw std::runtime_error("Cannot open file: " + fname);
		}
		std::string line;
		while (std::getline(f, line)) {
			if (!line.empty() && line.back() == '\r') {
				line.pop_back();
			}
			lines.push_back(line);
		}
	}
	return lines;
}

static std::string get_fasta_basename(const std::string& fasta) {
	auto basename = fs::path(fasta).filename().string();
	int n_suffixes = ends_with(fasta, ".gz") ? 2 : 1;
	std::vector<std::string> parts;
	std::stringstream ss(basename);
	std::string part;
	while (std::getline(ss, part, '.')) {
		parts.push_back(part);
	}
	if (static_cast<int>(parts.size()) <= n_suffixes) {
		return basename;
	}
	std::string out = parts[0];
	for (size_t i = 1; i + static_cast<size_t>(n_suffixes) < parts.size(); ++i) {
		out += "." + parts[i];
	}
	return out;
}

static std::string get_clean_sequence(const std::string& sequence) {
	std::string out;
	out.reserve(sequence.size());
	for (char c : sequence) {
		char u = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
		if (u == 'A' || u == 'C' || u == 'G' || u == 'T' || u == 'N') {
			out.push_back(u);
		} else {
			out.push_back('N');
		}
	}
	return out;
}

static uint64_t encode_kmer(const std::string& kmer) {
	uint64_t v = 0;
	for (char c : kmer) {
		v <<= 2;
		switch (c) {
			case 'A': break;
			case 'C': v |= 1; break;
			case 'G': v |= 2; break;
			case 'T': v |= 3; break;
			default: throw std::runtime_error("Invalid nucleotide in kmer");
		}
	}
	return v;
}

static std::vector<Region> condense_masked_regions(const std::vector<size_t>& masked) {
	std::vector<Region> masked_regions;
	if (masked.empty()) {
		return masked_regions;
	}
	size_t region_start = masked[0];
	for (size_t i = 0; i + 1 < masked.size(); ++i) {
		if (masked[i] + 1 != masked[i + 1]) {
			masked_regions.push_back({region_start, masked[i] + 1});
			region_start = masked[i + 1];
		}
	}
	masked_regions.push_back({region_start, masked.back() + 1});
	return masked_regions;
}

static std::string get_cigar(const std::vector<uint8_t>& seq_info) {
	static const std::vector<std::string> region_type_char_codes = {"X", "U", "O", "A", "L", "G"};
	std::string cigar;
	std::optional<uint8_t> curr_type;
	size_t curr_start_idx = 0;
	for (size_t i = 0; i < seq_info.size(); ++i) {
		if (!curr_type.has_value() || seq_info[i] != curr_type.value()) {
			if (curr_type.has_value()) {
				cigar += std::to_string(i - curr_start_idx) + region_type_char_codes[curr_type.value()];
			}
			curr_type = seq_info[i];
			curr_start_idx = i;
		}
	}
	if (curr_type.has_value()) {
		cigar += std::to_string(seq_info.size() - curr_start_idx) + region_type_char_codes[curr_type.value()];
	}
	return cigar;
}

static std::vector<Region> regions_from_ranges(const std::vector<std::pair<size_t, size_t>>& ranges) {
	std::vector<Region> out;
	for (const auto& r : ranges) {
		out.push_back({r.first, r.second});
	}
	return out;
}

static void type_check(const std::string& file) {
	const std::vector<std::string> allowed_suffixes = {".gz", ".fasta", ".fa", ".fna", ".txt", ".list"};
	bool good_suffix = false;
	for (const auto& suffix : allowed_suffixes) {
		if (ends_with(file, suffix)) {
			good_suffix = true;
			break;
		}
	}
	if (!good_suffix) {
		std::string msg = "Error: could not determine file type. Supported types are ";
		for (size_t i = 0; i < allowed_suffixes.size(); ++i) {
			if (i) msg += ", ";
			msg += allowed_suffixes[i];
		}
		throw std::runtime_error(msg);
	}
}

static void write_le64(std::ofstream& out, uint64_t v) {
	for (int i = 0; i < 8; ++i) {
		unsigned char b = static_cast<unsigned char>((v >> (8 * i)) & 0xff);
		out.write(reinterpret_cast<const char*>(&b), 1);
	}
}

static uint64_t read_le64(const unsigned char* p) {
	uint64_t v = 0;
	for (int i = 7; i >= 0; --i) {
		v <<= 8;
		v |= static_cast<uint64_t>(p[i]);
	}
	return v;
}

static void write_seen_kmers(const std::unordered_set<uint64_t>& seen_kmer_set, const std::string& file_basename) {
	std::string outfile = file_basename + ".kmers.bin";
	std::ofstream f(outfile, std::ios::binary);
	for (uint64_t n : seen_kmer_set) {
		write_le64(f, n);
	}
}

static std::unordered_set<uint64_t> read_seen_kmers(const std::string& seen_kmers_file) {
	if (!fs::exists(seen_kmers_file)) {
		throw std::runtime_error("Error: seen kmers file " + seen_kmers_file + " does not exist.");
	}
	std::ifstream f(seen_kmers_file, std::ios::binary);
	if (!f) {
		throw std::runtime_error("Error: seen kmers file " + seen_kmers_file + " does not exist.");
	}
	f.seekg(0, std::ios::end);
	std::streamsize size = f.tellg();
	f.seekg(0, std::ios::beg);
	std::vector<unsigned char> data(static_cast<size_t>(size));
	if (size > 0) {
		f.read(reinterpret_cast<char*>(data.data()), size);
	}
	if (data.size() % 8 != 0) {
		throw std::runtime_error("Error: seen kmers file has invalid length.");
	}
	std::unordered_set<uint64_t> out;
	for (size_t i = 0; i < data.size(); i += 8) {
		out.insert(read_le64(&data[i]));
	}
	return out;
}

static std::vector<std::pair<std::string, std::string>> read_fasta_records(const std::string& fasta) {
	std::vector<std::pair<std::string, std::string>> records;
	std::vector<std::string> lines = read_text_lines(fasta);
	std::string current_name;
	std::string current_seq;
	for (const auto& line : lines) {
		if (line.empty()) {
			continue;
		}
		if (line[0] == '>') {
			if (!current_name.empty()) {
				records.push_back({current_name, current_seq});
				current_seq.clear();
			}
			current_name = line.substr(1);
			auto space_pos = current_name.find(' ');
			if (space_pos != std::string::npos) {
				current_name = current_name.substr(0, space_pos);
			}
		} else if (!current_name.empty()) {
			current_seq += line;
		}
	}
	if (!current_name.empty()) {
		records.push_back({current_name, current_seq});
	}
	return records;
}

static std::unordered_set<uint64_t> compute_seen_kmers_from_samples_bed_and_fasta(const std::string& samples_bed, const std::string& fasta, int k) {
	std::unordered_set<uint64_t> seen_kmers;
	std::unordered_map<std::string, std::string> seq_dict;
	for (const auto& rec : read_fasta_records(fasta)) {
		seq_dict[rec.first] = get_clean_sequence(rec.second);
	}
	for (const auto& line : read_text_lines(samples_bed)) {
		if (line.empty()) {
			continue;
		}
		std::vector<std::string> fields;
		std::stringstream ss(line);
		std::string field;
		while (std::getline(ss, field, '\t')) {
			fields.push_back(field);
		}
		if (fields.size() < 3) {
			continue;
		}
		const std::string& seqname = fields[0];
		int start = std::stoi(fields[1]);
		int end = std::stoi(fields[2]);
		auto it = seq_dict.find(seqname);
		if (it == seq_dict.end()) {
			throw std::runtime_error("Error: sequence " + seqname + " found in bed file but not in fasta file.");
		}
		const std::string& sequence = it->second;
		for (int i = start; i <= end - k; ++i) {
			std::string kmer = sequence.substr(static_cast<size_t>(i), static_cast<size_t>(k));
			if (kmer.find('N') == std::string::npos) {
				seen_kmers.insert(encode_kmer(kmer));
			}
		}
	}
	return seen_kmers;
}

struct Contig {
	size_t start = 0;
	size_t end = 0;
	std::vector<std::pair<size_t, size_t>> internal_N_regions;
};

static std::pair<std::vector<Contig>, std::vector<std::pair<size_t, size_t>>> get_contigs_from_chromosome(const std::string& seq, int N_consecutive_allowed_Ns) {
	std::vector<size_t> N_idxs;
	for (size_t i = 0; i < seq.size(); ++i) {
		if (seq[i] == 'N') {
			N_idxs.push_back(i);
		}
	}
	std::vector<std::pair<size_t, size_t>> N_regions;
	if (!N_idxs.empty()) {
		size_t region_start = N_idxs[0];
		for (size_t i = 0; i + 1 < N_idxs.size(); ++i) {
			if (N_idxs[i] + 1 != N_idxs[i + 1]) {
				N_regions.push_back({region_start, N_idxs[i] + 1});
				region_start = N_idxs[i + 1];
			}
		}
		N_regions.push_back({region_start, N_idxs.back() + 1});
	}

	std::vector<Contig> contigs;
	size_t start_idx = 0;
	std::vector<std::pair<size_t, size_t>> internal_N_regions;
	for (const auto& nr : N_regions) {
		size_t N_start = nr.first;
		size_t N_end = nr.second;
		size_t N_region_length = N_end - N_start;
		if (static_cast<int>(N_region_length) > N_consecutive_allowed_Ns) {
			if (start_idx < N_start) {
				contigs.push_back({start_idx, N_start, internal_N_regions});
			}
			start_idx = N_end;
			internal_N_regions.clear();
		} else {
			internal_N_regions.push_back({N_start, N_end});
		}
	}
	if (start_idx < seq.size()) {
		contigs.push_back({start_idx, seq.size(), internal_N_regions});
	}
	return {contigs, N_regions};
}

class SeqInfo {
public:
	std::unordered_map<std::string, int> encoding_dict;
	std::vector<uint8_t> arr;
	size_t n;
	size_t current_sample_length;
	int k;
	size_t current_global_idx;
	std::unordered_map<uint64_t, size_t> sample_seen_kmers;
	size_t array_start;

	SeqInfo(size_t seqlen, int k_, std::optional<std::unordered_map<std::string, int>> encoding_dict_in = std::nullopt)
		: n(seqlen), current_sample_length(seqlen), k(k_), current_global_idx(0), array_start(0) {
		if (encoding_dict_in.has_value()) {
			encoding_dict = encoding_dict_in.value();
		} else {
			encoding_dict = {
				{"unannotated", 0},
				{"unique", 1},
				{"ignored", 2},
				{"ambiguous", 3},
				{"internal repeat", 4},
				{"global repeat", 5}
			};
		}
		arr.assign(seqlen, 0);
	}

	void record_kmer(uint64_t kmer_num, size_t relative_idx) {
		sample_seen_kmers[kmer_num] = current_global_idx + relative_idx;
	}

	long get_kmer_idx(uint64_t kmer_num) const {
		auto it = sample_seen_kmers.find(kmer_num);
		if (it == sample_seen_kmers.end()) {
			return -1;
		}
		return static_cast<long>(it->second - current_global_idx);
	}

	void reset_annotations(size_t n_bases) {
		set_range(0, n_bases, static_cast<uint8_t>(encoding_dict.at("unannotated")));
	}

	void purge_kmers(std::optional<size_t> min_idx = std::nullopt) {
		if (!min_idx.has_value()) {
			sample_seen_kmers.clear();
			return;
		}
		if (min_idx.value() > std::numeric_limits<size_t>::max()) {
			throw std::runtime_error("Error: min_idx must be an integer, but got invalid value");
		}
		for (auto it = sample_seen_kmers.begin(); it != sample_seen_kmers.end(); ) {
			if (it->second - current_global_idx < min_idx.value()) {
				it = sample_seen_kmers.erase(it);
			} else {
				++it;
			}
		}
	}

	size_t resolve_index(size_t idx) const {
		if (idx >= n) {
			throw std::out_of_range("index " + std::to_string(idx) + " is out of bounds for size " + std::to_string(n));
		}
		size_t adj_idx = array_start + idx;
		if (adj_idx >= n) {
			adj_idx -= n;
		}
		return adj_idx;
	}

	std::vector<size_t> resolve_slice(size_t start, size_t stop) const {
		if (start > n) start = n;
		if (stop > n) stop = n;
		std::vector<size_t> out;
		for (size_t i = start; i < stop; ++i) {
			out.push_back(resolve_index(i));
		}
		return out;
	}

	uint8_t get(size_t key) const {
		return arr[resolve_index(key)];
	}

	std::vector<uint8_t> get_range(size_t start, size_t stop) const {
		std::vector<uint8_t> out;
		for (size_t idx : resolve_slice(start, stop)) {
			out.push_back(arr[idx]);
		}
		return out;
	}

	void set(size_t key, uint8_t value) {
		arr[resolve_index(key)] = value;
	}

	void set_range(size_t start, size_t stop, uint8_t value) {
		for (size_t idx : resolve_slice(start, stop)) {
			arr[idx] = value;
		}
	}

	void move_internal_pointer(size_t n_bases, bool last_sample_accepted = false) {
		reset_annotations(n_bases);
		if (last_sample_accepted) {
			purge_kmers();
		} else {
			purge_kmers(n_bases);
		}
		array_start = resolve_index(n_bases);
	}

	std::string get_cigar(std::optional<size_t> n_bases = std::nullopt) const {
		static const std::vector<std::string> region_type_char_codes = {"X", "U", "O", "A", "L", "G"};
		std::vector<uint8_t> cigar_arr = get_range(0, n_bases.has_value() ? std::min(n_bases.value(), current_sample_length) : current_sample_length);
		std::string cigar;
		std::optional<uint8_t> curr_type;
		size_t curr_start_idx = 0;
		for (size_t i = 0; i < cigar_arr.size(); ++i) {
			if (!curr_type.has_value() || cigar_arr[i] != curr_type.value()) {
				if (curr_type.has_value()) {
					cigar += std::to_string(i - curr_start_idx) + region_type_char_codes[curr_type.value()];
				}
				curr_type = cigar_arr[i];
				curr_start_idx = i;
			}
		}
		if (curr_type.has_value()) {
			cigar += std::to_string(this->n - curr_start_idx) + region_type_char_codes[curr_type.value()];
		}
		return cigar;
	}

	size_t size() const {
		return n;
	}
};

struct CheckSampleResult {
	int sample_end_coord = 0;
	int duplicate_start_idx = -1;
	std::vector<std::pair<int, int>> ignored_regions;
	int next_start_offset = 0;
	std::optional<std::unordered_map<uint64_t, size_t>> sample_seen_kmers;
};

static CheckSampleResult check_sample(
	const std::string& seq,
	size_t sample_offset,
	const std::vector<std::pair<size_t, size_t>>& internal_N_regions,
	const std::unordered_set<uint64_t>& global_seen_kmers,
	int k,
	double dedup_parameter,
	int min_sample_len,
	int overlap,
	const std::string& evaluation_mode,
	SeqInfo* seq_info
) {
	CheckSampleResult result;
	result.sample_end_coord = static_cast<int>(seq.size());
	result.duplicate_start_idx = -1;
	result.next_start_offset = static_cast<int>(seq.size()) - overlap;
	result.sample_seen_kmers = seq_info == nullptr ? std::optional<std::unordered_map<uint64_t, size_t>>(std::unordered_map<uint64_t, size_t>{}) : std::nullopt;

	std::vector<std::pair<size_t, size_t>> valid_kmer_start_ranges;
	std::vector<std::pair<int, int>> ignored_regions_before_ambiguous_bases;
	if (seq_info != nullptr) {
		for (const auto& nr : internal_N_regions) {
			seq_info->set_range(nr.first, nr.second, static_cast<uint8_t>(seq_info->encoding_dict["ambiguous"]));
		}
	}

	long long valid_region_start = static_cast<long long>(sample_offset);
	for (const auto& nr : internal_N_regions) {
		long long N_start = static_cast<long long>(nr.first);
		long long N_end = static_cast<long long>(nr.second);
		long long final_kmer_start = N_start - k;
		long long ignored_region_start = std::max(valid_region_start, final_kmer_start + 1);
		if (ignored_region_start < N_start) {
			ignored_regions_before_ambiguous_bases.push_back({static_cast<int>(ignored_region_start), static_cast<int>(N_start)});
		}
		if (final_kmer_start >= valid_region_start) {
			valid_kmer_start_ranges.push_back({static_cast<size_t>(valid_region_start), static_cast<size_t>(final_kmer_start + 1)});
			valid_region_start = N_end;
		}
	}
	long long final_kmer_start = static_cast<long long>(seq.size()) - k;
	if (final_kmer_start >= valid_region_start) {
		valid_kmer_start_ranges.push_back({static_cast<size_t>(valid_region_start), static_cast<size_t>(final_kmer_start + 1)});
	}

	if (seq_info != nullptr) {
		for (const auto& ir : ignored_regions_before_ambiguous_bases) {
			seq_info->set_range(static_cast<size_t>(ir.first), static_cast<size_t>(ir.second), static_cast<uint8_t>(seq_info->encoding_dict["ignored"]));
		}
	}

	bool done_evaluating = false;
	bool found_first_duplicate = false;
	int first_duplicate_idx = -1;
	int first_duplicate_match_idx = -1;
	bool already_decided_agnostic = false;
	int n_duplicate_kmers = 0;
	int n_duplicate_bases = 0;
	int last_duplicate_end_idx = -1;
	int n_allowed_full_length_duplicate_bases = static_cast<int>(std::floor(dedup_parameter * static_cast<double>(seq.size())));
	bool below_allowed_duplication_threshold = true;
	int exceed_duplication_threshold_idx = -1;

	for (const auto& range : valid_kmer_start_ranges) {
		if (done_evaluating) {
			break;
		}
		for (size_t kmer_start_idx = range.first; kmer_start_idx < range.second; ++kmer_start_idx) {
			bool is_global_repeat = false;
			bool is_local_repeat = false;
			uint64_t kmer_encoding = static_cast<uint64_t>(-1);
			long local_match_idx = -1;

			if (seq_info == nullptr) {
				std::string kmer = seq.substr(kmer_start_idx, k);
				kmer_encoding = encode_kmer(kmer);
				auto it = result.sample_seen_kmers->find(kmer_encoding);
				if (it != result.sample_seen_kmers->end()) {
					local_match_idx = static_cast<long>(it->second);
					is_local_repeat = true;
				} else if (global_seen_kmers.find(kmer_encoding) != global_seen_kmers.end()) {
					is_global_repeat = true;
				}
			} else {
				uint8_t idx_info = seq_info->get(kmer_start_idx);
				if (idx_info == seq_info->encoding_dict["ambiguous"] || idx_info == seq_info->encoding_dict["unique"] || idx_info == seq_info->encoding_dict["ignored"]) {
					continue;
				}
				if (idx_info == seq_info->encoding_dict["global repeat"]) {
					is_global_repeat = true;
				} else {
					std::string kmer = seq.substr(kmer_start_idx, k);
					kmer_encoding = encode_kmer(kmer);
					local_match_idx = seq_info->get_kmer_idx(kmer_encoding);
					if (local_match_idx >= 0) {
						seq_info->set(kmer_start_idx, static_cast<uint8_t>(seq_info->encoding_dict["internal repeat"]));
						is_local_repeat = true;
					} else if (idx_info == seq_info->encoding_dict["unannotated"] && global_seen_kmers.find(kmer_encoding) != global_seen_kmers.end()) {
						seq_info->set(kmer_start_idx, static_cast<uint8_t>(seq_info->encoding_dict["global repeat"]));
						is_global_repeat = true;
					}
				}
			}

			if (!is_global_repeat && !is_local_repeat) {
				if (seq_info == nullptr) {
					(*result.sample_seen_kmers)[kmer_encoding] = kmer_start_idx;
				} else {
					seq_info->record_kmer(kmer_encoding, kmer_start_idx);
					seq_info->set(kmer_start_idx, static_cast<uint8_t>(seq_info->encoding_dict["unique"]));
				}
			} else if (!already_decided_agnostic) {
				if (evaluation_mode == "per_kmer") {
					if (dedup_parameter > 0.0) {
						std::uniform_real_distribution<double> dist(0.0, 1.0);
						if (dist(rng_gen) < dedup_parameter) {
							continue;
						}
					}
					done_evaluating = true;
					if (!found_first_duplicate) {
						found_first_duplicate = true;
						first_duplicate_idx = static_cast<int>(kmer_start_idx);
						first_duplicate_match_idx = static_cast<int>(local_match_idx);
					}
				} else if (evaluation_mode == "per_sample_agnostic") {
					if (dedup_parameter > 0.0) {
						std::uniform_real_distribution<double> dist(0.0, 1.0);
						if (dist(rng_gen) < dedup_parameter) {
							already_decided_agnostic = true;
							continue;
						}
					}
					done_evaluating = true;
					if (!found_first_duplicate) {
						found_first_duplicate = true;
						first_duplicate_idx = static_cast<int>(kmer_start_idx);
						first_duplicate_match_idx = static_cast<int>(local_match_idx);
					}
				} else if (evaluation_mode == "per_sample_threshold") {
					if (!found_first_duplicate) {
						found_first_duplicate = true;
						first_duplicate_idx = static_cast<int>(kmer_start_idx);
						first_duplicate_match_idx = static_cast<int>(local_match_idx);
					}
					++n_duplicate_kmers;
					if (static_cast<int>(kmer_start_idx) >= last_duplicate_end_idx) {
						n_duplicate_bases += k;
					} else {
						n_duplicate_bases += k - (last_duplicate_end_idx - static_cast<int>(kmer_start_idx));
					}
					last_duplicate_end_idx = static_cast<int>(kmer_start_idx) + k;
					if ((static_cast<double>(n_duplicate_bases) / static_cast<double>(kmer_start_idx + k)) <= dedup_parameter) {
						below_allowed_duplication_threshold = true;
						exceed_duplication_threshold_idx = -1;
					} else if (below_allowed_duplication_threshold) {
						below_allowed_duplication_threshold = false;
						exceed_duplication_threshold_idx = static_cast<int>(kmer_start_idx);
					}
					if (n_duplicate_bases > n_allowed_full_length_duplicate_bases) {
						done_evaluating = true;
					}
				}
			}

			if (done_evaluating) {
				break;
			}
		}
	}

	int truncating_duplicate_idx = evaluation_mode == "per_sample_threshold" ? exceed_duplication_threshold_idx : first_duplicate_idx;
	bool invalid_sample = truncating_duplicate_idx > -1 && truncating_duplicate_idx < min_sample_len;
	std::optional<std::pair<int, int>> skipped_region;

	if (invalid_sample) {
		result.sample_end_coord = -1;
		result.sample_seen_kmers = std::nullopt;
		bool is_local_duplicate = first_duplicate_match_idx > -1;
		if (is_local_duplicate) {
			result.next_start_offset = first_duplicate_match_idx + 1;
			skipped_region = std::make_pair(0, result.next_start_offset + k - 1);
		} else {
			result.next_start_offset = first_duplicate_idx + 1;
			result.duplicate_start_idx = first_duplicate_idx;
			if (first_duplicate_idx > 0) {
				skipped_region = std::make_pair(0, first_duplicate_idx + k - 1);
			}
		}
	} else if (truncating_duplicate_idx > -1) {
		result.sample_end_coord = truncating_duplicate_idx - 1 + k;
		result.duplicate_start_idx = truncating_duplicate_idx;
		result.next_start_offset = std::max(result.sample_end_coord - overlap, truncating_duplicate_idx + 1);
	}

	if (skipped_region.has_value()) {
		result.ignored_regions = {skipped_region.value()};
	} else {
		for (const auto& ir : ignored_regions_before_ambiguous_bases) {
			if (ir.first >= result.next_start_offset) {
				break;
			}
			result.ignored_regions.push_back({ir.first, std::min(ir.second, result.next_start_offset)});
		}
	}

	return result;
}

static std::vector<Region> pairs_to_regions(const std::vector<std::pair<int, int>>& pairs) {
	std::vector<Region> out;
	for (const auto& p : pairs) {
		out.push_back({static_cast<size_t>(p.first), static_cast<size_t>(p.second)});
	}
	return out;
}

using DedupSeqReturn = std::tuple<std::vector<Region>, std::vector<Region>, std::vector<Region>, std::vector<Region>, std::unordered_set<uint64_t>>;

static DedupSeqReturn deduplicate_seq(const std::string& seq, std::unordered_set<uint64_t> seen_kmers, bool retain_info, const Args& args) {
	int k = args.kmer;
	int sample_len = args.sample_len;
	int min_sample_len = args.min_sample_len;
	int overlap = args.overlap.value_or(k - 1);
	double dedup_param = args.dedup_param;
	std::string evaluation_method = args.evaluation_method;

	std::vector<Region> sample_regions;
	std::vector<size_t> masked_starts;
	std::vector<Region> ignored_regions;

	if (seen_kmers.empty() && args.seen_kmers.has_value()) {
		// no-op; the caller handles loading actual checkpoint sets
	}

	if (min_sample_len != -1 && static_cast<int>(seq.size()) < min_sample_len) {
		std::cout << "Warning: the sequence length is less than min_sample_len. Skipping this sequence.\n";
		return DedupSeqReturn{sample_regions, std::vector<Region>{}, std::vector<Region>{{Region{0, seq.size()}}}, std::vector<Region>{}, seen_kmers};
	}
	if (min_sample_len == -1 && static_cast<int>(seq.size()) < sample_len) {
		std::cout << "Warning: the sequence length is less than sample_len. Skipping this sequence.\n";
		return DedupSeqReturn{sample_regions, std::vector<Region>{}, std::vector<Region>{{Region{0, seq.size()}}}, std::vector<Region>{}, seen_kmers};
	}

	auto contig_info = get_contigs_from_chromosome(seq, args.allowed_consecutive_ambiguous_chars);
	const std::vector<Contig>& valid_contigs = contig_info.first;
	const std::vector<std::pair<size_t, size_t>>& ambiguous_regions = contig_info.second;
	bool allow_whole_contigs_override = args.allow_whole_contigs && static_cast<int>(seq.size()) < min_sample_len && valid_contigs.size() == 1;

	for (const auto& contig : valid_contigs) {
		SeqInfo seq_info(sample_len, k);
		SeqInfo* seq_info_ptr = retain_info ? &seq_info : nullptr;
		size_t sample_start = contig.start;
		size_t max_start_idx = contig.end >= static_cast<size_t>(min_sample_len) ? contig.end - static_cast<size_t>(min_sample_len) : 0;
		size_t sample_offset = 0;
		std::vector<std::pair<size_t, size_t>> contig_N_regions = contig.internal_N_regions;

		if (allow_whole_contigs_override && max_start_idx < sample_start) {
			max_start_idx = sample_start;
		}

		while (sample_start <= max_start_idx) {
			size_t sample_end = std::min(contig.end, sample_start + static_cast<size_t>(sample_len));

			size_t first_overlapping_region_idx = 0;
			bool found_first_overlapping_region = false;
			std::vector<std::pair<size_t, size_t>> sample_N_regions;
			for (size_t i = 0; i < contig_N_regions.size(); ++i) {
				auto [N_start, N_end] = contig_N_regions[i];
				if (N_start >= sample_end) {
					break;
				}
				if (N_end > sample_start + sample_offset) {
					if (!found_first_overlapping_region) {
						first_overlapping_region_idx = i;
						found_first_overlapping_region = true;
					}
					sample_N_regions.push_back({std::max(N_start, sample_start + sample_offset) - sample_start, std::min(N_end, sample_end) - sample_start});
				}
			}
			contig_N_regions = std::vector<std::pair<size_t, size_t>>(contig_N_regions.begin() + static_cast<std::ptrdiff_t>(first_overlapping_region_idx), contig_N_regions.end());

			if (seq_info_ptr != nullptr) {
				seq_info_ptr->current_global_idx = sample_start;
				seq_info_ptr->current_sample_length = sample_end - sample_start;
			}

			CheckSampleResult checked = check_sample(
				seq.substr(sample_start, sample_end - sample_start),
				sample_offset,
				sample_N_regions,
				seen_kmers,
				k,
				dedup_param,
				min_sample_len,
				std::min(overlap, static_cast<int>((sample_end - sample_start) - 1)),
				evaluation_method,
				seq_info_ptr
			);

			if (checked.sample_end_coord > -1) {
				sample_regions.push_back({sample_start, sample_start + static_cast<size_t>(checked.sample_end_coord)});
				if (seq_info_ptr == nullptr) {
					for (const auto& kv : *checked.sample_seen_kmers) {
						seen_kmers.insert(kv.first);
					}
				} else {
					for (const auto& kv : seq_info_ptr->sample_seen_kmers) {
						seen_kmers.insert(kv.first);
					}
				}
			}

			if (checked.duplicate_start_idx > -1) {
				masked_starts.push_back(sample_start + static_cast<size_t>(checked.duplicate_start_idx));
			}
			for (const auto& ir : checked.ignored_regions) {
				ignored_regions.push_back({sample_start + static_cast<size_t>(ir.first), sample_start + static_cast<size_t>(ir.second)});
			}

			sample_start = sample_start + static_cast<size_t>(checked.next_start_offset);
			if (checked.sample_end_coord == -1) {
				sample_offset = 0;
			} else {
				long long next_offset_calc = static_cast<long long>(sample_end) - static_cast<long long>(sample_start) - (k - 1);
				sample_offset = next_offset_calc > 0 ? static_cast<size_t>(next_offset_calc) : 0;
			}
			if (seq_info_ptr != nullptr) {
				seq_info_ptr->move_internal_pointer(static_cast<size_t>(checked.next_start_offset), checked.sample_end_coord > -1);
			}
		}

		if (sample_start < contig.end) {
			ignored_regions.push_back({sample_start, contig.end});
		}
	}

	std::vector<Region> masked_regions = condense_masked_regions(masked_starts);
	std::vector<Region> ambiguous_region_vec = regions_from_ranges(ambiguous_regions);
	return DedupSeqReturn{sample_regions, masked_regions, ignored_regions, ambiguous_region_vec, seen_kmers};
}

static std::pair<std::unordered_set<uint64_t>, GenomeAnnotations> deduplicate_genome(const std::string& fasta, std::unordered_set<uint64_t> seen_kmers, const Args& args) {
	GenomeAnnotations kmer_region_annotations;
	for (const auto& rec : read_fasta_records(fasta)) {
		const std::string& seqname = rec.first;
		std::cout << "	Processing contig " << seqname << std::endl;
		std::string clean_sequence = get_clean_sequence(rec.second);
		auto result = deduplicate_seq(clean_sequence, seen_kmers, args.retain_info, args);
		seen_kmers = std::move(std::get<4>(result));
		AnnotationMap ann;
		ann["samples"] = std::move(std::get<0>(result));
		ann["masks"] = std::move(std::get<1>(result));
		ann["ignored"] = std::move(std::get<2>(result));
		ann["ambiguous"] = std::move(std::get<3>(result));
		kmer_region_annotations.push_back({seqname, std::move(ann)});
	}
	return {seen_kmers, kmer_region_annotations};
}

static void write_region_beds(const GenomeAnnotations& kmer_region_annotations, const std::string& file_prefix, bool write_ambiguous = false, bool write_ignored = false, bool write_masks = false) {
	std::vector<std::string> write_bed_types = {"samples"};
	if (write_ambiguous) write_bed_types.push_back("ambiguous");
	if (write_ignored) write_bed_types.push_back("ignored");
	if (write_masks) write_bed_types.push_back("masks");
	for (const auto& bed_type : write_bed_types) {
		std::ofstream f(file_prefix + "." + bed_type + ".bed");
		for (const auto& entry : kmer_region_annotations) {
			const std::string& seqname = entry.first;
			auto it = entry.second.find(bed_type);
			if (it == entry.second.end()) {
				continue;
			}
			for (const auto& region : it->second) {
				f << seqname << '\t' << region.start << '\t' << region.end << '\n';
			}
		}
	}
}

static void json_write_string(std::ostream& os, const std::string& s) {
	os << '"';
	for (char c : s) {
		switch (c) {
			case '"': os << "\\\""; break;
			case '\\': os << "\\\\"; break;
			case '\n': os << "\\n"; break;
			case '\r': os << "\\r"; break;
			case '\t': os << "\\t"; break;
			default: os << c; break;
		}
	}
	os << '"';
}

static std::string args_to_json(const Args& args) {
	std::ostringstream oss;
	oss << "{\n";
	oss << "    \"input\": [";
	for (size_t i = 0; i < args.input.size(); ++i) {
		if (i) oss << ", ";
		json_write_string(oss, args.input[i]);
	}
	oss << "],\n";
	oss << "    \"allow_whole_contigs\": " << (args.allow_whole_contigs ? "true" : "false") << ",\n";
	oss << "    \"dedup_param\": " << args.dedup_param << ",\n";
	oss << "    \"evaluation_method\": "; json_write_string(oss, args.evaluation_method); oss << ",\n";
	oss << "    \"kmer\": " << args.kmer << ",\n";
	oss << "    \"sample_len\": " << args.sample_len << ",\n";
	oss << "    \"min_sample_len\": " << args.min_sample_len << ",\n";
	oss << "    \"allowed_consecutive_ambiguous_chars\": " << args.allowed_consecutive_ambiguous_chars << ",\n";
	oss << "    \"output_dir\": "; json_write_string(oss, args.output_dir); oss << ",\n";
	oss << "    \"seen_kmers\": ";
	if (args.seen_kmers.has_value()) json_write_string(oss, args.seen_kmers.value()); else oss << "null";
	oss << ",\n";
	oss << "    \"retain_info\": " << (args.retain_info ? "true" : "false") << ",\n";
	oss << "    \"save_every\": " << args.save_every << ",\n";
	oss << "    \"overlap\": ";
	if (args.overlap.has_value()) oss << args.overlap.value(); else oss << "null";
	oss << ",\n";
	oss << "    \"save_kmers_at_end\": " << (args.save_kmers_at_end ? "true" : "false") << ",\n";
	oss << "    \"write_ambiguous_beds\": " << (args.write_ambiguous_beds ? "true" : "false") << ",\n";
	oss << "    \"write_ignored_beds\": " << (args.write_ignored_beds ? "true" : "false") << ",\n";
	oss << "    \"write_masked_beds\": " << (args.write_masked_beds ? "true" : "false") << ",\n";
	oss << "    \"seed\": " << args.seed << "\n";
	oss << "}";
	return oss.str();
}

static void deduplicate(const Args& args) {
	if (!fs::is_directory(args.output_dir)) {
		fs::create_directories(args.output_dir);
	} else {
		std::cout << "Output directory already exists. Continue and potentially overwrite files? (y/n): ";
		std::string delete_check;
		std::getline(std::cin, delete_check);
		std::string lower = delete_check;
		std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
		if (lower != "y" && lower != "yes") {
			std::cout << "Exiting...\n";
			std::exit(1);
		}
	}

	{
		std::ofstream f(fs::path(args.output_dir) / "config.json");
		f << args_to_json(args) << std::endl;
	}

	std::vector<std::string> fastas;
	if (args.input.size() > 1 || ends_with(args.input[0], ".fa") || ends_with(args.input[0], ".fasta") || ends_with(args.input[0], ".fasta.gz") || ends_with(args.input[0], ".fna") || ends_with(args.input[0], ".fna.gz")) {
		fastas = args.input;
	} else if (args.input.size() == 1 && (ends_with(args.input[0], ".txt") || ends_with(args.input[0], ".list"))) {
		for (const auto& line : read_text_lines(args.input[0])) {
			fastas.push_back(line);
		}
	}

	std::vector<std::string> valid_fastas;
	for (const auto& fasta : fastas) {
		if (fs::is_regular_file(fasta)) {
			valid_fastas.push_back(fasta);
		}
	}
	std::vector<std::string> invalid_fastas;
	for (const auto& fasta : fastas) {
		if (!fs::is_regular_file(fasta)) {
			invalid_fastas.push_back(fasta);
		}
	}
	if (!invalid_fastas.empty()) {
		std::cout << "Warning: could not find the following fastas: ";
		for (size_t i = 0; i < invalid_fastas.size(); ++i) {
			if (i) std::cout << ", ";
			std::cout << invalid_fastas[i];
		}
		std::cout << "\n";
	}

	{
		std::ofstream f(fs::path(args.output_dir) / "basename_fasta_match.txt");
		for (const auto& fasta : valid_fastas) {
			f << get_fasta_basename(fasta) << '\t' << fasta << '\n';
		}
	}

	std::vector<size_t> saved_kmers_idxs;
	std::vector<std::string> needed_beds = {"samples"};
	if (args.write_ambiguous_beds) needed_beds.push_back("ambiguous");
	if (args.write_ignored_beds) needed_beds.push_back("ignored");
	if (args.write_masked_beds) needed_beds.push_back("masks");
	int last_beds_idx = -1;
	for (size_t i = 0; i < valid_fastas.size(); ++i) {
		std::string fasta_basename = get_fasta_basename(valid_fastas[i]);
		std::string saved_kmers_file = (fs::path(args.output_dir) / (fasta_basename + ".kmers.bin")).string();
		if (fs::exists(saved_kmers_file)) {
			saved_kmers_idxs.push_back(i);
		}
		bool all_beds_exist = true;
		for (const auto& bed_type : needed_beds) {
			std::string bed_file = (fs::path(args.output_dir) / (fasta_basename + "." + bed_type + ".bed")).string();
			if (!fs::exists(bed_file)) {
				all_beds_exist = false;
				break;
			}
		}
		if (all_beds_exist) {
			last_beds_idx = static_cast<int>(i);
		}
	}

	int last_checkpoint_idx = -1;
	for (size_t saved_kmers_idx : saved_kmers_idxs) {
		if (static_cast<int>(saved_kmers_idx) <= last_beds_idx) {
			last_checkpoint_idx = static_cast<int>(saved_kmers_idx);
		}
	}

	std::optional<std::string> last_checkpoint_kmers_file;
	if (last_checkpoint_idx > -1) {
		std::string last_checkpoint_fasta = valid_fastas[static_cast<size_t>(last_checkpoint_idx)];
		std::string last_checkpoint_basename = get_fasta_basename(last_checkpoint_fasta);
		last_checkpoint_kmers_file = (fs::path(args.output_dir) / (last_checkpoint_basename + ".kmers.bin")).string();
	}

	size_t next_fasta_idx = static_cast<size_t>(last_checkpoint_idx + 1);
	std::vector<std::string> process_fastas(valid_fastas.begin() + static_cast<std::ptrdiff_t>(next_fasta_idx), valid_fastas.end());

	std::unordered_set<uint64_t> seen_kmers;
	if (!args.seen_kmers.has_value()) {
		if (!last_checkpoint_kmers_file.has_value()) {
			seen_kmers = {};
		} else {
			seen_kmers = read_seen_kmers(last_checkpoint_kmers_file.value());
		}
	} else {
		if (!last_checkpoint_kmers_file.has_value()) {
			seen_kmers = read_seen_kmers(args.seen_kmers.value());
		} else if (fs::weakly_canonical(args.seen_kmers.value()) == fs::weakly_canonical(last_checkpoint_kmers_file.value())) {
			seen_kmers = read_seen_kmers(args.seen_kmers.value());
		} else {
			std::cout << "Warning: a seen kmers file was supplied, but it does not match the latest checkpoint kmer set. Taking the union of the two sets of kmers and proceeding with deduplication.\n";
			auto seen_kmers_1 = read_seen_kmers(args.seen_kmers.value());
			auto seen_kmers_2 = read_seen_kmers(last_checkpoint_kmers_file.value());
			seen_kmers = seen_kmers_1;
			seen_kmers.insert(seen_kmers_2.begin(), seen_kmers_2.end());
		}
	}

	for (size_t i = 0; i < process_fastas.size(); ++i) {
		const std::string& fasta = process_fastas[i];
		std::string fasta_basename = get_fasta_basename(fasta);
		bool all_beds_exist = true;
		for (const auto& bed_type : needed_beds) {
			std::string bed_file = (fs::path(args.output_dir) / (fasta_basename + "." + bed_type + ".bed")).string();
			if (!fs::exists(bed_file)) {
				all_beds_exist = false;
				break;
			}
		}
		std::string out_prefix = (fs::path(args.output_dir) / fasta_basename).string();

		if (all_beds_exist) {
			std::cout << "Found existing output for " << fasta << ", skipping deduplication and computing seen kmers from existing samples\n";
			std::string samples_bed = (fs::path(args.output_dir) / (fasta_basename + ".samples.bed")).string();
			auto fasta_seen_kmers = compute_seen_kmers_from_samples_bed_and_fasta(samples_bed, fasta, args.kmer);
			seen_kmers.insert(fasta_seen_kmers.begin(), fasta_seen_kmers.end());
		} else {
			std::cout << "Deduplicating " << fasta << std::endl;
			auto result = deduplicate_genome(fasta, seen_kmers, args);
			seen_kmers = std::move(result.first);
			write_region_beds(result.second, out_prefix, args.write_ambiguous_beds, args.write_ignored_beds, args.write_masked_beds);
		}

		bool save_kmers_to_file = args.save_every > 0 && ((next_fasta_idx + i + 1) % static_cast<size_t>(args.save_every) == 0);
		if (save_kmers_to_file) {
			write_seen_kmers(seen_kmers, out_prefix);
		}
	}

	if (args.save_kmers_at_end) {
		std::string final_seen_kmer_file_basename = "final";
		int idx_suffix = 1;
		while (fs::exists(fs::path(args.output_dir) / (final_seen_kmer_file_basename + ".kmers.bin"))) {
			final_seen_kmer_file_basename = "final_" + std::to_string(idx_suffix);
			++idx_suffix;
		}
		write_seen_kmers(seen_kmers, (fs::path(args.output_dir) / final_seen_kmer_file_basename).string());
	}
}

static void print_usage(const char* program_name) {
	std::cout << "Usage: " << program_name << " [options] <input_files...>\n\n"
	          << "Options:\n"
	          << "  --allow_whole_contigs\n"
	          << "  -d, --dedup_param <float>\n"
	          << "  -e, --evaluation_method <str>\n"
	          << "  -k, --kmer <int>\n"
	          << "  -l, --sample_len <int>\n"
	          << "  -m, --min_sample_len <int>\n"
	          << "  -n, --allowed_consecutive_ambiguous_chars <int>\n"
	          << "  -o, --output_dir <path>\n"
	          << "  -p, --seen_kmers <path>\n"
	          << "  -r, --retain_info\n"
	          << "  -s, --save_every <int>\n"
	          << "  -v, --overlap <int>\n"
	          << "  --save_kmers_at_end\n"
	          << "  --write_ambiguous_beds\n"
	          << "  --write_ignored_beds\n"
	          << "  --write_masked_beds\n"
	          << "  -seed, --seed <int>\n";
}

static Args parse_args(int argc, char* argv[]) {
	Args args;
	for (int i = 1; i < argc; ++i) {
		std::string arg = argv[i];
		if (arg == "-h" || arg == "--help") {
			print_usage(argv[0]);
			std::exit(0);
		} else if (arg == "--allow_whole_contigs") {
			args.allow_whole_contigs = true;
		} else if (arg == "-d" || arg == "--dedup_param") {
			if (i + 1 < argc) args.dedup_param = std::stod(argv[++i]);
		} else if (arg == "-e" || arg == "--evaluation_method") {
			if (i + 1 < argc) args.evaluation_method = argv[++i];
		} else if (arg == "-k" || arg == "--kmer") {
			if (i + 1 < argc) args.kmer = std::stoi(argv[++i]);
		} else if (arg == "-l" || arg == "--sample_len") {
			if (i + 1 < argc) args.sample_len = std::stoi(argv[++i]);
		} else if (arg == "-m" || arg == "--min_sample_len") {
			if (i + 1 < argc) args.min_sample_len = std::stoi(argv[++i]);
		} else if (arg == "-n" || arg == "--allowed_consecutive_ambiguous_chars") {
			if (i + 1 < argc) args.allowed_consecutive_ambiguous_chars = std::stoi(argv[++i]);
		} else if (arg == "-o" || arg == "--output_dir") {
			if (i + 1 < argc) args.output_dir = argv[++i];
		} else if (arg == "-p" || arg == "--seen_kmers") {
			if (i + 1 < argc) args.seen_kmers = std::string(argv[++i]);
		} else if (arg == "-r" || arg == "--retain_info") {
			args.retain_info = true;
		} else if (arg == "-s" || arg == "--save_every") {
			if (i + 1 < argc) args.save_every = std::stoi(argv[++i]);
		} else if (arg == "-v" || arg == "--overlap") {
			if (i + 1 < argc) args.overlap = std::stoi(argv[++i]);
		} else if (arg == "--save_kmers_at_end") {
			args.save_kmers_at_end = true;
		} else if (arg == "--write_ambiguous_beds") {
			args.write_ambiguous_beds = true;
		} else if (arg == "--write_ignored_beds") {
			args.write_ignored_beds = true;
		} else if (arg == "--write_masked_beds") {
			args.write_masked_beds = true;
		} else if (arg == "-seed" || arg == "--seed") {
			if (i + 1 < argc) args.seed = std::stoi(argv[++i]);
		} else if (!arg.empty() && arg[0] != '-') {
			args.input.push_back(arg);
		}
	}
	if (args.min_sample_len == -1) {
		args.min_sample_len = args.sample_len;
	}
	if (!args.overlap.has_value()) {
		args.overlap = args.kmer - 1;
	}
	return args;
}

int main(int argc, char* argv[]) {
	try {
		if (argc < 2) {
			print_usage(argv[0]);
			return 1;
		}

		Args args = parse_args(argc, argv);
		rng_gen.seed(args.seed);

		for (const auto& f : args.input) {
			if (!fs::is_regular_file(f)) {
				throw std::runtime_error("Error: could not find supplied list of fasta files");
			}
			type_check(f);
		}

		if (args.kmer < 1) {
			throw std::runtime_error("Error: kmer size must be a positive integer");
		}
		if (args.kmer < 16) {
			std::cout << "Warning: small kmer sizes will result in very strict deduplication. Consider increasing the kmer size in the range 16-32 (default: 32)\n";
		}
		if (args.kmer > 32) {
			std::cout << "Warning: kmer sizes above 32 may not work or may result in very slow evaluation. Consider decreasing the kmer size in the range 16-32 (default: 32)\n";
		}
		if (args.evaluation_method != "per_kmer" && args.evaluation_method != "per_sample_agnostic" && args.evaluation_method != "per_sample_threshold") {
			throw std::runtime_error("Error: evaluation method must be one of per_kmer, per_sample_agnostic, or per_sample_threshold");
		}
		if (args.dedup_param < 0.0 || args.dedup_param > 1.0) {
			throw std::runtime_error("Error: deduplication parameter must be between 0.0 and 1.0");
		}
		if (args.min_sample_len > args.sample_len) {
			std::cout << "Warning: min sample length cannot be bigger than the default sample length; defaulting to equal the standard sample length\n";
			args.min_sample_len = args.sample_len;
		}
		if (args.sample_len < args.kmer) {
			std::cout << "Warning: sample len is smaller than k, meaning no deduplication will occur\n";
		}

		deduplicate(args);
		return 0;
	} catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		return 1;
	}
}