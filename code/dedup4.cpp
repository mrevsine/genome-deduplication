/**
 * Genome Deduplication Tool - C++ Implementation (Version 3)
 * 
 * Deduplicates and samples genomic sequences simultaneously
 * Only considers a kmer if it is within a sample
 * 
 * Basic algorithm:
 *     1. Sample sequence to first `sample_len` bases
 *     2. Iterate over all kmers in sample
 *     3. If kmer not in seen_kmers, add it to seen_kmers
 *     4. If kmer in seen_kmers, jump to next possible sample and repeat
 * 
 * The program divides a fasta into three categories:
 *     1. Samples - valid regions with no repeats, on which we train
 *     2. Masked regions - regions containing duplicates
 *     3. Ignored regions - regions that we threw out because we found a
 *        duplicate later on in the sample region
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <cstdint>
#include <random>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <cstring>
#include <cctype>
#include <zlib.h>
#include <chrono>
#include <ctime>

namespace fs = std::filesystem;

// Structure to hold genomic regions (start inclusive, end exclusive)
struct Region {
    size_t start;
    size_t end;
    
    Region(size_t s, size_t e) : start(s), end(e) {}
};

// Structure to hold deduplication results for a sequence
struct DeduplicationResult {
    std::vector<Region> sample_regions;
    std::vector<Region> masked_regions;
    std::vector<Region> skipped_regions;
    std::vector<Region> ambiguous_regions;
};

// Structure to hold program arguments
struct Args {
    std::vector<std::string> input;
    int kmer = 32;
    int sample_len = 1000;
    int min_sample_len = -1; // Will be set to sample_len if not specified
    std::string output_dir = "dedup_out";
    std::string seen_kmers_file;
    // Per-kmer mode retain percent (0-1)
    double per_kmer_retain_pct = 0.0;
    // Per-sample agnostic retain percent (0-1)
    double agnostic_retain_pct = 0.0;
    // Per-sample threshold settings (interpreted same as Python workspace)
    double ambiguous_base_threshold = 0.0;
    double duplicate_base_threshold = 0.0;
    // Evaluation method: per_kmer | per_sample_agnostic | per_sample_threshold
    std::string evaluation_method = "per_kmer";
    int save_every = 0;
    bool no_overlap = false;
    bool no_save_kmers_at_end = false;
    int seed = 123;
    bool yes = false;
    // Config output controls
    bool print_config = false;
    bool write_config = true;
};

// Global random number generator
std::mt19937 rng_gen;

/**
 * Encode a k-mer string into a 64-bit integer
 * A=0, C=1, G=2, T=3
 */
uint64_t encode_kmer(const std::string& kmer) {
    uint64_t kmer_num = 0;
    for (char c : kmer) {
        kmer_num = (kmer_num << 2);
        switch(c) {
            case 'A': kmer_num |= 0; break;
            case 'C': kmer_num |= 1; break;
            case 'G': kmer_num |= 2; break;
            case 'T': kmer_num |= 3; break;
            default: 
                throw std::runtime_error("Invalid nucleotide in k-mer: " + std::string(1, c));
        }
    }
    return kmer_num;
}

/**
 * Decode a k-mer integer back to string
 */
std::string decode_kmer(uint64_t kmer_num, int k = 32) {
    std::string kmer(k, 'A');
    for (int i = k - 1; i >= 0; --i) {
        int nucleotide_code = kmer_num & 3;
        switch(nucleotide_code) {
            case 0: kmer[i] = 'A'; break;
            case 1: kmer[i] = 'C'; break;
            case 2: kmer[i] = 'G'; break;
            case 3: kmer[i] = 'T'; break;
        }
        kmer_num >>= 2;
    }
    return kmer;
}

/**
 * Get the basename of a fasta file (without path and extensions)
 */
std::string get_fasta_basename(const std::string& fasta) {
    fs::path p(fasta);
    std::string filename = p.filename().string();
    
    // Remove .gz if present
    if (filename.size() > 3 && filename.substr(filename.size() - 3) == ".gz") {
        filename = filename.substr(0, filename.size() - 3);
    }
    
    // Remove fasta extensions
    size_t last_dot = filename.find_last_of('.');
    if (last_dot != std::string::npos) {
        return filename.substr(0, last_dot);
    }
    return filename;
}

/**
 * Condense a list of indices into contiguous regions
 * Updated version: does not add k to region end (workspace version change)
 * e.g. [2,3,4,7,8,20] -> [(2,5), (7,9), (20,21)]
 */
std::vector<Region> condense_masked_regions(const std::vector<size_t>& masked) {
    std::vector<Region> masked_regions;
    if (masked.empty()) {
        return masked_regions;
    }
    
    size_t region_start = masked[0];
    for (size_t i = 0; i < masked.size() - 1; ++i) {
        // If the next masked position is not contiguous with the current one
        if (masked[i] + 1 != masked[i + 1]) {
            // Workspace version: don't add k, just add 1
            masked_regions.emplace_back(region_start, masked[i] + 1);
            region_start = masked[i + 1];
        }
    }
    // Workspace version: don't add k, just add 1
    masked_regions.emplace_back(region_start, masked.back() + 1);
    return masked_regions;
}

/**
 * Duplicate base count across seq_info considering k
 * Marks k bases for each index with duplicate codes (4,5)
 */
static size_t get_duplicate_base_count(const std::vector<int>& seq_info, int k, const std::vector<int>& duplicate_codes = {4,5}) {
    std::vector<uint8_t> dup(seq_info.size(), 0);
    for (size_t i = 0; i < seq_info.size(); ++i) {
        int n = seq_info[i];
        bool is_dup = false;
        for (int code : duplicate_codes) {
            if (n == code) { is_dup = true; break; }
        }
        if (is_dup) {
            size_t end = std::min(seq_info.size(), i + static_cast<size_t>(k));
            for (size_t j = i; j < end; ++j) dup[j] = 1;
        }
    }
    size_t sum = 0;
    for (uint8_t v : dup) sum += v;
    return sum;
}

struct CheckOverallResult {
    int checked_sample_len;
    int duplicate_start_idx;
    int ambiguous_idx;
    std::pair<int,int> ignored_region; // -1,-1 if none
    size_t next_start_offset;
    std::unordered_map<uint64_t, size_t> sample_seen_kmers;
    bool valid_sample_kmers;
    std::vector<int> seq_info_out; // updated per-sample seq_info
};

/**
 * Overall sample evaluation (per_sample_agnostic / per_sample_threshold)
 * Faithful translation of dedup2_workspace.py::check_sample_overall
 */
static CheckOverallResult check_sample_overall(
    const std::string& seq,
    const std::vector<int>& seq_info_in,
    size_t seq_info_offset,
    const std::unordered_set<uint64_t>& seen_kmers,
    int k,
    int min_sample_len,
    bool no_overlap_samples,
    const Args& args
) {
    // Region type codes
    const int RT_UNANNOTATED = 0;
    const int RT_UNIQUE = 1;
    const int RT_IGNORED = 2;
    const int RT_AMBIGUOUS = 3;
    const int RT_INTERNAL = 4;
    const int RT_GLOBAL = 5;

    CheckOverallResult R;
    R.checked_sample_len = -1;
    R.duplicate_start_idx = -1;
    R.ambiguous_idx = -1;
    R.ignored_region = {-1,-1};
    R.valid_sample_kmers = false;
    R.seq_info_out.assign(seq_info_in.begin(), seq_info_in.end());

    const int final_start_idx = static_cast<int>(seq.size()) - k; // last kmer start
    int sample_end_coord = static_cast<int>(seq.size());
    R.next_start_offset = no_overlap_samples ? seq.size() : (seq.size() - k + 1);

    // Build first draft of sample_seen_kmers from prior uniques up to min(offset, final_start_idx+1)
    size_t unique_limit = std::min(seq_info_offset, static_cast<size_t>(final_start_idx + 1));
    for (size_t i = 0; i < unique_limit; ++i) {
        if (R.seq_info_out[i] == RT_UNIQUE) {
            std::string kmer = seq.substr(i, k);
            uint64_t num = encode_kmer(kmer);
            if (R.sample_seen_kmers.find(num) == R.sample_seen_kmers.end()) {
                R.sample_seen_kmers[num] = static_cast<size_t>(i);
            }
        }
    }

    // Re-check prior internal repeats before offset
    for (size_t i = 0; i < seq_info_offset && i + k <= seq.size(); ++i) {
        if (R.seq_info_out[i] == RT_INTERNAL) {
            std::string kmer = seq.substr(i, k);
            uint64_t num = encode_kmer(kmer);
            if (R.sample_seen_kmers.find(num) == R.sample_seen_kmers.end()) {
                R.sample_seen_kmers[num] = i;
                R.seq_info_out[i] = RT_UNIQUE;
            } else if (seen_kmers.find(num) != seen_kmers.end()) {
                R.seq_info_out[i] = RT_GLOBAL;
            }
        }
    }

    // Mark new ambiguous bases (Ns) after offset
    for (size_t idx = seq_info_offset; idx < seq.size(); ++idx) {
        if (seq[idx] == 'N') {
            R.seq_info_out[idx] = RT_AMBIGUOUS;
        }
    }

    // Determine contig boundaries (non-ambiguous runs)
    std::vector<std::pair<int,int>> contigs;
    int contig_start = static_cast<int>(seq_info_offset);
    for (size_t idx = seq_info_offset; idx < seq.size(); ++idx) {
        if (R.seq_info_out[idx] == RT_AMBIGUOUS) {
            if (static_cast<int>(idx) > contig_start) {
                contigs.emplace_back(contig_start, static_cast<int>(idx));
            }
            contig_start = static_cast<int>(idx) + 1;
        }
    }
    if (contig_start < static_cast<int>(seq.size())) {
        contigs.emplace_back(contig_start, static_cast<int>(seq.size()));
    }

    // Process contigs: mark last k-1 bases as ignored, iterate kmers
    for (const auto& ce : contigs) {
        int cstart = ce.first;
        int cend = ce.second; // exclusive
        int ignored_start = std::max(cstart, cend - k + 1);
        for (int i = ignored_start; i < cend; ++i) {
            R.seq_info_out[i] = RT_IGNORED;
        }
        for (int kstart = cstart; kstart <= cend - k; ++kstart) {
            // kmer should not contain N by construction
            std::string kmer = seq.substr(kstart, k);
            uint64_t num = encode_kmer(kmer);
            if (R.sample_seen_kmers.find(num) == R.sample_seen_kmers.end() &&
                seen_kmers.find(num) == seen_kmers.end()) {
                R.sample_seen_kmers[num] = static_cast<size_t>(kstart);
                R.seq_info_out[kstart] = RT_UNIQUE;
            } else {
                if (R.sample_seen_kmers.find(num) != R.sample_seen_kmers.end()) {
                    R.seq_info_out[kstart] = RT_INTERNAL;
                } else {
                    R.seq_info_out[kstart] = RT_GLOBAL;
                }
            }
        }
    }

    // Compute sequence-level metrics
    size_t n_ambiguous_bases = 0;
    for (int v : R.seq_info_out) if (v == RT_AMBIGUOUS) ++n_ambiguous_bases;
    size_t n_duplicate_bases = get_duplicate_base_count(R.seq_info_out, k);

    // Decide whether to keep sample
    bool valid_sequence = true;
    if (args.evaluation_method == "per_sample_agnostic") {
        double p = args.agnostic_retain_pct;
        if (n_ambiguous_bases > 0) {
            valid_sequence = false; // TODO - allow for a more complex Ns policy, like some allowed # of consecutive Ns a la Evo-2
        } else if (n_duplicate_bases > 0) {
            std::uniform_real_distribution<double> dist(0.0, 1.0);
            if (dist(rng_gen) >= p) {
                valid_sequence = false;
            }
        }
    } else if (args.evaluation_method == "per_sample_threshold") {
        // Faithful to Python: thresholds compared as counts
        double allowed_amb = args.ambiguous_base_threshold * seq.size();
        double allowed_dup = args.duplicate_base_threshold * seq.size();
        if (static_cast<double>(n_ambiguous_bases) > allowed_amb || static_cast<double>(n_duplicate_bases) > allowed_dup) {
            valid_sequence = false;
        }
    }

    if (!valid_sequence) {
        // Find first problematic base
        std::pair<int,int> first_issue = {-1,-1};
        for (size_t i = 0; i < R.seq_info_out.size(); ++i) {
            int n = R.seq_info_out[i];
            if (n == RT_AMBIGUOUS || n == RT_INTERNAL || n == RT_GLOBAL) {
                first_issue = {static_cast<int>(i), n};
                break;
            }
        }
        if (first_issue.second == RT_AMBIGUOUS) {
            R.ambiguous_idx = first_issue.first;
            if (R.ambiguous_idx > 0) {
                R.ignored_region = {0, R.ambiguous_idx};
            }
            R.next_start_offset = static_cast<size_t>(R.ambiguous_idx + 1);
        } else if (first_issue.second == RT_INTERNAL) {
            // Use original kmer position
            std::string kmer = seq.substr(first_issue.first, k);
            uint64_t num = encode_kmer(kmer);
            size_t orig = R.sample_seen_kmers[num];
            R.ignored_region = {0, static_cast<int>(orig + 1)};
            R.next_start_offset = orig + 1;
        } else if (first_issue.second == RT_GLOBAL) {
            R.duplicate_start_idx = first_issue.first;
            if (R.duplicate_start_idx > 0) {
                R.ignored_region = {0, R.duplicate_start_idx};
            }
            R.next_start_offset = static_cast<size_t>(R.duplicate_start_idx + 1);
        }
        R.checked_sample_len = -1;
        R.valid_sample_kmers = false;
        R.sample_seen_kmers.clear();
    } else {
        R.checked_sample_len = sample_end_coord;
        R.valid_sample_kmers = true;
    }

    return R;
}

/**
 * Check a sample for duplicates
 */
struct CheckSampleResult {
    int checked_sample_len;
    int duplicate_start_idx;
    int ambiguous_idx;
    std::pair<int, int> ignored_region; // -1, -1 if none
    size_t next_start_offset;
    std::unordered_map<uint64_t, size_t> sample_seen_kmers;
    bool valid_sample_kmers;
};

CheckSampleResult check_sample(
    const std::string& seq,
    const std::unordered_set<uint64_t>& seen_kmers,
    int k,
    double allowed_duplicate_rate,
    int min_sample_len,
    bool no_overlap_samples,
    const std::string& N_policy = "allow_none",
    bool debug = false
) {
    CheckSampleResult result;
    result.checked_sample_len = -1;
    result.duplicate_start_idx = -1;
    result.ambiguous_idx = -1;
    result.ignored_region = {-1, -1};
    result.valid_sample_kmers = false;
    
    size_t sample_end_coord = seq.length();
    result.next_start_offset = no_overlap_samples ? seq.length() : seq.length() - k + 1;
    
    // Check for Ns
    if (N_policy == "allow_none") {
        size_t N_index = seq.find('N');
        if (N_index != std::string::npos) {
            if (N_index < static_cast<size_t>(min_sample_len)) {
                result.ambiguous_idx = N_index;
                if (N_index > 0) {
                    // Workspace version: no k added here
                    result.ignored_region = {0, static_cast<int>(N_index)};
                }
                result.next_start_offset = N_index + 1;
                result.checked_sample_len = -1;
                return result;
            } else {
                sample_end_coord = N_index;
            }
        }
    }
    
    // Loop through all kmers in this possible sample
    if (sample_end_coord >= static_cast<size_t>(min_sample_len)) {
        for (size_t kmer_start_idx = 0; kmer_start_idx <= sample_end_coord - k; ++kmer_start_idx) {
            std::string kmer = seq.substr(kmer_start_idx, k);
            uint64_t kmer_num = encode_kmer(kmer);
            
            if (debug && kmer_start_idx == 0) {
                std::cout << "    First kmer: " << kmer << std::endl;
                std::cout << "    Encoded: " << kmer_num << std::endl;
                std::cout << "    In seen_kmers: " << (seen_kmers.find(kmer_num) != seen_kmers.end()) << std::endl;
            }
            
            // If we haven't seen this kmer before, record it
            if (seen_kmers.find(kmer_num) == seen_kmers.end() &&
                result.sample_seen_kmers.find(kmer_num) == result.sample_seen_kmers.end()) {
                result.sample_seen_kmers[kmer_num] = kmer_start_idx;
            } else {
                if (debug) {
                    std::cout << "    Found dup at kmer_start_idx=" << kmer_start_idx 
                              << " kmer=" << kmer.substr(0,8) << "..."
                              << " in_seen=" << (seen_kmers.find(kmer_num) != seen_kmers.end())
                              << " in_sample=" << (result.sample_seen_kmers.find(kmer_num) != result.sample_seen_kmers.end())
                              << " sample_size=" << result.sample_seen_kmers.size()
                              << std::endl;
                    
                    // Check if internal
                    auto it_internal = result.sample_seen_kmers.find(kmer_num);
                    if (it_internal != result.sample_seen_kmers.end()) {
                        std::cout << "      -> INTERNAL match, original at " << it_internal->second << std::endl;
                    }
                }
                // Found a duplicate - check if we should allow it
                if (allowed_duplicate_rate > 0) {
                    std::uniform_real_distribution<double> dist(0.0, 1.0);
                    if (dist(rng_gen) < allowed_duplicate_rate) {
                        continue;
                    }
                }
                
                // Process this duplicate
                if (kmer_start_idx < static_cast<size_t>(min_sample_len)) {
                    result.checked_sample_len = -1;
                    
                    // Check if internal duplicate
                    auto it = result.sample_seen_kmers.find(kmer_num);
                    if (it != result.sample_seen_kmers.end()) {
                        // Internal match - workspace version changes
                        result.next_start_offset = it->second + 1;
                        // Workspace version: ignored region is (0, next_start_offset)
                        result.ignored_region = {0, static_cast<int>(result.next_start_offset)};
                    } else {
                        // Global match - workspace version changes
                        result.next_start_offset = kmer_start_idx + 1;
                        result.duplicate_start_idx = kmer_start_idx;
                        if (kmer_start_idx > 0) {
                            // Workspace version: ignored region is (0, kmer_start_idx) without k
                            result.ignored_region = {0, static_cast<int>(kmer_start_idx)};
                        }
                    }
                    result.valid_sample_kmers = false;
                } else {
                    // Offending kmer past min_sample_len
                    result.checked_sample_len = kmer_start_idx;
                    result.duplicate_start_idx = kmer_start_idx;
                    result.next_start_offset = kmer_start_idx + 1;
                    result.valid_sample_kmers = true;
                }
                break;
            }
        }
        
        // If no duplicate found (valid_sample_kmers is still false from initialization if no duplicate was found)
        if (result.checked_sample_len == -1 && result.duplicate_start_idx == -1 && result.valid_sample_kmers == false && result.ignored_region.first == -1) {
            result.checked_sample_len = sample_end_coord;
            result.valid_sample_kmers = true;
        }
    }
    
    return result;
}

/**
 * Deduplicate a single sequence
 */
DeduplicationResult deduplicate_seq(
    const std::string& seq,
    std::unordered_set<uint64_t>& seen_kmers,
    const Args& args,
    bool debug = false
) {
    DeduplicationResult result;
    
    // Check validity of inputs
    if (seq.length() < static_cast<size_t>(args.min_sample_len)) {
        std::cerr << "Warning: sequence length is less than min_sample_len. Skipping this sequence." << std::endl;
        result.skipped_regions.emplace_back(0, seq.length());
        return result;
    }
    
    std::vector<size_t> masked_starts;
    std::vector<size_t> ambiguous_positions;
    
    size_t max_start_idx = seq.length() - args.min_sample_len;
    size_t sample_start = 0;
    
    // Evaluation-method specific state for per-sample modes
    std::vector<int> sample_seq_info; // empty by default
    size_t seq_info_offset = 0;

    // Investigate every possible sample
    while (sample_start <= max_start_idx) {
        size_t sample_end = std::min(seq.length(), sample_start + args.sample_len);
        std::string sample_seq = seq.substr(sample_start, sample_end - sample_start);

        if (args.evaluation_method == "per_kmer") {
            // Determine duplication allowance based on evaluation method
            double allowed_duplicate_rate = args.per_kmer_retain_pct;

            CheckSampleResult check_result = check_sample(
                sample_seq,
                seen_kmers,
                args.kmer,
                allowed_duplicate_rate,
                args.min_sample_len,
                args.no_overlap,
                "allow_none",
                debug && sample_start < 200
            );

            if (debug && sample_start < 200) {
                std::cout << "  sample_start=" << sample_start
                          << " checked_len=" << check_result.checked_sample_len
                          << " dup_idx=" << check_result.duplicate_start_idx
                          << " ignored=(" << check_result.ignored_region.first << "," << check_result.ignored_region.second << ")"
                          << " next_offset=" << check_result.next_start_offset
                          << " sample_kmers=" << check_result.sample_seen_kmers.size()
                          << std::endl;
            }

            if (check_result.checked_sample_len > -1) {
                result.sample_regions.emplace_back(sample_start, sample_start + check_result.checked_sample_len);
            }
            if (check_result.duplicate_start_idx > -1) {
                masked_starts.push_back(sample_start + check_result.duplicate_start_idx);
            }
            if (check_result.ambiguous_idx > -1) {
                ambiguous_positions.push_back(sample_start + check_result.ambiguous_idx);
            }
            if (check_result.ignored_region.first != -1) {
                result.skipped_regions.emplace_back(sample_start + check_result.ignored_region.first,
                                                   sample_start + check_result.ignored_region.second);
            }
            if (check_result.valid_sample_kmers) {
                for (const auto& kv : check_result.sample_seen_kmers) {
                    seen_kmers.insert(kv.first);
                }
            }
            sample_start = sample_start + check_result.next_start_offset;
        } else {
            // Prepare seq_info for this sample window
            size_t seqlen = sample_end - sample_start;
            std::vector<int> seq_info(seqlen, 0); // 0 = unannotated
            // Copy prior info into front
            for (size_t i = 0; i < sample_seq_info.size() && i < seq_info.size(); ++i) {
                seq_info[i] = sample_seq_info[i];
            }

            CheckOverallResult o = check_sample_overall(
                sample_seq,
                seq_info,
                seq_info_offset,
                seen_kmers,
                args.kmer,
                args.min_sample_len,
                args.no_overlap,
                args
            );

            if (o.checked_sample_len > -1) {
                result.sample_regions.emplace_back(sample_start, sample_start + o.checked_sample_len);
            }
            if (o.duplicate_start_idx > -1) {
                masked_starts.push_back(sample_start + o.duplicate_start_idx);
            }
            if (o.ambiguous_idx > -1) {
                ambiguous_positions.push_back(sample_start + o.ambiguous_idx);
            }
            if (o.ignored_region.first != -1) {
                result.skipped_regions.emplace_back(sample_start + o.ignored_region.first,
                                                   sample_start + o.ignored_region.second);
            }
            if (o.valid_sample_kmers) {
                for (const auto& kv : o.sample_seen_kmers) {
                    seen_kmers.insert(kv.first);
                }
            }

            sample_start = sample_start + o.next_start_offset;

            // Update per-sample seq_info tracking per Python logic
            size_t next_iteration_amount = std::min(static_cast<size_t>(args.sample_len), o.next_start_offset + static_cast<size_t>(args.kmer) - 1);
            seq_info_offset = o.seq_info_out.size() > next_iteration_amount ? (o.seq_info_out.size() - next_iteration_amount) : 0;
            // Slice sample_seq_info = sample_seq_info[next_start_offset:]
            sample_seq_info.clear();
            if (o.next_start_offset < o.seq_info_out.size()) {
                sample_seq_info.insert(sample_seq_info.end(), o.seq_info_out.begin() + static_cast<long>(o.next_start_offset), o.seq_info_out.end());
            }
        }
    }
    
    // Handle possible skipped sequence at the end
    // Workspace version: simplified - just add entire remaining region as skipped
    if (sample_start < seq.length()) {
        result.skipped_regions.emplace_back(sample_start, seq.length());
    }
    
    // Convert masked starting indices and ambiguous positions to regions
    // Workspace version: condense_masked_regions doesn't use k parameter
    result.masked_regions = condense_masked_regions(masked_starts);
    result.ambiguous_regions = condense_masked_regions(ambiguous_positions);
    
    return result;
}

/**
 * Write out BED file regions
 */
void writeout_bed(
    const std::string& name,
    const std::vector<Region>& regions,
    std::ofstream& bedfile
) {
    for (const auto& region : regions) {
        bedfile << name << "\t" << region.start << "\t" << region.end << "\n";
    }
}

/**
 * Output all bed files
 */
void output_dump(
    const std::vector<std::pair<std::string, DeduplicationResult>>& local_dict,
    const std::string& file_basename
) {
    std::ofstream sample_regions(file_basename + ".samples.bed");
    std::ofstream masked_regions(file_basename + ".masks.bed");
    std::ofstream skipped_regions(file_basename + ".ignored.bed");
    std::ofstream ambiguous_regions(file_basename + ".ambiguous.bed");
    
    for (const auto& entry : local_dict) {
        const std::string& seqname = entry.first;
        const DeduplicationResult& result = entry.second;
        
        writeout_bed(seqname, result.sample_regions, sample_regions);
        writeout_bed(seqname, result.masked_regions, masked_regions);
        writeout_bed(seqname, result.skipped_regions, skipped_regions);
        writeout_bed(seqname, result.ambiguous_regions, ambiguous_regions);
    }
}

/**
 * Simple FASTA parser with gzip support
 */
class FastaParser {
private:
    gzFile gz_file;
    std::string current_name;
    std::string current_seq;
    bool has_next;
    bool is_gzipped;
    std::string pending_line;
    bool has_pending;
    
    bool getline(std::string& line) {
        if (has_pending) {
            line = pending_line;
            has_pending = false;
            return true;
        }
        
        line.clear();
        char buffer[4096];
        
        if (gzgets(gz_file, buffer, sizeof(buffer)) != nullptr) {
            line = buffer;
            // Remove trailing newline
            if (!line.empty() && line.back() == '\n') {
                line.pop_back();
            }
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
            return true;
        }
        return false;
    }
    
    void putback_line(const std::string& line) {
        pending_line = line;
        has_pending = true;
    }
    
public:
    FastaParser(const std::string& filename) : has_next(true), has_pending(false) {
        // Check if file is gzipped
        is_gzipped = (filename.size() > 3 && filename.substr(filename.size() - 3) == ".gz");
        
        gz_file = gzopen(filename.c_str(), "rb");
        if (gz_file == nullptr) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        advance();
    }
    
    ~FastaParser() {
        if (gz_file != nullptr) {
            gzclose(gz_file);
        }
    }
    
    void advance() {
        current_seq.clear();
        current_name.clear();
        
        std::string line;
        bool found_header = false;
        
        while (getline(line)) {
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                if (found_header) {
                    // Put the line back for next iteration
                    putback_line(line);
                    has_next = true;
                    return;
                }
                current_name = line.substr(1);
                // Extract just the ID (first word)
                size_t space_pos = current_name.find(' ');
                if (space_pos != std::string::npos) {
                    current_name = current_name.substr(0, space_pos);
                }
                found_header = true;
            } else if (found_header) {
                current_seq += line;
            }
        }
        
        has_next = found_header && !current_seq.empty();
    }
    
    bool hasNext() const {
        return has_next;
    }
    
    std::string getName() const {
        return current_name;
    }
    
    std::string getSequence() const {
        return current_seq;
    }
};

/**
 * Clean sequence - convert to uppercase and replace ambiguous chars with N
 */
std::string clean_sequence(const std::string& seq) {
    std::string cleaned;
    cleaned.reserve(seq.length());
    
    for (char c : seq) {
        char upper_c = std::toupper(c);
        if (upper_c == 'A' || upper_c == 'C' || upper_c == 'G' || upper_c == 'T' || upper_c == 'N') {
            cleaned += upper_c;
        } else {
            cleaned += 'N';
        }
    }
    
    return cleaned;
}

/**
 * Estimate memory size of seen_kmers set
 * Rough approximation to match Python's sys.getsizeof
 */
size_t estimate_set_memory(const std::unordered_set<uint64_t>& seen_kmers) {
    // Python's sys.getsizeof includes object overhead + hash table overhead
    // Rough estimate: base overhead (232 bytes) + (element_count * 8 bytes per uint64_t * load_factor_overhead)
    // unordered_set typically has load factor around 1.0 and needs extra space for hash table
    return 232 + (seen_kmers.size() * 8 * 2); // multiply by 2 for hash table overhead
}

/**
 * Deduplicate a genome file
 */
void deduplicate_genome(
    const std::string& fasta,
    std::unordered_set<uint64_t>& seen_kmers,
    bool save_kmers_to_file,
    const Args& args
) {
    // Use vector of pairs to maintain insertion order
    std::vector<std::pair<std::string, DeduplicationResult>> local_dict;
    
    std::string fasta_basename = get_fasta_basename(fasta);
    std::string out_prefix = args.output_dir + "/" + fasta_basename;
    
    FastaParser parser(fasta);
    
    while (parser.hasNext()) {
        std::string seqname = parser.getName();
        std::string sequence = parser.getSequence();
        
        std::string clean_seq = clean_sequence(sequence);
        
        //std::cout << seqname << ": " << sequence.length() << " bp" << std::endl;
        
        DeduplicationResult result = deduplicate_seq(clean_seq, seen_kmers, args);
        
        // Workspace version: print memory estimate
        size_t mem_estimate = estimate_set_memory(seen_kmers);
        std::cout << "Seen kmer size: " << mem_estimate << std::endl;
        
        local_dict.push_back({seqname, result});
        
        parser.advance();
    }
    
    output_dump(local_dict, out_prefix);
    
    // Note: Saving kmers to binary file not implemented in this version
    // Would require serialization library like Boost.Serialization or custom implementation
}

/**
 * Main deduplication function
 */
void deduplicate(const Args& args) {
    //std::cout << "Starting deduplication..." << std::endl;
    
    // Read input files
    std::vector<std::string> fastas;
    
    // Check if input is a list file or direct fasta files
    if (args.input.size() == 1) {
        std::string input = args.input[0];
        if (input.find(".txt") != std::string::npos || input.find(".list") != std::string::npos) {
            // Read list file
            std::ifstream list_file(input);
            std::string line;
            while (std::getline(list_file, line)) {
                if (!line.empty()) {
                    fastas.push_back(line);
                }
            }
        } else {
            fastas = args.input;
        }
    } else {
        fastas = args.input;
    }
    
    // Filter to valid files
    std::vector<std::string> valid_fastas;
    for (const auto& fasta : fastas) {
        if (fs::exists(fasta)) {
            valid_fastas.push_back(fasta);
        } else {
            std::cerr << "Warning: could not find fasta: " << fasta << std::endl;
        }
    }
    
    // Write basename to file map
    std::string basename_fasta_file = args.output_dir + "/basename_fasta_match.txt";
    std::ofstream basename_file(basename_fasta_file);
    for (const auto& fasta : valid_fastas) {
        std::string basename = get_fasta_basename(fasta);
        basename_file << basename << "\t" << fasta << "\n";
    }
    basename_file.close();
    
    // Initialize seen kmers
    std::unordered_set<uint64_t> seen_kmers;
    // Note: Loading from pickle file would require custom deserialization
    
    // Process each fasta
    for (size_t i = 0; i < valid_fastas.size(); ++i) {
        const auto& fasta = valid_fastas[i];

        std::string fasta_basename = get_fasta_basename(fasta);
        std::string samples_bed = args.output_dir + "/" + fasta_basename + ".samples.bed";

        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        std::cout << "Working on " + fasta_basename + " at " + std::ctime(&now_time);
        
        if (std::filesystem::exists(samples_bed)) {

            std::cout << "Found existing output for " << fasta << "; skipping deduplication and computing seen kmers from existing samples" << std::endl;
            // Load existing samples and compute seen kmers
            int ret = system(("./code/kmc.sh " + samples_bed + " " + std::to_string(args.kmer) + " 16 64").c_str());
            if (ret != 0) {
                std::cerr << "Error: failed to compute seen kmers from existing samples for " << fasta << std::endl;
                exit(1);
            }

            // If ret is 0, the above was successful, meaning there is now a fasta_basename.kmers.txt file
            std::string kmers_filename = args.output_dir + "/kmc/" + fasta_basename + ".kmers.txt";
            std::ifstream file(kmers_filename);

            // Read kmers from file into seen_kmers
            if (!file.is_open()) {
                throw std::runtime_error("Cannot open file: " + kmers_filename);
            }
            std::string kmer;
            uint32_t count;
            while (file >> kmer >> count) {
                seen_kmers.insert(encode_kmer(kmer));
            }
            file.close();

            // Clean up by removing the kmers file
            std::filesystem::remove(kmers_filename);

        } else {

            std::cout << "Deduplicating " << fasta << std::endl;
            bool save_kmers = args.save_every > 0 && (i + 1) % args.save_every == 0;
            deduplicate_genome(fasta, seen_kmers, save_kmers, args);

        }

    }
    
    //std::cout << "Deduplication complete!" << std::endl;
}

/**
 * Print usage information
 */
void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [options] <input_files...>\n\n"
              << "Options:\n"
              << "  -e, --evaluation-method <str>  Evaluation method: per_kmer | per_sample_agnostic | per_sample_threshold (default: per_kmer)\n"
              << "  -a, --agnostic-retain-pct <float>  Likelihood to retain a sample with any duplication/ambiguity in per_sample_agnostic mode [0,1] (default: 0.0)\n"
              << "  -k, --kmer <int>           K-mer size [1,32] (default: 32)\n"
              << "  -l, --sample-len <int>     Sample length (default: 1000)\n"
              << "  -m, --min-sample-len <int> Minimum sample length (default: same as sample_len)\n"
              << "  -o, --output-dir <path>    Output directory (default: dedup_out/)\n"
              << "  -p, --seen-kmers <path>    Pickle file containing seen kmers (default: None)\n"
              << "  -r, --per-kmer-retain-pct <float>  Likelihood to allow duplicate kmer [0,1] in per_kmer mode (default: 0.0)\n"
              << "  -b, --ambiguous-base-threshold <float>  Allowed ambiguous base threshold in per_sample_threshold mode [0,1] (default: 0.0)\n"
              << "  -d, --duplicate-base-threshold <float>  Allowed duplicate base threshold in per_sample_threshold mode [0,1] (default: 0.0)\n"
              << "  -s, --save-every <int>     Save seen kmers every n samples (default: 0)\n"
              << "  --no-overlap               Keep neighboring samples discrete\n"
              << "  --no-save-kmers-at-end     Don't save kmers at program end\n"
              << "  --print-config             Print run arguments as JSON to stdout\n"
              << "  --no-write-config          Do not write run arguments to config.json\n"
              << "  --seed <int>               Random seed (default: 123)\n"
              << "  -h, --help                 Show this help message\n";
}

/**
 * Parse command line arguments
 */
Args parse_args(int argc, char* argv[]) {
    Args args;
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            exit(0);
        } else if (arg == "-e" || arg == "--evaluation-method") {
            if (i + 1 < argc) {
                args.evaluation_method = argv[++i];
            }
        } else if (arg == "-k" || arg == "--kmer") {
            if (i + 1 < argc) {
                args.kmer = std::stoi(argv[++i]);
            }
        } else if (arg == "-a" || arg == "--agnostic-retain-pct" || arg == "--agnostic_retain_pct") {
            if (i + 1 < argc) {
                args.agnostic_retain_pct = std::stod(argv[++i]);
            }
        } else if (arg == "-l" || arg == "--sample-len") {
            if (i + 1 < argc) {
                args.sample_len = std::stoi(argv[++i]);
            }
        } else if (arg == "-m" || arg == "--min-sample-len") {
            if (i + 1 < argc) {
                args.min_sample_len = std::stoi(argv[++i]);
            }
        } else if (arg == "-o" || arg == "--output-dir") {
            if (i + 1 < argc) {
                args.output_dir = argv[++i];
            }
        } else if (arg == "-p" || arg == "--seen-kmers") {
            if (i + 1 < argc) {
                args.seen_kmers_file = argv[++i];
            }
        } else if (arg == "-r" || arg == "--per-kmer-retain-pct" || arg == "--per_kmer_retain_pct") {
            if (i + 1 < argc) {
                args.per_kmer_retain_pct = std::stod(argv[++i]);
            }
        } else if (arg == "-b" || arg == "--ambiguous-base-threshold" || arg == "--ambiguous_base_threshold") {
            if (i + 1 < argc) {
                args.ambiguous_base_threshold = std::stod(argv[++i]);
            }
        } else if (arg == "-d" || arg == "--duplicate-base-threshold" || arg == "--duplicate_base_threshold") {
            if (i + 1 < argc) {
                args.duplicate_base_threshold = std::stod(argv[++i]);
            }
        } else if (arg == "-s" || arg == "--save-every") {
            if (i + 1 < argc) {
                args.save_every = std::stoi(argv[++i]);
            }
        } else if (arg == "--seed") {
            if (i + 1 < argc) {
                args.seed = std::stoi(argv[++i]);
            }
        } else if (arg == "--no-overlap") {
            args.no_overlap = true;
        } else if (arg == "--no-save-kmers-at-end") {
            args.no_save_kmers_at_end = true;
        } else if (arg == "--print-config") {
            args.print_config = true;
        } else if (arg == "--no-write-config") {
            args.write_config = false;
        } else if (arg[0] != '-') {
            args.input.push_back(arg);
        }
    }
    
    // Set min_sample_len to sample_len if not specified
    if (args.min_sample_len == -1) {
        args.min_sample_len = args.sample_len;
    }
    
    return args;
}

// Escape JSON string content
static std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default: out += c; break;
        }
    }
    return out;
}

// Serialize Args to JSON
static std::string args_to_json(const Args& args) {
    std::ostringstream oss;
    oss << "{\n";
    // input array
    oss << "  \"input\": [";
    for (size_t i = 0; i < args.input.size(); ++i) {
        if (i) oss << ", ";
        oss << '"' << json_escape(args.input[i]) << '"';
    }
    oss << "],\n";
    // scalars
    oss << "  \"kmer\": " << args.kmer << ",\n";
    oss << "  \"sample_len\": " << args.sample_len << ",\n";
    oss << "  \"min_sample_len\": " << args.min_sample_len << ",\n";
    oss << "  \"output_dir\": \"" << json_escape(args.output_dir) << "\",\n";
    oss << "  \"seen_kmers_file\": \"" << json_escape(args.seen_kmers_file) << "\",\n";
    oss << "  \"per_kmer_retain_pct\": " << args.per_kmer_retain_pct << ",\n";
    oss << "  \"agnostic_retain_pct\": " << args.agnostic_retain_pct << ",\n";
    oss << "  \"ambiguous_base_threshold\": " << args.ambiguous_base_threshold << ",\n";
    oss << "  \"duplicate_base_threshold\": " << args.duplicate_base_threshold << ",\n";
    oss << "  \"evaluation_method\": \"" << json_escape(args.evaluation_method) << "\",\n";
    oss << "  \"save_every\": " << args.save_every << ",\n";
    oss << "  \"no_overlap\": " << (args.no_overlap ? "true" : "false") << ",\n";
    oss << "  \"no_save_kmers_at_end\": " << (args.no_save_kmers_at_end ? "true" : "false") << ",\n";
    oss << "  \"seed\": " << args.seed << "\n";
    oss << "}";
    return oss.str();
}

/**
 * Main function
 */
int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            print_usage(argv[0]);
            return 1;
        }
        
        Args args = parse_args(argc, argv);
        
        // Validate arguments
        if (args.input.empty()) {
            std::cerr << "Error: No input files provided" << std::endl;
            return 1;
        }
        
        if (args.kmer < 1) {
            std::cerr << "Error: kmer size must be a positive integer" << std::endl;
            return 1;
        }
        
        // Validate evaluation method
        if (args.evaluation_method != "per_kmer" &&
            args.evaluation_method != "per_sample_agnostic" &&
            args.evaluation_method != "per_sample_threshold") {
            std::cerr << "Error: evaluation method must be one of per_kmer, per_sample_agnostic, per_sample_threshold" << std::endl;
            return 1;
        }
        // Only validate per_kmer retain pct if in per_kmer mode
        if (args.evaluation_method == "per_kmer" && (args.per_kmer_retain_pct < 0.0 || args.per_kmer_retain_pct > 1.0)) {
            std::cerr << "Error: per_kmer retain pct must be between 0.0 and 1.0" << std::endl;
            return 1;
        }
        if (args.evaluation_method == "per_sample_agnostic" && (args.agnostic_retain_pct < 0.0 || args.agnostic_retain_pct > 1.0)) {
            std::cerr << "Error: agnostic retain pct must be between 0.0 and 1.0" << std::endl;
            return 1;
        }
        if (args.evaluation_method == "per_sample_threshold" && (args.ambiguous_base_threshold < 0.0 || args.ambiguous_base_threshold > 1.0)) {
            std::cerr << "Error: ambiguous base threshold must be between 0.0 and 1.0" << std::endl;
            return 1;
        }
        if (args.evaluation_method == "per_sample_threshold" && (args.duplicate_base_threshold < 0.0 || args.duplicate_base_threshold > 1.0)) {
            std::cerr << "Error: duplicate base threshold must be between 0.0 and 1.0" << std::endl;
            return 1;
        }
        
        if (args.sample_len < args.kmer) {
            std::cerr << "Warning: sample_len is smaller than k, no deduplication will occur" << std::endl;
        }
        
        if (args.min_sample_len > args.sample_len) {
            std::cerr << "Warning: min_sample_len > sample_len, setting min_sample_len = sample_len" << std::endl;
            args.min_sample_len = args.sample_len;
        }
        
        if (args.kmer < 16) {
            std::cerr << "Warning: small kmer sizes will result in very strict deduplication" << std::endl;
        }
        
        // Workspace version: removed -y/--yes flag handling for directory overwrite prompt
        // Now always creates or uses existing directory without prompting
        
        // Create output directory if it doesn't exist
        if (!fs::exists(args.output_dir)) {
            fs::create_directories(args.output_dir);
        }

        // Output config: write to file and/or print to stdout
        const std::string json = args_to_json(args);
        if (args.write_config) {
            std::ofstream cfg(args.output_dir + "/config.json");
            cfg << json << std::endl;
        }
        if (args.print_config) {
            std::cout << json << std::endl;
        }
        
        // Initialize random number generator
        rng_gen.seed(args.seed);
        
        // Run deduplication
        deduplicate(args);
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
