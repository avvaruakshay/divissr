/*
 * Repeat identification in fastq files
 * @file: divissr_fastq.cpp
 * @author: Akshay Avvaru
 * @version: 0.1 05/11/2020
*/

#include <fstream>
#include <iostream>
#include <math.h>
#include <string.h>
#include <unordered_map>
#include <vector>

#include "utils.h"
#include "fastq_reader.h"

using namespace std;

unordered_map<string, string> rclass_map;
unordered_map<string, uint> repeats_count;
unordered_map<string, unordered_map<uint, uint>> repeats_length_freq;

/* 
 * Update repeats
 * @param repeat_class         Repeat class
 * @param rlen                 Length of the repeat
 * @param repeats_count        Reference for the unordered_map tracking repeats count
 * @param repeats_length_freq  Reference for the unordered_map tracking length wise frequency
*/
void update_result( string repeat_class, uint rlen,
                    unordered_map<string, uint> &repeats_count, 
                    unordered_map<string, unordered_map<uint, uint>> &repeats_length_freq )
{
    if (repeats_count.find(repeat_class) != repeats_count.end()) {
        repeats_count[repeat_class] += 1;
        if (repeats_length_freq[repeat_class].find(rlen) != repeats_length_freq[repeat_class].end()) {
            repeats_length_freq[repeat_class][rlen] += 1;
        }
        else {
            repeats_length_freq[repeat_class][rlen] = 1;
        }
    }
    else {
        repeats_count[repeat_class] = 1;
        repeats_length_freq[repeat_class][rlen] = 1;
    }
}


int main(int argc, char *argv[]) {
    ios_base::sync_with_stdio(false);

    string fin, fout;
    uint m = 0, M = 0, cutoff = 0, d = 0;
    if (argc == 1) { utils::print_help(); exit (EXIT_FAILURE); }
    else if (argc > 1) {
        utils::parse_arguments(argc, argv, fin, fout, m, M, cutoff);
        utils::length_cutoff_error(M, cutoff);
    }

    fastq::FastqFile ins; ins.parse(fin);
    utils::input_file_error(ins.ins.good(), fin);
    
    ofstream out(fout); 
    ofstream filterOut(fin + "_divissr.filtered.fastq");
    utils::bitSeqWindow window;

    cout << endl << "Searching for tandem repeats in " << fin << endl;
    cout << "Min-motif: " << m << "\t Max-motif: " << M;
    cout << "\t Length-cutoff: " << cutoff <<  endl;

    uint64_t start_time = duration_cast<milliseconds>(
        system_clock::now().time_since_epoch()
    ).count();

    // integer tracking the start of the repeat
    // -1 indicates no repeat is found 
    int start = -1; 
    uint atomicity, end, rlen;
    string motif, repeat_class, orientation, seq_name;
    int window_repeat_check = 0;  // current window repeat check

    // bitstring to retrieve window sequence with AND operation
    const uint64_t NORM = ~(0ull) >> 2*(32-cutoff); 

    // non-redundant list of motifs used for checks
    vector<uint> motif_checks = utils::get_motif_sizes(m, M);
    const uint N = motif_checks.size();
    uint extract_first[N];
    uint64_t extract_second[N];
    for (uint i = 0; i < N; i++) {
        extract_first[i] = 2 * motif_checks[i]; // shift twice of motif size
        extract_second[i] = NORM >> ( 2 * ( cutoff - motif_checks[i] ));
    }

    fastq::Read curr_read = ins.fetch();
    
    while(curr_read.valid) {
        window.reset();
        bool filter = false;
        for(const auto c: curr_read.sequence) {
            switch(c) {
                case 'a': case 'A': break;
                case 'c': case 'C': window.seq |= 1ull; break;
                case 'g': case 'G': window.seq |= 2ull; break;
                case 't': case 'T': window.seq |= 3ull; break;
                case 'N': case 'n': 
                    // resets the window when an N is encountered
                    window.seq = 0; window.cutoff = -1;
                    if (start != -1) {
                        filter = true;
                        end = window.count - 1; rlen = end - start;
                        update_result(repeat_class, rlen, repeats_count, repeats_length_freq);
                    }
                    start = -1;
                    break;
                default: continue;
            }
            window.count += 1;
            window.cutoff += 1;
            window.seq &= NORM;

            if (window.cutoff >= cutoff) {
                window_repeat_check = 0;
                for (int i=0; i<N; ++i){
                    for (int i=0; i<N; ++i){
                        if ( (window.seq % divisor[i]) == ( window.seq >> rem_shift[i]) ) {
                            if (start == -1) { 
                                filter = true;
                                start = window.count - cutoff; 
                                atomicity = utils::check_atomicity(window.seq, cutoff, motif_checks[i]);
                                motif = utils::bit2base(window.seq, cutoff, atomicity);
                                if (rclass_map.find(motif) != rclass_map.end()) {
                                    repeat_class = rclass_map[motif];
                                    orientation = repeat_class.substr(atomicity, 1);
                                    repeat_class = repeat_class.substr(0, atomicity);
                                }
                                else {
                                    repeat_class = utils::get_repeat_class(window.seq, cutoff, atomicity, rclass_map);
                                    orientation = repeat_class.substr(atomicity, 1);
                                    repeat_class = repeat_class.substr(0, atomicity);
                                }
                            }
                            window_repeat_check = 1; break;
                        }
                    }
                }
                if (window_repeat_check == 0 & start != -1) {
                    filter = true;
                    end = window.count - 1; rlen = end - start;
                    update_result(repeat_class, rlen, repeats_count, repeats_length_freq);
                    start = -1;
                }
            }
            window.seq <<= 2;
        }

        // if reat ends in a repeat
        if (start != -1) {
            filter = true;
            end = window.count - 1; rlen = end - start;
            update_result(repeat_class, rlen, repeats_count, repeats_length_freq);
            start = -1;
        }

        // if repeat is found put read in filtered fastq file
        if (filter) {
            filterOut << curr_read.identifier << "\n";
            filterOut << curr_read.sequence   << "\n";
            filterOut << curr_read.separator  << "\n";
            filterOut << curr_read.baseQual   << "\n";
        }
        curr_read = ins.fetch();
    }


    out << "#Total reads: " << ins.currentCount() << "\n";
    for (auto it = repeats_count.begin(); it != repeats_count.end(); it++) {
        uint reads = ins.currentCount();
        out << it->first << "\t" << it->second << "\t" << (float(it->second)/float(reads))*100000 << "\t";
        for (auto jt = repeats_length_freq[it->first].begin(); jt != repeats_length_freq[it->first].end(); jt++) {
            out << jt->first << "-" << jt->second << ";";
        }
        out << "\n";
    }

    out.close();
    filterOut.close();
    return EXIT_SUCCESS;
}