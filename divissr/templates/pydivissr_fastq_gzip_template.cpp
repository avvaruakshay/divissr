/*
 * Repeat identification in fastq files
 * @file: hamSTR_fastq.cpp
 * @author: Akshay Avvaru
 * @version: 0.1 05/11/2020
*/

#include <fstream>
#include <iostream>
#include <math.h>
#include <string.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <chrono>
#include <bits/stdc++.h> 

#include "utils.h"
#include "fastq_gzip_reader.h"

using namespace std;
using namespace std::chrono;

unordered_map<string, string> rclass_map;
unordered_map<string, uint> repeats_count, repeats_reads, repeats_bases;
unordered_map<string, unordered_map<uint, uint>> repeats_length_freq;
unordered_map<uint, uint> read_length_dist;
uint total_repeats, total_repeat_reads, total_repeat_bases;

/* 
 * Update repeats
 * @param repeat_class         Repeat class
 * @param rlen                 Length of the repeat
 * @param repeats_count        Reference for the unordered_map tracking repeats count
 * @param repeats_length_freq  Reference for the unordered_map tracking length wise frequency
*/
void update_result( string repeat_class, uint rlen,
                    unordered_map<string, uint> &repeats_count, unordered_map<string, uint> &repeats_bases,
                    unordered_map<string, unordered_map<uint, uint>> &repeats_length_freq )
{
    if (repeats_count.find(repeat_class) != repeats_count.end()) {
        repeats_count[repeat_class] += 1;
        repeats_bases[repeat_class] += rlen;
        if (repeats_length_freq[repeat_class].find(rlen) != repeats_length_freq[repeat_class].end()) {
            repeats_length_freq[repeat_class][rlen] += 1;
        }
        else {
            repeats_length_freq[repeat_class][rlen] = 1;
        }
    }
    else {
        repeats_count[repeat_class] = 1;
        repeats_bases[repeat_class] = rlen;
        repeats_length_freq[repeat_class][rlen] = 1;
    }
}


int main(int argc, char *argv[]) {
    ios_base::sync_with_stdio(false);

    string fout = argv[1];
    ofstream filterOut;
    bool filter_reads = false;
    if (argv[2]) {
        filter_reads = true;
        filterOut.open(argv[2]);
    }

    fastq::Input ins = fastq::Input();
    
    ofstream out(fout); 
    utils::bitSeqWindow window;

    $ python_input;

    // cout << endl << "Searching for tandem repeats in " << fin << endl;
    cout << "Min-motif: " << m << "\t Max-motif: " << M;
    cout << "\t Length-cutoff: " << cutoff <<  endl;

    uint64_t start_time = duration_cast<milliseconds>(
        system_clock::now().time_since_epoch()
    ).count();

    // integer tracking the start of the repeat
    // -1 indicates no repeat is found 
    int start = -1; 
    uint end, rlen, atomicity;
    int window_repeat_check = 0;  // bool tracking if the window sequence is a repeat
    string seq_name, motif, repeat_class, strand;

    // bitstring to retrieve window sequence with AND operation
    const uint64_t NORM = ~(0ull) >> 2*(32-cutoff); 
    fastq::Read curr_read = ins.fetch();
    
    while(curr_read.valid) {
        window.reset();
        bool filter = false;
        uint read_length = curr_read.sequence.length();
        if (read_length_dist.find(read_length) != read_length_dist.end()) {
            read_length_dist[read_length] += 1;
        } else { read_length_dist[read_length] = 1; }
        vector<string> read_repeats;

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
                        end = window.count; rlen = end - start;
                        total_repeat_bases += rlen;
                        update_result(repeat_class, rlen, repeats_count, repeats_bases, repeats_length_freq);
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
                                atomicity = utils::check_atomicity(window.seq, cutoff, motif_checks[i]);

                                // atomicity should be greater than 
                                // minimum motif-size
                                if (atomicity >= m) {
                                    filter = true;
                                    total_repeats += 1;
                                    start = window.count - cutoff; 
                                    motif = utils::bit2base(window.seq, cutoff, atomicity);
                                    if (rclass_map.find(motif) != rclass_map.end()) {
                                        repeat_class = rclass_map[motif];
                                    } else {
                                        repeat_class = utils::get_repeat_class(window.seq, cutoff, atomicity, rclass_map);
                                    }
                                    repeat_class = repeat_class.substr(0, atomicity);
                                    vector<string>::iterator found = find(read_repeats.begin(), read_repeats.end(), repeat_class);
                                    if (found == read_repeats.end()) { read_repeats.push_back(repeat_class); }
                                    strand = repeat_class.substr(atomicity, 1);
                                }
                            }
                            window_repeat_check = 1; break;
                        }
                    }
                }
                if (window_repeat_check == 0 & start != -1) {
                    filter = true;
                    end = window.count - 1; rlen = end - start;
                    total_repeat_bases += rlen;
                    update_result(repeat_class, rlen, repeats_count, repeats_bases, repeats_length_freq);

                    start = -1;
                }
            }
            window.seq <<= 2;
        }

        // if read ends in a repeat
        if (start != -1) {
            filter = true;
            end = window.count; rlen = end - start;
            total_repeat_bases += rlen;
            update_result(repeat_class, rlen, repeats_count, repeats_bases, repeats_length_freq);
            start = -1;
        }

        // if repeat is found put read in filtered fastq file
        if (filter) {
            total_repeat_reads += 1;
            if (filter_reads) {
                filterOut << curr_read.identifier << "\n";
                filterOut << curr_read.sequence   << "\n";
                filterOut << curr_read.separator  << "\n";
                filterOut << curr_read.baseQual   << "\n";
            }
        }

        if (ins.currentCount() % 50000 == 0 ) {
            uint64_t now = duration_cast<milliseconds>(
                system_clock::now().time_since_epoch()
            ).count();
            float total_time = float(now-start_time)/1000.0;

            std::cout << fixed << setprecision(2) <<"Time elapsed: " << total_time << "secs|";
            std::cout << fixed << setprecision(2) <<"Reads processed: " << ins.currentCount() << "|";
            std::cout << fixed << setprecision(2) <<"Reads per sec: " << ins.currentCount()/total_time;
            std::cout << "          \r";

            std::cout.flush();
        }

        for(vector<string>::iterator it=read_repeats.begin(); it!=read_repeats.end(); it++) {
            if (repeats_reads.find(*it) != repeats_reads.end()) { repeats_reads[*it] += 1; }
            else { repeats_reads[*it] = 1; }
        }

        curr_read = ins.fetch();
    }

    uint64_t now = duration_cast<milliseconds>(
        system_clock::now().time_since_epoch()
    ).count();
    float total_time = float(now-start_time)/1000.0;
    uint total_reads = ins.currentCount();
    uint total_bases = ins.currentBases();

    std::cout << fixed << setprecision(2) << "Time elapsed: " << total_time << "secs|";
    std::cout << fixed << setprecision(2) << "Reads processed: " << total_reads << "|";
    std::cout << fixed << setprecision(2) << "Reads per sec: " << total_reads/total_time << "\n";


    out << "#TotalReads: "          << total_reads          << "\n";
    out << "#TotalBases: "          << total_bases          << "\n";
    out << "#RepeatReads: "         << total_repeat_reads   << "\n";
    out << "#RepeatBases: "         << total_repeat_bases   << "\n";
    out << "#TotalRepeats: "        << total_repeats        << "\n";
    out << "#PercentRepeatReads: "  << (float(total_repeat_reads)/float(total_reads))*100 << "\n";
    out << "#NumRepClasses: "       << repeats_count.size() << "\n";

    out << "#ReadLengthDist: ";
    vector<uint> read_lengths;
    for(auto it = read_length_dist.begin(); it != read_length_dist.end(); it++) {
        uint length = it->first;
        read_lengths.push_back(length);
    }
    sort(read_lengths.begin(), read_lengths.end(), [](uint a, uint b) { return a < b; });
    for (auto l = read_lengths.begin(); l != read_lengths.end(); l++) {
        if (l == read_lengths.end() - 1) {
            out << *l << '-' << read_length_dist[*l];
        } else {
            out << *l << '-' << read_length_dist[*l] << ';';
        }
    }
    out << "\n";

    for (auto it = repeats_count.begin(); it != repeats_count.end(); it++) {
        string rclass = it->first;
        uint freq = it->second;
        uint rreads = repeats_reads[rclass];
        uint rbases = repeats_bases[rclass];
        out << rclass << "\t" << freq << "\t" << rreads << "\t" \
            << rbases << "\t";
        vector<uint> rlengths;
        for(auto it=repeats_length_freq[rclass].begin(); it!=repeats_length_freq[rclass].end(); it++) {
            rlengths.push_back(it->first);
        }
        sort(rlengths.begin(), rlengths.end(), [](uint a, uint b) { return a < b; });
        for (auto jt = rlengths.begin(); jt != rlengths.end(); jt++) {
            if (jt == rlengths.end() - 1) {
                out << *jt << "-" << repeats_length_freq[rclass][*jt];
            } else {
                out << *jt << "-" << repeats_length_freq[rclass][*jt] << ";";
            }
        }
        out << "\n";
    }

    out.close();
    if (filter_reads) { filterOut.close(); }
    return EXIT_SUCCESS;
}