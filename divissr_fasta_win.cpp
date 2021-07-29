/*
    DiviSSR: DNA tandem repeat identification tool
    @file divissr_fasta_win.cpp
    @author Akshay Avvaru
    @version 0.1 06/08/2020
*/

#include <fstream>
#include <iostream>
#include <math.h>
#include <string.h>
#include <unordered_map>
#include <vector>
#include <chrono>
#include <assert.h>

#include "utils_win.h"

using namespace std;
using namespace std::chrono;

/* Main function of DiviSSR */
int main(int argc, char* argv[]) {
    ios_base::sync_with_stdio(false);

    string fin, fout, comp_fout;
    unsigned int m = 0, M = 0, cutoff = 0;
    int compound = 0, overlap_d = 0;
    int analyse_flag = 0;
    if (argc == 1) { utils::print_help(); exit (EXIT_FAILURE); }
    else if (argc > 1) { 
        utils::parse_arguments(argc, argv, fin, fout, m, M, cutoff, \
                                compound, overlap_d, comp_fout, analyse_flag);
        utils::length_cutoff_error(M, cutoff);
    }

    unsigned long long int gsize = 0, GC = 0;
    int sequences = 0;
    utils::count_seq(fin, sequences); // total number of sequences
    
    ifstream ins(fin); // input fasta file
    utils::input_file_error(ins.good(), fin);
    ofstream out(fout); // output file
    ofstream comp_out(comp_fout);

    unordered_map<string, string> rclass_map;
    string line;
    utils::bitSeqWindow window;
    utils::compoundRepeat compound_repeat;

    cout << endl << "Searching for tandem repeats in " << fin << endl;
    cout << "Min-motif: " << m << "\t Max-motif: " << M;
    cout << "\t Length-cutoff: " << cutoff <<  endl << endl;

    unsigned long long int start_time = duration_cast<milliseconds>(
        system_clock::now().time_since_epoch()
    ).count();

    // integer tracking the start of the repeat
    // -1 indicates no repeat is found 
    long long int start = -1; 
    long long int end;
    int rlen, atomicity;
    int window_repeat_check = 0;  // bool tracking if the window sequence is a repeat
    int numseq = 0; // current sequence number
    
    string seq_name, motif, repeat_class, strand;

    // NORM is require to fetch the current window sequence
    unsigned long long int const NORM = ~(0ull) >> 2*(32-cutoff);
    // non-redundant list of motifs used for checks
    vector<unsigned int> motif_checks = utils::get_motif_sizes(m, M);
    const unsigned int N = motif_checks.size();
    unsigned long long int divisor[N]; // list of divisors
    unsigned int rem_shift[N]; // list of remainder sizes
    for (int i=0; i<N; i++) {
        unsigned int d = cutoff / motif_checks[i];
        unsigned int r = cutoff % motif_checks[i];
        unsigned long long int D = 0ull;
        for (int j=0; j<d; j++) { D = D << (2*motif_checks[i]); D += 1; }
        D = D << (2*r);
        divisor[i] = D;
        rem_shift[i] = 2*(cutoff - r);
    }

    while(getline(ins, line)) {
        if (line[0] == '>') {
            float progress = ((float) numseq) / ((float) sequences);
            if (start != -1) {
                end = window.count; rlen = end - start;
                if (compound) {
                    compound_repeat.end = end;
                    if (compound_repeat.motif.size() > 1) {
                        compound_repeat.report();
                        comp_out << compound_repeat.output << '\n';
                    }
                }
                out << seq_name << "\t" << start << "\t" << end << "\t" \
                    << repeat_class << "\t" << rlen << "\t" \ 
                    << strand << "\t" << rlen/atomicity << "\t" << motif << '\n';
            }
            compound_repeat.reset();
            seq_name = line.substr(1, line.find(' ')-1);
            window.reset(); start = -1;
            
            utils::update_progress_bar(start_time, numseq, sequences);
            numseq++;
        }
        else {
            for(const auto c: line) {
                switch(c) {
                    case 'a': case 'A': break;
                    case 'c': case 'C': window.seq |= 1ull; GC += 1; break;
                    case 'g': case 'G': window.seq |= 2ull; GC += 1; break;
                    case 't': case 'T': window.seq |= 3ull; break;
                    default: 
                        window.seq = 0; window.cutoff = -1;
                        if (start != -1) {
                            end = window.count; rlen = end - start;
                            if (compound) { 
                                compound_repeat.end = end;
                                if (compound_repeat.motif.size() > 1) {
                                    compound_repeat.report();
                                    comp_out << compound_repeat.output << '\n';
                                }
                            }
                            out << seq_name << "\t" << start << "\t" << end << "\t" \
                                << repeat_class << "\t" << rlen << "\t" \ 
                                << strand << "\t" << rlen/atomicity << "\t" << motif << '\n';
                        }
                        compound_repeat.reset();
                        start = -1;
                        break;
                }
                gsize += 1;
                window.count += 1;
                window.cutoff += 1;
                window.seq &= NORM;

                if (window.cutoff >= cutoff) { // To be optimized
                    window_repeat_check = 0;
                    for (int i=0; i<N; ++i){
                        if ( (window.seq % divisor[i]) == ( window.seq >> rem_shift[i]) ) {
                            if (start == -1) {
                                atomicity = utils::check_atomicity(window.seq, cutoff, motif_checks[i]);
                                
                                // atomicity should be greater than 
                                // minimum motif-size
                                if (atomicity >= m) {

                                    start = window.count - cutoff;
                                    motif = utils::bit2base(window.seq, cutoff, atomicity);
                                    if (rclass_map.find(motif) != rclass_map.end()) { 
                                        repeat_class = rclass_map[motif];
                                    } else {
                                        repeat_class = utils::get_repeat_class(window.seq, cutoff, atomicity, rclass_map);
                                    }
                                    strand = repeat_class.substr(atomicity, 1);
                                    repeat_class = repeat_class.substr(0, atomicity);

                                    if (compound) {
                                        if ( start <= compound_repeat.end + overlap_d) {
                                            compound_repeat.overlap.push_back(start-compound_repeat.end);
                                        }
                                        else {
                                            if (compound_repeat.motif.size() > 1) {
                                                compound_repeat.report();
                                                comp_out << compound_repeat.output << '\n';
                                            }
                            
                                            compound_repeat.reset();
                                            compound_repeat.start = start;
                                            compound_repeat.seq_name = seq_name;
                                        }
                                        compound_repeat.repeat_class.push_back(repeat_class);
                                        compound_repeat.strand.push_back(strand);
                                        compound_repeat.motif.push_back(motif);
                                    }

                                }
                            }
                            window_repeat_check = 1; break;
                        }
                    }
                    if (window_repeat_check == 0 & start != -1) {
                        end = window.count - 1; rlen = end - start;
                        if (compound) {
                            compound_repeat.end = end;
                            compound_repeat.rlen.push_back(rlen);
                        }
                        out << seq_name << "\t" << start << "\t" << end << "\t" \
                            << repeat_class << "\t" << rlen << "\t" \ 
                            << strand << "\t" << rlen/atomicity << "\t" << motif << '\n';
                        start = -1;
                    }
                }
                window.seq <<= 2;
            }
        }
    }
    if (start != -1) {
        end = window.count; rlen = end - start;
        if (compound) {
            compound_repeat.end = end;
            if (compound_repeat.motif.size() > 1) {
                compound_repeat.report();
                comp_out << compound_repeat.output << '\n';
            }
        }
        compound_repeat.reset();
        out << seq_name << "\t" << start << "\t" << end << "\t" \
            << repeat_class << "\t" << rlen << "\t" \ 
            << strand << "\t" << rlen/atomicity << "\t" << motif << '\n';
    }

    if (analyse_flag) {
        float gc_percent = (float(GC) / float(gsize))*100;
        out << "#FileName: " << fin << '\n';
        out << "#GenomeSize: " << gsize << '\n';
        out << "#GC: " << gc_percent << '\n';
        out << "#NumSeq: " << sequences << '\n';
    }
    
    unsigned long long int end_time = duration_cast<milliseconds>(
        system_clock::now().time_since_epoch()
    ).count();
    float total_time = float(end_time - start_time)/1000.0;
    utils::update_progress_bar(start_time, numseq, sequences);

    ins.close(); out.close(); comp_out.close();
    return EXIT_SUCCESS;
}