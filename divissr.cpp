/*
    DiviSSR: DNA tandem repeat identification tool
    @file divissr.cpp
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

#include "utils.h"

using namespace std;
using namespace std::chrono;

// global variable
unordered_map<string, string> rClassMap;


/* Main function of LOOPER */
int main(int argc, char* argv[]) {
    ios_base::sync_with_stdio(false);

    string fin, fout;
    uint m = 0, M = 0, cutoff = 0;
    if (argc == 1) { utils::print_help(); exit (EXIT_FAILURE); }
    else if (argc > 1) { 
        utils::parse_arguments(argc, argv, fin, fout, m, M, cutoff);
        utils::length_cutoff_error(M, cutoff);
    }

    uint64_t gsize = 0, GC = 0;
    int sequences = 0;
    utils::count_seq(fin, sequences, gsize, GC); // total number of sequences
    ifstream ins(fin); // input fasta file
    utils::input_file_error(ins.good(), fin);
    ofstream out(fout); // output file
    for (int i=0; i<200; i++) { out << " "; }; out << "\n";
    string line;
    utils::bitSeqWindow window;

    cout << endl << "Searching for tandem repeats in " << fin << endl;
    cout << "Min-motif: " << m << "\t Max-motif: " << M;
    cout << "\t Length-cutoff: " << cutoff <<  endl << endl;

    uint64_t start_time = duration_cast<milliseconds>(
        system_clock::now().time_since_epoch()
    ).count();

    // integer tracking the start of the repeat
    // -1 indicates no repeat is found 
    int start = -1; 
    uint atomicity;           // atomicity of the repeat
    uint end;                 // end of the repeat    
    uint rlen;                // length of the repeat
    int repeat_check;         // track if the window sequence is a repeat
    int numseq          = 0;  // current sequence number
    int const BAR_WIDTH = 60; // progress of the progress bar
    string motif, repeat_class, seq_name;

    // NORM is require to fetch the current window sequence
    uint64_t const NORM = ~(0ull) >> 2*(32-cutoff);

    // non-redundant list of motifs used for checks
    vector<uint> motif_checks = utils::get_motif_sizes(m, M);

    const uint N = motif_checks.size();
    uint64_t divisor[N];    // list of divisors
    uint rem_shift[N];      // list of remainder sizes
    for (int i=0; i<N; i++) {
        uint d = cutoff / motif_checks[i];
        uint r = cutoff % motif_checks[i];
        uint64_t D = 0ull;
        for (int j=0; j<d; j++) { D = D << (2*motif_checks[i]); D += 1; }
        D = D << (2*r);
        divisor[i]   = D;
        rem_shift[i] = 2*(cutoff - r);
    }

    while(getline(ins, line)) {
        if (line[0] == '>') {
            float progress = ((float) numseq) / ((float) sequences);
            if (start != -1) {
                end = window.count; rlen = end - start;
                out << seq_name << "\t" << start << "\t" << end \
                << "\t" << repeat_class.substr(0, atomicity) << "\t" << rlen \
                << "\t" << repeat_class.substr(atomicity, 1) << "\t" \
                << rlen/atomicity << "\t" << motif << '\n';
            }
            seq_name = line.substr(1, line.find(' ')-1);
            window.reset(); start = -1;
            utils::update_progress_bar(start_time, numseq, sequences);
            cout.flush();
            numseq++;
        }
        else {
            for(const auto c: line) {
                switch(c) {
                    case 'a': case 'A': break;
                    case 'c': case 'C': window.seq |= 1ull; GC += 1; break;
                    case 'g': case 'G': window.seq |= 2ull; GC += 1; break;
                    case 't': case 'T': window.seq |= 3ull; break;
                    case 'N': case 'n': 
                        window.seq = 0; window.cutoff = -1;
                        if (start != -1) {
                            end = window.count; rlen = end - start;
                            out << seq_name << "\t" << start << "\t" << end \
                            << "\t" << repeat_class.substr(0, atomicity) << \
                            "\t" << rlen << "\t" << \
                            repeat_class.substr(atomicity, 1) << "\t" \
                            << rlen/atomicity << "\t" << motif << '\n';
                        }
                        start = -1;
                        break;
                    default: continue;
                }
                window.count    += 1;
                window.cutoff   += 1;
                window.seq &= NORM;
                gsize += 1;

                if (window.cutoff >= cutoff) { // To be optimized
                    repeat_check = 0;
                    for (int i=0; i<N; ++i){
                        if ( (window.seq % divisor[i]) == ( window.seq >> rem_shift[i]) ) {
                            if (start == -1) {
                                atomicity = utils::check_atomicity(window.seq, cutoff, motif_checks[i]);
                                
                                // atomicity should be greater than 
                                // minimum motif-size
                                if (atomicity >= m) {
                                    start = window.count - cutoff;
                                    motif = utils::bit2base(window.seq, cutoff, atomicity);
                                    if (rClassMap.find(motif) != rClassMap.end()) { 
                                        repeat_class = rClassMap[motif];
                                    } else {
                                        repeat_class = utils::get_repeat_class(window.seq, cutoff, atomicity, rClassMap);
                                    }
                                }
                            }
                            repeat_check = 1; break;
                        }
                    }
                    if (repeat_check == 0 & start != -1) {
                        end = window.count - 1; rlen = end - start;
                        out << seq_name << "\t" << start << "\t" << end \
                        << "\t" << repeat_class.substr(0, atomicity) << "\t" \
                        << rlen << "\t" << repeat_class.substr(atomicity, 1) \
                        << "\t" << rlen/atomicity << "\t" << motif << '\n';
                        start = -1;
                    }
                }
                window.seq <<= 2;
            }
        }
    }
    if (start != -1) {
        end = window.count; rlen = end - start;
        out << seq_name << "\t" << start << "\t" << end \
        << "\t" << repeat_class.substr(0, atomicity) << "\t" << rlen \
        << "\t" << repeat_class.substr(atomicity, 1) << "\t" \
        << rlen/atomicity << "\t" << motif << '\n';
    }

    uint64_t end_time = duration_cast<milliseconds>(
        system_clock::now().time_since_epoch()
    ).count();
    float total_time = float(end_time - start_time)/1000.0;
    utils::update_progress_bar(start_time, numseq, sequences);
    cout << "Total time taken: " << total_time << " seconds" << endl;

    out.clear(); out.seekp(0);
    out << "#FileName: " << fin << endl;
    out << "#GenomeSize: " << gsize << endl;
    out << "#GC%: " << (float(GC) / float(gsize))*100 << endl;
    out << "#NumSeq: " << sequences;
    ins.close(); out.close();
    return EXIT_SUCCESS;
}