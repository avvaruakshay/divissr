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

/* Main function of DiviSSR */
int main(int argc, char* argv[]) {
    ios_base::sync_with_stdio(false);

    string fout = argv[1];

    uint64_t gsize = 0, GC = 0, N_bases = 0;
    ofstream out(fout);

    unordered_map<string, string> rclass_map;
    
    string line;
    utils::bitSeqWindow window;
    utils::compoundRepeat compound_repeat;
    
    ofstream comp_out(argv[2]);
    
    $ fasta_gzip;

    $ python_input;

    cout << '\n' << "Searching for tandem repeats in " << fin << '\n';
    cout << "Min-motif: " << m << "\t Max-motif: " << M;
    cout << "\t Length-cutoff: " << cutoff <<  "\n\n";

    uint64_t start_time = duration_cast<milliseconds>(
        system_clock::now().time_since_epoch()
    ).count();
    

    // integer tracking the start of the repeat
    // -1 indicates no repeat is found 
    int start = -1; 
    int end, rlen, atomicity;
    int window_repeat_check = 0;  // bool tracking if the window sequence is a repeat
    int numseq = 0;    // current sequence number
    string seq_name, motif, repeat_class, strand;

    // NORM is require to fetch the current window sequence
    uint64_t const NORM = ~(0ull) >> 2*(32-cutoff);

    for (int i=0; i<2000; i++) { out << " "; }

    while(getline(cin, line)) {
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
                    compound_repeat.reset();
                }
                out << seq_name << "\t" << start << "\t" << end << "\t" \
                    << repeat_class << "\t" << rlen << "\t" \ 
                    << strand << "\t" << rlen/atomicity << "\t" << motif << '\n';
            }
            seq_name = line.substr(1, line.find(" ")-1);
            // cout << seq_name << endl;
            window.reset(); start = -1;

            utils::update_progress_bar(start_time, numseq, sequences);
            numseq++;
        }
        else {
            for(const auto c: line) {
                gsize += 1;
                switch(c) {
                    case 'a': case 'A': break;
                    case 'c': case 'C': window.seq |= 1ull; GC += 1; break;
                    case 'g': case 'G': window.seq |= 2ull; GC += 1; break;
                    case 't': case 'T': window.seq |= 3ull; break;
                    case 'N': case 'n': 
                        N_bases += 1;
                        window.seq = 0; window.cutoff = -1;
                        if (start != -1) {
                            end = window.count; rlen = end - start;
                            if (compound) { 
                                compound_repeat.end = end;
                                if (compound_repeat.motif.size() > 1) {
                                    compound_repeat.report();
                                    comp_out << compound_repeat.output << '\n';
                                }
                                compound_repeat.reset();
                            }
                            out << seq_name << "\t" << start << "\t" << end << "\t" \
                                << repeat_class << "\t" << rlen << "\t" \ 
                                << strand << "\t" << rlen/atomicity << "\t" << motif << '\n';
                        }
                        start = -1;
                        break;
                    default: continue;
                }
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
            compound_repeat.edn = end;
            if (compound_repeat.motif.size() > 1) {
                compound_repeat.report();
                comp_out << compound_repeat.output << '\n';
            }
            compound_repeat.reset();
        }
        out << seq_name << "\t" << start << "\t" << end << "\t" \
            << repeat_class << "\t" << rlen << "\t" \ 
            << strand << "\t" << rlen/atomicity << "\t" << motif << '\n';
    }
    
    string analyse_flag = argv[2];
    if (analyse_flag == "1") {
        float gc_percent = (float(GC) / float(gsize))*100;
        out << "#FileName: " << fin << '\n';
        out << "#GenomeSize: " << gsize << '\n';
        out << "#GC: " << gc_percent << '\n';
        out << "#NumSeq: " << sequences << '\n';
    }

    uint64_t end_time = duration_cast<milliseconds>(
        system_clock::now().time_since_epoch()
    ).count();
    float total_time = float(end_time - start_time)/1000.0;
    utils::update_progress_bar(start_time, numseq, sequences);

    out.close();
    return EXIT_SUCCESS;
}