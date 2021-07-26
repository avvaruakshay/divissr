/*
    A fastq object
    @file fastq_reader.h
    @author Akshay Avvaru
    @version 0.1 02/11/2020
*/

#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

namespace fastq {

    struct Read {
        uint count;
        string identifier;
        string sequence;
        string separator;
        string baseQual;
        bool valid=true;
        void reset() { identifier=""; sequence=""; separator=""; baseQual=""; }
    };
    
    /* 
     * An object for parsing fastq file 
     * Format of fastq files
     * Each read is represented by four lines 
     * 1. Identifier - Contains read name and information. Starts with "@". 
     * 2. Sequence   - Sequence of the read. 
     * 3. Separator  - Separator usually "+". 
     * 4. Qual       - Base call quality string. 
    */
    class Input {

        private:
            Read read;
            uint read_count = 0, bases = 0, total_reads = 0;
        public:
            uint line_count = 1;
        
        public:


            /*
            * Fetches a read data
            * @returns Read a information of the current read
            */
            Read fetch() {
                read.valid = false;
                bool line_read = false;
                uint line_count = 1;
                string line;
                while (line_count % 5 != 0 && std::getline(std::cin, line) ) {
                    switch (line_count) {
                        case 1: read.identifier = line; break;
                        case 2: read.sequence = line; break;
                        case 3: read.separator = line; break;
                        case 4: read.baseQual = line; break;
                    }
                    line_count++;
                    line_read = true;
                }
                if (read.identifier != "" && line_read) {
                    read.valid = true;
                    read_count++; bases += read.sequence.length();
                }
                return read;
            }

            /* Returns current read count */
            uint currentCount() { return read_count; }

            /* Returns total bases so far */
            uint currentBases() { return bases; }
    };
}