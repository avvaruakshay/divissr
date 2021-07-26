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
        unsigned int count;
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
    class FastqFile {

        private:
            std::string file_name {};
            Read read;
            unsigned int read_count = 0, bases = 0, total_reads = 0;
        public:
            std::ifstream ins; // made public to check for input file error
        
        public:
            /*
            * Opens fastq_file and stores the stream in "ins"
            * @param name name of the fastq file
            */
            void parse(std::string name) {
                file_name = name;
                ins.open(file_name);
            }


            /*
            * Fetches a read data
            * @returns Read a information of the current read
            */
            Read fetch() {
                read.valid = false;
                if (ins.eof()) {
                    // std::cout << "End of the file reached! Total reads: " << read_count << '\n';
                    ins.close(); // close the file
                    // reset read to have empty members and return read
                    read.reset(); read.valid = false;
                }
                else {
                    std::getline(ins, read.identifier);
                    std::getline(ins, read.sequence);
                    std::getline(ins, read.separator);
                    std::getline(ins, read.baseQual);
                    if (read.identifier != "") {
                        read.valid = true;
                        read_count++; bases += read.sequence.length();
                    }
                }
                return read;
            }

            /* Returns current read count */
            unsigned int currentCount() { return read_count; }

            /* Returns total bases so far */
            unsigned int currentBases() { return bases; }
    };
}