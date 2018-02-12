//
//  evo_permute_codons.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 12/02/2018.
//  Copyright © 2018 Milan Malinsky. All rights reserved.
//

#include "evo_permute_codons.h"
#include "process_vcf_coding_sequences.h"
#include "process_vcf_utils.h"

#define SUBPROGRAM "permuteCodons"

#define DEBUG 1

static const char *PERMUTE_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] <-l list_of_multiple_aligment_files.txt>\n"
"Given multiple aligments of gene sequences (e.g. from the output of 'evo getCodingSeq')\n"
"this programs generates new alignment files with a random permutation of codons.\n"
"This can be used to calculate a null distribution of various statistics on coding sequences.\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -p,   --ploidy <d|h>                    d: diploid (default, sequences for haplotypes 1 and 2 are assumed to be interleaved in the alignment files); h: haploid\n"
"       -l,   --listOfFiles LIST.txt            a list with multiple alignment filenames, one per line (either -a or -l required)\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hp:a:l:t:n";

static const struct option longopts[] = {
    { "ploidy",   required_argument, NULL, 'p' },
    { "listOfFiles",   required_argument, NULL, 'l' },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string alignmentListFile = "";
    static char ploidy = 'd';
}


int permuteCodons(int argc, char** argv) {
    parsePermuteCodonsOptions(argc, argv);
    std::vector<string> allAligmentFiles; string statsFileName; string pcaVectorsFileName;
    
    std::cerr << "Permuting gene coding sequences" << std::endl;
    
    
    std::ifstream* alignmentList = new std::ifstream(opt::alignmentListFile.c_str());
    std::cerr << "for the genes in: " << opt::alignmentListFile << std::endl;
    std::string line;
    while (getline(*alignmentList, line)) {
        allAligmentFiles.push_back(line);
    }
    
        
    // Loop over the mutiple alignment files for the first time to load the codons:
    std::vector<string> allSeqs; std::vector<string> sampleNames; std::vector<int> geneLengths;
    for (std::vector<std::string>::size_type i = 0; i != allAligmentFiles.size(); i++) {
        std::ifstream* alignment = new std::ifstream(allAligmentFiles[i].c_str());
        std::string line; int lineNum = 1; int j = 0;
        while (getline(*alignment, line)) {
            if (lineNum % 2 == 1) assert(line[0] == '>');
            lineNum++;
            if (line[0] == '>') {
                if (i == 0)
                    sampleNames.push_back(line);
                continue;
            }
            if (i == 0) {
                allSeqs.push_back(line);
                allSeqs[j].reserve(1000000);
            } else {
                allSeqs[j] += line;
            }
            if (j == 0)
                geneLengths.push_back((int)line.length());
            j++;
        } alignment->close();
        
        assert(allSeqs[0].length() == vector_sum(geneLengths));
        std::cerr << "loaded seqs for: " << allAligmentFiles[i] << std::endl;
        std::cerr << "allSeqs[0].size(): " << allSeqs[0].length() << "; equals " << allSeqs[0].length()/3 << "AA" << std::endl;
    }
    
    // Loop over the mutiple alignment files for the second time to print out the permuted sequences:
    int totalAA = (int)allSeqs[0].length()/3;
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_int_distribution<> randomAApos(0, totalAA-1); // define the range
    for (std::vector<std::string>::size_type i = 0; i != allAligmentFiles.size(); i++) {
        std::string permutedFileName = stripExtension(allAligmentFiles[i].c_str()) + "_permuted.txt";
        std::ofstream* permutedFile = new std::ofstream(permutedFileName.c_str());
        int thisLengthAA = geneLengths[i]/3;
        std::vector<string> thisSeqs(sampleNames.size(),"");
        for (int j = 0; j < thisLengthAA; j++) {
            int pos = randomAApos(eng);
            for (int k = 0; k < allSeqs.size(); k++) {
                thisSeqs[k] += allSeqs[k].substr(pos*3,3);
            }
        }
        for (int k = 0; k < allSeqs.size(); k++) {
            *permutedFile << sampleNames[k] << std::endl;
            *permutedFile << thisSeqs[k] << std::endl;
            thisSeqs[k] = "";
        }
    }
    return 0;
}


void parsePermuteCodonsOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'p': arg >> opt::ploidy; break;
            case 'l': arg >> opt::alignmentListFile; break;
            case 'h':
                std::cout << PERMUTE_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (opt::alignmentListFile == "") {
        std::cerr << "The -l option must be specified\n";
        die = true;
    }
    
    
    if (opt::ploidy != 'h' && opt::ploidy != 'd') {
        std::cerr << "The -p (--ploidy) option can only have the values 'h' or 'd'\n";
        std::cerr << "At the moment I don't support other than haploid and diploid species\n";
        die = true;
    }
    
    if (argc - optind < 0) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 0)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << PERMUTE_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}
