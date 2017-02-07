//
//  evo_codingStats_from_alignment.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 17/01/2017.
//  Copyright © 2017 Milan Malinsky. All rights reserved.
//

#include "evo_codingStats_from_alignment.h"
#include "process_vcf_coding_sequences.h"


#define SUBPROGRAM "codingStats"

#define DEBUG 1

static const char *CODINGSTATS_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] <-a multiple_alignment.fa | -l list_of_multiple_aligment_files.txt>\n"
"Calculate statistics out of multiple aligments of gene sequences (e.g. from the output of 'evo getCodingSeq')\n"
"No checks are performed on the sequences\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -p,   --ploidy <d|h>                    d: diploid (default, sequences for haplotypes 1 and 2 are assumed to be interleaved in the alignment files); h: haploid\n"
"       -t,   --tStV RATIO                      observed genome-wide tS/tV ratio in the dataset\n"
"       -a,   --alignment FILE.fa               a multiple alignment file (either -a or -l required)\n"
"       -l,   --listOfFiles LIST.txt            a list with multiple alignment filenames, one per line (either -a or -l required)\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hp:a:l:t:";



static const struct option longopts[] = {
    { "ploidy",   required_argument, NULL, 'p' },
    { "alignment",   required_argument, NULL, 'a' },
    { "tStV",   required_argument, NULL, 't' },
    { "listOfFiles",   required_argument, NULL, 'l' },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string alignmentFile = "";
    static string alignmentListFile = "";
    static char ploidy = 'd';
    static double tStVratio = 0.5;      // Under equal mutation probability,
                                        //half transitions will be observed compared with transversions
}


int getCodingStats(int argc, char** argv) {
    parseCodingStatsOptions(argc, argv);
    std::vector<string> allAligmentFiles; string statsFileName;
    
    std::cerr << "Calculating gene coding statistics" << std::endl;
    
    if (opt::alignmentFile != "") {
        allAligmentFiles.push_back(opt::alignmentFile);
        // look at suffix
        size_t suffixPos = opt::alignmentFile.find_last_of('.');
        if (suffixPos !=  std::string::npos) {
            if (opt::alignmentFile.substr(suffixPos) == ".fa" || opt::alignmentFile.substr(suffixPos) == ".fasta") {
                statsFileName = stripExtension(opt::alignmentFile) + "_stats.txt";
                std::cerr << "for the gene: " << stripExtension(opt::alignmentFile) << std::endl;
            } else {
                statsFileName = opt::alignmentFile + "_stats.txt";
                std::cerr << "for the gene: " << opt::alignmentFile << std::endl;
            }
        }
    } else {
        std::ifstream* alignmentList = new std::ifstream(opt::alignmentListFile.c_str());
        statsFileName = stripExtension(opt::alignmentListFile) + "_stats.txt";
        std::cerr << "for the genes in: " << opt::alignmentListFile << std::endl;
        std::string line;
        while (getline(*alignmentList, line)) {
            allAligmentFiles.push_back(line);
        }
    }
    
    std::ofstream* statsFile = new std::ofstream(statsFileName.c_str());
    if (opt::ploidy == 'd') {
        std::cout << "transcript" << "\t" << "ntLengh" << "\t" << "pN" << "\t" << "pS" << "\t" << "hetN" << "\t" << "hetS" << "\t" << "pNstdErr" << "\t" << "pSstdErr" << "\t" << "pN-pS_stdErr" << std::endl;
        *statsFile << "transcript" << "\t" << "ntLengh" << "\t" << "pN" << "\t" << "pS" << "\t" << "hetN" << "\t" << "hetS" << "\t" << "pNstdErr" << std::endl;
    } else {
        std::cout << "transcript" << "\t" << "ntLengh" << "\t" << "pN" << "\t" << "pS" << std::endl;
        *statsFile << "transcript" << "\t" << "ntLengh" << "\t" << "pN" << "\t" << "pS" << std::endl;
    }
    // Loop over the mutiple alignment files:
    for (std::vector<std::string>::size_type i = 0; i != allAligmentFiles.size(); i++) {
        std::ifstream* alignment = new std::ifstream(allAligmentFiles[i].c_str());
        std::vector<string> allSeqs; std::vector<string> allSeqsH2;
        std::string line; int lineNum = 1;
        while (getline(*alignment, line)) {
            if (lineNum % 2 == 1) assert(line[0] == '>');
            lineNum++;
            if (line[0] == '>') continue;
            if (opt::ploidy == 'd' && allSeqsH2.size() == (allSeqs.size() - 1))
                allSeqsH2.push_back(line);
            else
                allSeqs.push_back(line);
        } alignment->close();
        
        //std::cerr << "loaded seqs for: " << allAligmentFiles[i] << std::endl;
        //std::cerr << "allSeqs.size(): " << allSeqs.size() << std::endl;
        if (allSeqs.size() > 0) {
            assert(allSeqs[0].length() % 3 == 0); // The gene length must be divisible by three
            if (opt::ploidy == 'd')
                assert(allSeqs[0].length() == allSeqsH2[0].length());
            
            std::vector<string> statsThisGene; statsThisGene.push_back(allAligmentFiles[i]);
            if (opt::ploidy == 'd') {
                // std::cerr << "getting stats for: " << allAligmentFiles[i] << std::endl;
                getStatsBothPhasedHaps(allSeqs, allSeqsH2, statsThisGene, opt::tStVratio);
                print_vector_stream(statsThisGene, std::cout);
                print_vector(statsThisGene, *statsFile);
            } else {
                std::vector<double> pNpS = getPhasedPnPs(allSeqs);
                assert(pNpS.size() == 2);
                statsThisGene.push_back(numToString(allSeqs[0].length()));
                statsThisGene.push_back(numToString(pNpS[0]));
                statsThisGene.push_back(numToString(pNpS[1]));
                print_vector_stream(statsThisGene, std::cout);
                print_vector(statsThisGene, *statsFile);
            }
        } else{
            std::cout << allAligmentFiles[i] << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << std::endl;
            *statsFile << allAligmentFiles[i] << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << std::endl;
        }
    }
    return 0;
}


void parseCodingStatsOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'p': arg >> opt::ploidy; break;
            case 'a': arg >> opt::alignmentFile; break;
            case 'l': arg >> opt::alignmentListFile; break;
            case 't': arg >> opt::tStVratio; break;
            case 'h':
                std::cout << CODINGSTATS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (opt::alignmentFile == "" && opt::alignmentListFile == "") {
        std::cerr << "Either -a or -l options must be specified\n";
        die = true;
    }
    
    if (opt::alignmentFile != "" && opt::alignmentListFile != "") {
        std::cerr << "The -a and -l options can't both be specified\n";
        die = true;
    }
    
    if (opt::ploidy != 'h' && opt::ploidy != 'd') {
        std::cerr << "The -p (--ploidy) option can only have the values 'h' or 'd'\n";
        std::cerr << "At the moment I don't support other than haploid and diploid species\n";
        die = true;
    }
    
    if (opt::tStVratio <= 0) {
        std::cerr << "The tS/tV ratio must be positive\n";
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
        std::cout << "\n" << CODINGSTATS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

