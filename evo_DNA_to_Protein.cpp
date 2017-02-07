//
//  evo_DNA_to_Protein.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 01/02/2017.
//  Copyright © 2017 Milan Malinsky. All rights reserved.
//

#include "evo_DNA_to_Protein.h"

#define SUBPROGRAM "DNAtoProtein"

#define DEBUG 1

static const char *DNATOPROTEIN_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] <-a multiple_alignment.fa | -l list_of_multiple_aligment_files.txt>\n"
"Translate multiple aligments of gene sequences (e.g. from the output of 'evo getCodingSeq')\n"
"from DNA into protein, optionally editing labels to work with PHYLIP and Haploviewer\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -a,   --alignment FILE.fa               a multiple alignment file (either -a or -l required)\n"
"       -l,   --listOfFiles LIST.txt            a list with multiple alignment filenames, one per line (either -a or -l required)\n"
"       --newLabels=FILE_WITH_LABELS.txt        one per line, with as many records as there are in the fasta files\n"
"       --outFolder=FOLDER                      put the output in FOLDER (default .)\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "ha:l:";

enum { OPT_NEW_LABELS, OPT_FOLDER };

static const struct option longopts[] = {
    { "alignment",   required_argument, NULL, 'a' },
    { "listOfFiles",   required_argument, NULL, 'l' },
    { "newLabels", required_argument, NULL, OPT_NEW_LABELS },
    { "outFolder", required_argument, NULL, OPT_FOLDER },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string alignmentFile = "";
    static string alignmentListFile = "";
    static string newLabelFile = "";
    static string outFolder = "";
}


int DNAtoProtein(int argc, char** argv) {
    parseDNAtoProteinOptions(argc, argv);
    std::vector<string> allAligmentFiles; string AAFileName;
    std::vector<string> allNewLabels;
    
    std::cerr << "Translate multiple aligments from DNA to protein" << std::endl;
    
    if (opt::alignmentFile != "") {
        allAligmentFiles.push_back(opt::alignmentFile);
        // look at suffix
        size_t suffixPos = opt::alignmentFile.find_last_of('.');
        if (suffixPos !=  std::string::npos) {
            if (opt::alignmentFile.substr(suffixPos) == ".fa" || opt::alignmentFile.substr(suffixPos) == ".fasta")
                std::cerr << "for the gene: " << stripExtension(opt::alignmentFile) << std::endl;
            else
                std::cerr << "for the gene: " << opt::alignmentFile << std::endl;
        } else {
            std::cerr << "for the gene: " << opt::alignmentFile << std::endl;
        }
    } else {
        std::ifstream* alignmentList = new std::ifstream(opt::alignmentListFile.c_str());
        std::cerr << "for the genes in: " << opt::alignmentListFile << std::endl;
        std::string line;
        while (getline(*alignmentList, line)) {
            allAligmentFiles.push_back(line);
        }
    }
    if (opt::newLabelFile != "") {
        std::cerr << "using labels from: " << opt::newLabelFile << std::endl;
        std::ifstream* newLabels = new std::ifstream(opt::newLabelFile.c_str());
        std::string line;
        while (getline(*newLabels, line)) {
            allNewLabels.push_back(line);
        }
    }
    
    // Loop over the mutiple alignment files:
    for (std::vector<std::string>::size_type i = 0; i != allAligmentFiles.size(); i++) {
        std::ifstream* alignment = new std::ifstream(allAligmentFiles[i].c_str());
        std::vector<string> allSeqs; std::vector<string> allSeqsAA;
        // Load the DNA alignment
        std::string line; int lineNum = 1;
        while (getline(*alignment, line)) {
            if (lineNum % 2 == 1) assert(line[0] == '>');
            lineNum++;
            if (line[0] == '>' && opt::newLabelFile != "") continue;
            if (line[0] == '>' && opt::newLabelFile == "") { allNewLabels.push_back(line); continue; }
            allSeqs.push_back(line);
        } alignment->close();
        
       // std::cerr << "loaded seqs for: " << allAligmentFiles[i] << std::endl;
        //std::cerr << "allSeqs.size(): " << allSeqs.size() << std::endl;
        // Translate to Amino Acid sequences:
        std::string::size_type geneLengthNt = allSeqs[0].length();
        int numSamples = (int)allSeqs.size();
        std::vector<string> altCodons; altCodons.resize(numSamples);
        allSeqsAA.resize(numSamples);
        if (numSamples > 0) {
            assert(geneLengthNt % 3 == 0); // The gene length must be divisible by three
            for (std::vector<std::string>::size_type j = 0; j != numSamples; j++) {
                allSeqsAA[j].reserve(geneLengthNt/3);
            }
            for (string::size_type k = 0; k != geneLengthNt; k++) {
                for (std::vector<std::string>::size_type j = 0; j != numSamples; j++) {
                    altCodons[j] += allSeqs[j][k];
                }
                // Find the types of mutation we are dealing with
                if ((k+1)%3 == 0) {
                    for (std::vector<std::string>::size_type j = 0; j != numSamples; j++) {
                        assert(allSeqs[j].substr(k-2,3) == altCodons[j]);
                        allSeqsAA[j] += getAminoAcidOneLetterCode(altCodons[j]);
                    }
                    for (std::vector<std::string>::size_type j = 0; j != numSamples; j++) {
                        altCodons[j] = "";
                    }
                }
            }
        } else{
            std::cerr << allAligmentFiles[i] << "has no alignment.." << std::endl;
        }
        //std::cerr << "processed seqs for: " << allAligmentFiles[i] << std::endl;
        //std::cerr << "going to output: " << allAligmentFiles[i] << std::endl;
        size_t suffixPos = allAligmentFiles[i].find_last_of('.');
        if (suffixPos != std::string::npos) {
            if (allAligmentFiles[i].substr(suffixPos) == ".fa" || allAligmentFiles[i].substr(suffixPos) == ".fasta")
                AAFileName = opt::outFolder + stripExtension(allAligmentFiles[i]) + "_AA.fasta";
            else
                AAFileName = opt::outFolder + allAligmentFiles[i] + "_AA.fasta";
        } else {
            AAFileName = opt::outFolder + allAligmentFiles[i] + "_AA.fasta";
        }
        //std::cerr << "into: " << AAFileName << std::endl;
        std::ofstream* AAalignment = new std::ofstream(AAFileName.c_str());
        for (std::vector<std::string>::size_type j = 0; j != numSamples; j++) {
            *AAalignment << allNewLabels[j] << std::endl;
            *AAalignment << allSeqsAA[j] << std::endl;
        } AAalignment->close();
    }
    return 0;
}


void parseDNAtoProteinOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'a': arg >> opt::alignmentFile; break;
            case 'l': arg >> opt::alignmentListFile; break;
            case OPT_NEW_LABELS: arg >> opt::newLabelFile; break;
            case OPT_FOLDER: arg >> opt::outFolder; break;
            case 'h':
                std::cout << DNATOPROTEIN_USAGE_MESSAGE;
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
        std::cout << "\n" << DNATOPROTEIN_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}
