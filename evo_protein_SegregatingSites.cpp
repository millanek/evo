//
//  evo_protein_SegregatingSites.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 03/02/2017.
//  Copyright © 2017 Milan Malinsky. All rights reserved.
//

#include "evo_protein_SegregatingSites.h"

#define SUBPROGRAM "ProteinSs"

#define DEBUG 1

static const char *DNATOPROTEIN_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] <-a multiple_alignment.fa | -l list_of_multiple_aligment_files.txt>\n"
"Output segregating sites from a (protein) multiple alignment\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -a,   --alignment FILE.fa               a multiple alignment file (either -a or -l required)\n"
"       -l,   --listOfFiles LIST.txt            a list with multiple alignment filenames, one per line (either -a or -l required)\n"
"       --hapLabels=FILE_WITH_LABELS.txt        one per line, with as many records as there are haplotypes\n"
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


int ProteinSs(int argc, char** argv) {
    parseProteinSsOptions(argc, argv);
    std::vector<string> allAligmentFiles; string AAFileName;
    std::vector<string> allNewLabels;
    
    
    std::cerr << "Outputting segregating sites" << std::endl;
    
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
        std::vector<int> ssPos;
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
        std::string::size_type geneLength = allSeqs[0].length();
        int numSamples = (int)allSeqs.size();
        std::vector<string> altCodons; altCodons.resize(numSamples);
        allSeqsAA.resize(numSamples);
        std::vector<string> uniqueHpt; std::vector<string> uniqueNames;
        std::vector<string> otherSamples;
        if (numSamples > 0) {
            for (std::vector<std::string>::size_type j = 0; j != numSamples; j++) {
                allSeqsAA[j].reserve(geneLength);
            }
            for (string::size_type k = 0; k != geneLength; k++) {
                // Find if it's a segregating site:
                bool segregating = false;
                for (std::vector<std::string>::size_type j = 0; j != numSamples - 1; j++) {
                    for (std::vector<std::string>::size_type l = j + 1; l != numSamples; l++) {
                        if (allSeqs[j][k] != allSeqs[l][k]) {
                            segregating = true; break;
                        }
                    }
                    if (segregating) break;
                }
                if (segregating) {
                    ssPos.push_back((int)k+1);
                    for (std::vector<std::string>::size_type j = 0; j != numSamples; j++) {
                        allSeqsAA[j] += allSeqs[j][k];
                    }
                }
            }
            
            // Only print unique haplotypes (and keep a list of other samples with the same):
            // otherSamples.resize(numSamples);
            for (std::vector<std::string>::size_type j = 0; j != numSamples; j++) {
                bool unique = true;
                for (std::vector<std::string>::size_type l = 0; l != uniqueHpt.size(); l++) {
                    if (allSeqsAA[j] == uniqueHpt[l]) {
                        unique = false;
                        if (otherSamples[l] != "") otherSamples[l] += ',';
                        otherSamples[l] += allNewLabels[j].substr(1);
                    }
                }
                if (unique) {
                    uniqueHpt.push_back(allSeqsAA[j]);
                    uniqueNames.push_back(allNewLabels[j]);
                    otherSamples.push_back("");
                }
            }
            
            string ancSeq = allSeqsAA[0];  // For now we don't know the ancestral states at the segregating sites
                                            // so just arbitrarily taking the first one as the reference
            for (string::size_type k = 0; k != ancSeq.length(); k++) {
                for (std::vector<std::string>::size_type j = 1; j != numSamples; j++) {
                    if (allSeqsAA[j][k] == ancSeq[k]) {
                        allSeqsAA[j][k] = '.';
                    }
                }
                for (std::vector<std::string>::size_type j = 1; j != uniqueHpt.size(); j++) {
                    if (uniqueHpt[j][k] == ancSeq[k]) {
                        uniqueHpt[j][k] = '.';
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
                AAFileName = opt::outFolder + stripExtension(allAligmentFiles[i]) + "_ssPos.txt";
            else
                AAFileName = opt::outFolder + allAligmentFiles[i] + "_ssPos.txt";
        } else {
            AAFileName = opt::outFolder + allAligmentFiles[i] + "_ssPos.txt";
        }
        //std::cerr << "into: " << AAFileName << std::endl;
        std::ofstream* AAalignment = new std::ofstream(AAFileName.c_str());
        *AAalignment << "Amino Acid positions:" << "\t"; print_vector(ssPos, *AAalignment,' ');
        /*
        for (std::vector<std::string>::size_type j = 0; j != numSamples; j++) {
            if (allNewLabels[j].length() < 16)
                *AAalignment << allNewLabels[j] << "\t\t" << allSeqsAA[j] << std::endl;
            else
                *AAalignment << allNewLabels[j] << "\t" << allSeqsAA[j] << std::endl;
        } AAalignment->close(); */
        for (std::vector<std::string>::size_type j = 0; j != uniqueHpt.size(); j++) {
            if (uniqueNames[j].length() < 16)
                *AAalignment << uniqueNames[j] << "\t\t" << uniqueHpt[j] << "\t" << otherSamples [j] << std::endl;
            else
                *AAalignment << uniqueNames[j] << "\t" << uniqueHpt[j] << "\t" << otherSamples [j] << std::endl;
        } AAalignment->close();
        
        
        
    }
    return 0;
}


void parseProteinSsOptions(int argc, char** argv) {
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
