//
//  process_vcf_get_sequences.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 17/07/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include <iostream>
#include "process_vcf_utils.h"
#include "process_vcf_IUPAC.h"
#include "process_vcf_get_sequences.h"
#include "process_vcf_print_routines.h"

/* 
 TO DO:
 High priority:
 - deal with the possibility that a scaffold may not have any variants (DONE - apart from scaffold_0)
 
 */


#define SUBPROGRAM "getWGSeq"

#define DEBUG 1

static const char *GETSEQ_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] VCF_FILE GENOME_SEQUENCE.fa ANNOTATION.gffExtract\n"
"Obtain full genome sequences from a VCF file (e.g. for multiple alignment and phylogenetic analyses), output to STD_OUT\n"
"\n"
"       -h, --help                                  display this help and exit\n"
"       --by-scaffold                               output by scaffold/LG (each scaffold/LG) has its own file with sequences\n"
"                                                   for all samples\n"
"       --whole-genome                              output is one file with the whole genome concatenated for all samples\n"
"       --split NUM                                 split output into sequences containing approx. NUM segregating sites\n"
"                                                   each file contains sequences for all samples; this is intended for phylogenetic analyses\n"
"                                                   incompatible with --by-scaffold\n"
"       --LDhat                                     generate output sequences for the LDhat program (compatible with v2.1)\n"
"                                                   can be used in conjunction with --split\n"
"       --mtDNA                                     get only scaffold_747 and scaffold_2036 (corresponding to M.zebra mtDNA\n"
"       -s SAMPLES.txt, --samples=SAMPLES.txt       supply a file of sample identifiers to be used for the output\n"
"                                                   (default: sample ids from the vcf file are used)\n"
"       --incl-Pn=Mz_coords.PNsequence.NoIndels.fa  Also include a P.nyererei (outgroup) sequence (for now works only with --split)\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum { OPT_LDHAT, OPT_BY_SCAFFOLD, OPT_SPLIT, OPT_WG, OPT_PN };

static const char* shortopts = "hpws:c";

static const struct option longopts[] = {
    { "samples",   required_argument, NULL, 's' },
    { "by-scaffold",   no_argument, NULL, OPT_BY_SCAFFOLD },
    { "whole-genome",   no_argument, NULL, OPT_WG },
    { "LDhat",   no_argument, NULL, OPT_LDHAT },
    { "split",   required_argument, NULL, OPT_SPLIT },
    { "incl-Pn",   required_argument, NULL, OPT_PN },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string genomeFile;
    static string outgroupFile;
    static bool bLDhat = false;
    static bool bByScaffold = false;
    static bool bWholeGenome = false; // Print the whole genome into one file
    static string sampleNameFile;
    static int splitNum = 0;

}



int getSeqMain(int argc, char** argv) {
    string line;
    parseGetSeqOptions(argc, argv);
    string vcfFileRoot = stripExtension(opt::vcfFile);
    
    std::cerr << "Generating full genome sequences using variants from: " << opt::vcfFile << std::endl;
    std::cerr << "and the reference genome: " << opt::genomeFile << std::endl;
    if (opt::splitNum > 0)
        std::cerr << "with splits at every " << opt::splitNum << " variants" << std::endl;
    
    // Open connections to read from the vcf and reference genome files
    std::ifstream* vcfFile = new std::ifstream(opt::vcfFile.c_str());
    std::ifstream* genomeFile = new std::ifstream(opt::genomeFile.c_str());
    
    std::map<string, string> outgroupSeqs;
    if (!opt::outgroupFile.empty()) { outgroupSeqs = readMultiFastaToMap(opt::outgroupFile); }
    
    // Prepare an object for output files
    std::ofstream* wgFiles;
    
    string currentScaffoldNum = "";
    string currentScaffoldReference;
    string::size_type inStrPos;
    size_t numSamples;
    std::vector<string> sampleNames;
    string thisScaffoldName;
    std::vector<string> scaffoldStrings;
    unsigned int processedVariantCounter = 0;
    std::vector<string::size_type> splits;
    
    while (getline(*vcfFile, line)) {
        if (line[0] == '#' && line[1] == '#') 
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            std::vector<std::string> fields = split(line, '\t');
            numSamples = fields.size()-NUM_NON_GENOTYPE_COLUMNS;
            //std::cerr << "numSamples: " << numSamples << std::endl;
            // Initialize vectors
            scaffoldStrings.resize(numSamples);
            for (std::vector<string>::size_type i = 0; i != scaffoldStrings.size(); i++) {
                scaffoldStrings[i].reserve(30000000);
                scaffoldStrings[i] = "";
            }
            
            // Open output files
            if (opt::bWholeGenome) wgFiles = new std::ofstream[numSamples];
            for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                if (opt::sampleNameFile.empty())
                    sampleNames.push_back(fields[i]);
                else 
                    sampleNames = readSampleNamesFromTextFile(opt::sampleNameFile);
                if (opt::bWholeGenome) {
                    wgFiles[i-NUM_NON_GENOTYPE_COLUMNS].open(sampleNames[i-NUM_NON_GENOTYPE_COLUMNS].c_str());
                    wgFiles[i-NUM_NON_GENOTYPE_COLUMNS] << ">" << sampleNames[i-NUM_NON_GENOTYPE_COLUMNS] << std::endl;
                }
            }
        } else {
            processedVariantCounter++;
            std::vector<std::string> fields = split(line, '\t');
            std::vector<std::string> info = split(fields[7], ';');
            if (fields[0] != currentScaffoldNum) {
                if (currentScaffoldNum != "") {
                    for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                        scaffoldStrings[i].append(currentScaffoldReference.substr(inStrPos, string::npos));
                    }
                    
                    #ifdef DEBUG
                    if (scaffoldStrings[0].length() != currentScaffoldReference.length()) {
                        std::cerr << "Error!!! Reference scaffold/LG length: " << currentScaffoldReference.length() << " vcf scaffold length: " << scaffoldStrings[0].length() << std::endl;
                    }
                    #endif
                    
                    std::ofstream* scaffoldFile;
                    if (opt::splitNum == 0 && !opt::bWholeGenome) scaffoldFile = new std::ofstream(currentScaffoldNum.c_str());
                    
                    
                    std::cerr << currentScaffoldNum << " processed. Total variants: " << processedVariantCounter << " Writing output files..." << std::endl;
                    
                    if (opt::splitNum > 0) {
                        if (!opt::outgroupFile.empty()) {
                            print_split_incl_outgroup(currentScaffoldNum, splits, sampleNames, numSamples, scaffoldStrings, processedVariantCounter, outgroupSeqs, "Pnyererei");
                        } else {
                            print_split(currentScaffoldNum, splits, sampleNames, numSamples, scaffoldStrings, processedVariantCounter);
                        }
                    } else {
                        if (opt::bLDhat || opt::bByScaffold) {
                            if (opt::bLDhat)
                                *scaffoldFile << numSamples << "\t" << scaffoldStrings[0].length() << "\t" << "2" << std::endl;
                            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                                *scaffoldFile << ">" << sampleNames[i] << std::endl;
                                print80bpPerLineFile(scaffoldFile, scaffoldStrings[i]);
                                scaffoldStrings[i] = "";
                            }
                        } else if (opt::bWholeGenome) {
                            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                                print80bpPerLine(wgFiles, i, scaffoldStrings[i]);
                                scaffoldStrings[i] = "";

                            }
                        }
                    }
                    splits.clear();
                    processedVariantCounter = 1;
                    currentScaffoldNum = fields[0];
                    
                    while (currentScaffoldNum != thisScaffoldName) {
                        std::cerr << "Starting to read " << thisScaffoldName << std::endl;
                        std::cerr << "No variants in " << thisScaffoldName << std::endl;
                        std::cerr << "currentScaffoldNum " << currentScaffoldNum << std::endl;
                        currentScaffoldReference = readScaffold(genomeFile, thisScaffoldName);
                        if (opt::bWholeGenome) {
                            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                                print80bpPerLine(wgFiles, i, currentScaffoldReference);
                            }
                        }
                        std::cerr << "Finished reading" << std::endl;
                        thisScaffoldName.erase(0,1);
                        
                    }
                } else {
                    getline(*genomeFile, thisScaffoldName);
                    thisScaffoldName.erase(0,1);
                    currentScaffoldNum = fields[0];
                }
                
                // if (opt::bWholeGenome) printInAllOutputs(wgFiles, numSamples, thisScaffoldName);
                inStrPos = 0;
                             
                std::cerr << "Starting to read " << thisScaffoldName << std::endl;
                currentScaffoldReference = readScaffold(genomeFile, thisScaffoldName);
                thisScaffoldName.erase(0,1);
                std::cerr << "Finished reading" << std::endl;
                std::cerr << "Generating sequences with variants from the VCF file..." << std::endl;
            }
            if (info[0] != "INDEL") {
                for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                    //std::cerr << "Going through genotypes1:" << i << std::endl;
                    //std::cerr << scaffoldStrings.size() << " " << inStrPos << " " << fields[1] << " " << currentScaffoldReference.size() << std::endl;
                    std::vector<string> genotypeFields = split(fields[i], ':');
                    std::vector<string> genotype; genotype.push_back(numToString(genotypeFields[0][0])); genotype.push_back(numToString(genotypeFields[0][2]));
                    if (opt::bLDhat) {
                        for (int j = 0; j != ((atoi(fields[1].c_str()) - 1)-inStrPos); j++) {
                            scaffoldStrings[i- NUM_NON_GENOTYPE_COLUMNS].append("0");
                        }
                        if (genotype[0] == "0" && genotype[1] == "0")
                            scaffoldStrings[i- NUM_NON_GENOTYPE_COLUMNS].append("0");
                        else if (genotype[0] == "1" && genotype[1] == "1")
                            scaffoldStrings[i- NUM_NON_GENOTYPE_COLUMNS].append("1");
                        else {
                            string ambiguityBase = getAmbiguityCode(fields[3], fields[4]);
                            scaffoldStrings[i- NUM_NON_GENOTYPE_COLUMNS].append("2");
                        }
                    } else {
                        scaffoldStrings[i- NUM_NON_GENOTYPE_COLUMNS].append(currentScaffoldReference.substr(inStrPos, (atoi(fields[1].c_str()) - 1)-inStrPos));
                        appendGenotypeBaseToString(scaffoldStrings[i- NUM_NON_GENOTYPE_COLUMNS], fields[3], fields[4], genotype);
                    }
                }
                inStrPos = atoi(fields[1].c_str());

                #ifdef DEBUG
                if (currentScaffoldReference[inStrPos-1] != fields[3][0]) {
                    std::cerr << "Error!!! Sequence: " << currentScaffoldReference[inStrPos-1] << " vcf-ref: " << fields[3][0] << std::endl;
                }
                #endif
            }
            if (opt::splitNum > 0) {
                if (processedVariantCounter % opt::splitNum == 0) {
                    splits.push_back(inStrPos);
                    std::cerr << processedVariantCounter << " variants processed..." << std::endl;
                    std::cerr << "Split at bp: " << inStrPos << std::endl;
                }
            } else {
                if (processedVariantCounter % 10000 == 0)
                    std::cerr << processedVariantCounter << " variants processed..." << std::endl;
            }
        }
    }
    
    // Also the final scaffold
    for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
        scaffoldStrings[i].append(currentScaffoldReference.substr(inStrPos, string::npos));
    }
    std::ofstream* scaffoldFile;
    if (opt::splitNum == 0 && !opt::bWholeGenome) scaffoldFile = new std::ofstream(currentScaffoldNum.c_str());
    
    
    std::cerr << currentScaffoldNum << " processed. Total variants: " << processedVariantCounter << " Writing output files..." << std::endl;
    
    if (opt::splitNum > 0) {
        if (!opt::outgroupFile.empty()) {
            print_split_incl_outgroup(currentScaffoldNum, splits, sampleNames, numSamples, scaffoldStrings, processedVariantCounter, outgroupSeqs, "Pnyererei");
        } else {
            print_split(currentScaffoldNum, splits, sampleNames, numSamples, scaffoldStrings, processedVariantCounter);
        }
    } else {
        if (opt::bLDhat || opt::bByScaffold) {
            if (opt::bLDhat)
                *scaffoldFile << numSamples << "\t" << scaffoldStrings[0].length() << "\t" << "2" << std::endl;
            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                *scaffoldFile << ">" << sampleNames[i] << std::endl;
                print80bpPerLineFile(scaffoldFile, scaffoldStrings[i]);
                scaffoldStrings[i] = "";
            }
        } else if (opt::bWholeGenome) {
            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                print80bpPerLine(wgFiles, i, scaffoldStrings[i]);
                scaffoldStrings[i] = "";
                
            }
        }
    }
    
    
    /*
    if (opt::bWholeGenome) {
        if (system(NULL)) puts ("Ok");
        else exit (EXIT_FAILURE);
     
        string catAll = "cat ";
     
        for (std::vector<std::string>::size_type i = 0; i != sampleNames.size(); i++) {
            catAll += sampleNames[i] + " ";
        }
        catAll += "> " + vcfFileRoot + "whole_genome_all.fa";
        system (catAll.c_str());
            
    }*/
    
    return 0;
}

void print_split(const std::string& currentScaffoldNum, const std::vector<string::size_type>& splits, const std::vector<std::string>& sampleNames, const size_t numSamples, std::vector<std::string>& scaffoldStrings, const unsigned int totalProcessedVariants) {
    //for (std::vector<int>::size_type j = 0; j != splits.size(); j++) {
    
    std::cerr << "There are: " << splits.size() << " splits of size " << opt::splitNum << std::endl;
    
    std::ofstream* scaffoldFile;
    if (splits.size() == 0) {
        // Only print the whole scaffold if it contains a rasonble number of variants
        if (totalProcessedVariants % opt::splitNum > (opt::splitNum * 0.8)) {
            scaffoldFile = new std::ofstream(currentScaffoldNum.c_str());
            if (opt::bLDhat)
                *scaffoldFile << numSamples << "\t" << scaffoldStrings[0].length() << "\t" << "2" << std::endl;
            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                *scaffoldFile << ">" << sampleNames[i] << std::endl;
                print80bpPerLineFile(scaffoldFile, scaffoldStrings[i]);
                scaffoldStrings[i] = "";
            }
            scaffoldFile->close();
        } else {
            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                scaffoldStrings[i] = "";
            }
          
        }
    } else {
        string firstSplitFileName = currentScaffoldNum + "_1_" + numToString(splits[0]);
        scaffoldFile = new std::ofstream(firstSplitFileName.c_str());
        if (opt::bLDhat)
            *scaffoldFile << numSamples << "\t" << splits[0] << "\t" << "2" << std::endl;
        for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
            *scaffoldFile << ">" << sampleNames[i] << std::endl;
            print80bpPerLineFile(scaffoldFile, scaffoldStrings[i].substr(0,splits[0]));
        }
        scaffoldFile->close();
        
        for (std::vector<int>::size_type j = 1; j != splits.size(); j++) {
            string splitFileName = currentScaffoldNum + "_" + numToString(splits[j-1]+1) + "_" + numToString(splits[j]);
            scaffoldFile = new std::ofstream(splitFileName.c_str());
            if (opt::bLDhat)
                *scaffoldFile << numSamples << "\t" << (splits[j] - splits[j-1]) << "\t" << "2" << std::endl;
            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                *scaffoldFile << ">" << sampleNames[i] << std::endl;
                print80bpPerLineFile(scaffoldFile, scaffoldStrings[i].substr(splits[j-1], splits[j] - splits[j-1]));
            }
            scaffoldFile->close();
        }
        
        // Print the last part of the scaffold only if there are a reasonble number of variants
        if (totalProcessedVariants % opt::splitNum > (opt::splitNum * 0.8)) {
            string lastSplitFileName = currentScaffoldNum + "_" + numToString(splits[splits.size()-1]+1) + "_" + numToString(scaffoldStrings[0].length());
            
            scaffoldFile = new std::ofstream(lastSplitFileName.c_str());
            if (opt::bLDhat)
                *scaffoldFile << numSamples << "\t" << (scaffoldStrings[0].length() - splits[splits.size()-1]) << "\t" << "2" << std::endl;
            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                *scaffoldFile << ">" << sampleNames[i] << std::endl;
                print80bpPerLineFile(scaffoldFile, scaffoldStrings[i].substr(splits[splits.size()-1],std::string::npos));
                scaffoldStrings[i] = "";
            }
            scaffoldFile->close();
        } else {
            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                scaffoldStrings[i] = "";
            }
        }
    }
    
}

void print_split_incl_outgroup(const std::string& currentScaffoldNum, const std::vector<string::size_type>& splits, const std::vector<std::string>& sampleNames, const size_t numSamples, std::vector<std::string>& scaffoldStrings, const unsigned int totalProcessedVariants, std::map<string, string>& outgroupSeqs, const std::string& outgroupName) {
    //for (std::vector<int>::size_type j = 0; j != splits.size(); j++) {
    
    std::cerr << "There are: " << splits.size() << " splits of size " << opt::splitNum << std::endl;
    
    std::ofstream* scaffoldFile;
    if (splits.size() == 0) {
        // Only print the whole scaffold if it contains a rasonble number of variants
        if (totalProcessedVariants % opt::splitNum > (opt::splitNum * 0.8)) {
            scaffoldFile = new std::ofstream(currentScaffoldNum.c_str());
            if (opt::bLDhat)
                *scaffoldFile << numSamples << "\t" << scaffoldStrings[0].length() << "\t" << "2" << std::endl;
            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                *scaffoldFile << ">" << sampleNames[i] << std::endl;
                print80bpPerLineFile(scaffoldFile, scaffoldStrings[i]);
                scaffoldStrings[i] = "";
            }
            *scaffoldFile << ">" << outgroupName << std::endl;
            print80bpPerLineFile(scaffoldFile, outgroupSeqs[currentScaffoldNum]);
            scaffoldFile->close();
        } else {
            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                scaffoldStrings[i] = "";
            }
            
        }
    } else {
        string firstSplitFileName = currentScaffoldNum + "_1_" + numToString(splits[0]);
        scaffoldFile = new std::ofstream(firstSplitFileName.c_str());
        if (opt::bLDhat)
            *scaffoldFile << numSamples << "\t" << splits[0] << "\t" << "2" << std::endl;
        for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
            *scaffoldFile << ">" << sampleNames[i] << std::endl;
            print80bpPerLineFile(scaffoldFile, scaffoldStrings[i].substr(0,splits[0]));
        }
        *scaffoldFile << ">" << outgroupName << std::endl;
        print80bpPerLineFile(scaffoldFile, outgroupSeqs[currentScaffoldNum].substr(0,splits[0]));
        scaffoldFile->close();
        
        for (std::vector<int>::size_type j = 1; j != splits.size(); j++) {
            string splitFileName = currentScaffoldNum + "_" + numToString(splits[j-1]+1) + "_" + numToString(splits[j]);
            scaffoldFile = new std::ofstream(splitFileName.c_str());
            if (opt::bLDhat)
                *scaffoldFile << numSamples << "\t" << (splits[j] - splits[j-1]) << "\t" << "2" << std::endl;
            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                *scaffoldFile << ">" << sampleNames[i] << std::endl;
                print80bpPerLineFile(scaffoldFile, scaffoldStrings[i].substr(splits[j-1], splits[j] - splits[j-1]));
            }
            *scaffoldFile << ">" << outgroupName << std::endl;
            print80bpPerLineFile(scaffoldFile, outgroupSeqs[currentScaffoldNum].substr(splits[j-1], splits[j] - splits[j-1]));
            scaffoldFile->close();
        }
        
        // Print the last part of the scaffold only if there are a reasonble number of variants
        if (totalProcessedVariants % opt::splitNum > (opt::splitNum * 0.8)) {
            string lastSplitFileName = currentScaffoldNum + "_" + numToString(splits[splits.size()-1]+1) + "_" + numToString(scaffoldStrings[0].length());
            
            scaffoldFile = new std::ofstream(lastSplitFileName.c_str());
            if (opt::bLDhat)
                *scaffoldFile << numSamples << "\t" << (scaffoldStrings[0].length() - splits[splits.size()-1]) << "\t" << "2" << std::endl;
            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                *scaffoldFile << ">" << sampleNames[i] << std::endl;
                print80bpPerLineFile(scaffoldFile, scaffoldStrings[i].substr(splits[splits.size()-1],std::string::npos));
                scaffoldStrings[i] = "";
            }
            *scaffoldFile << ">" << outgroupName << std::endl;
            print80bpPerLineFile(scaffoldFile, outgroupSeqs[currentScaffoldNum].substr(splits[splits.size()-1],std::string::npos));
            scaffoldFile->close();
        } else {
            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                scaffoldStrings[i] = "";
            }
        }
    }
    
}





void printInAllOutputs(std::ofstream*& outFiles, size_t numSamples, string toPrint) {
    for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
        outFiles[i] << toPrint << std::endl;
    }
}

void print80bpPerLine(std::ofstream*& outFiles, std::vector<std::string>::size_type i, string toPrint) {
    string::size_type lines = toPrint.length() / 80;
    for (string::size_type j = 0; j <= lines; j++) {
        outFiles[i] << toPrint.substr(j*80,80) << std::endl;
    }
}

void parseGetSeqOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case '?': die = true; break;
            case 's': arg >> opt::sampleNameFile; break;
            case OPT_LDHAT: opt::bLDhat = true; break;
            case OPT_BY_SCAFFOLD: opt::bByScaffold = true; break;
            case OPT_SPLIT: arg >> opt::splitNum; break;
            case OPT_WG: opt::bWholeGenome = true; break;
            case OPT_PN: arg >> opt::outgroupFile; break;
            case 'h':
                std::cout << GETSEQ_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    if (argc - optind < 2) {
        std::cerr << "missing arguments\n";
        die = true;
    } 
    else if (argc - optind > 2)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (opt::splitNum > 0 && opt::bByScaffold) {
        std::cerr << "The option --split is incompatible with --by-scaffold\n";
        die = true;
    }
    
    if (opt::bWholeGenome && (opt::splitNum > 0 || opt::bByScaffold)) {
        std::cerr << "The option --whole-genome is incompatible with --by-scaffold and --split\n";
        die = true;
    }
    
    if (opt::splitNum < 0) {
        std::cerr << "The argument to --split cannot be negative\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << GETSEQ_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::genomeFile = argv[optind++];
}

