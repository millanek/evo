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
#include "process_vcf_annotation_tools.h"

/* 
 TO DO:
 High priority:
 - deal with the possibility that a scaffold may not have any variants (DONE - apart from scaffold_0)
 
 */


#define SUBPROGRAM "getWGSeq"

#define DEBUG 1

static const char *GETSEQ_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] VCF_FILE GENOME_SEQUENCE.fa\n"
"Obtain full genome sequences from a VCF file (e.g. for multiple alignment and phylogenetic analyses), output to STD_OUT\n"
"\n"
"       -h, --help                                  display this help and exit\n"
"       --by-scaffold                               output by scaffold/LG (each scaffold/LG) has its own file with sequences\n"
"                                                   for all samples\n"
"       --whole-genome                              output is one file with the whole genome concatenated for all samples\n"
"       --methylome                                 This is for Greg's study - C->T and G->A and VCF with .fa can be revesrse strands\n"
"       -H,   --het-treatment <r|p|b|i>             r: assign het bases randomly (default); p: use the phase information in a VCF outputting\n"
"                                                   haplotype 1 for each individual; b: use both haplotypes as phased; i: use IUPAC codes\n"
"       --split NUM                                 split output into sequences containing approx. NUM segregating sites\n"
"                                                   each file contains sequences for all samples; this is intended for phylogenetic analyses\n"
"                                                   incompatible with --by-scaffold\n"
"       --LDhat                                     generate output sequences for the LDhat program (compatible with v2.1)\n"
"                                                   can be used in conjunction with --split\n"
"       --mtDNA                                     get only scaffold_747 and scaffold_2036 (corresponding to M.zebra mtDNA\n"
"       -s SAMPLES.txt, --samples=SAMPLES.txt       supply a file of sample identifiers to be used for the output\n"
"                                                   (default: sample ids from the vcf file are used)\n"
"       --incl-Pn=Mz_coords.PNsequence.NoIndels.fa  Also include an outgroup sequence (for now works only with --split)\n"
"       --accessibleGenomeBED=BEDfile.bed           (optional) a bed file specifying the regions of the genome where we could call SNPs\n"
"       --makeSVDinput                              Generates input for SVDquartets to STD_OUT\n"
"       --makeBootstrapSeqs=FILENAME_ROOT           Generate bootstrap sequence replicates for SVDquartets\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum { OPT_LDHAT, OPT_BY_SCAFFOLD, OPT_SPLIT, OPT_WG, OPT_PN, OPT_ACC_GEN_BED, OPT_SVD, OPT_SVD_BOOT, OPT_METH };

static const char* shortopts = "hs:H:";

static const struct option longopts[] = {
    { "samples",   required_argument, NULL, 's' },
    { "by-scaffold",   no_argument, NULL, OPT_BY_SCAFFOLD },
    { "het-treatment",   required_argument, NULL, 'H' },
    { "whole-genome",   no_argument, NULL, OPT_WG },
    { "LDhat",   no_argument, NULL, OPT_LDHAT },
    { "split",   required_argument, NULL, OPT_SPLIT },
    { "incl-Pn",   required_argument, NULL, OPT_PN },
    { "accessibleGenomeBED", required_argument, NULL, OPT_ACC_GEN_BED },
    { "makeSVDinput", no_argument, NULL, OPT_SVD },
    { "methylome", no_argument, NULL, OPT_METH },
    { "makeBootstrapSeqs", required_argument, NULL, OPT_SVD_BOOT },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string genomeFile = "";
    static string outgroupFile;
    static bool bLDhat = false;
    static bool bByScaffold = false;
    static bool bWholeGenome = false; // Print the whole genome into one file
    static string sampleNameFile;
    static int splitNum = 0;
    static char hetTreatment = 'r';
    static bool bSVD = false;
    static string bootSVDnameRoot;
    static string accesibleGenBedFile;
    static bool methylome = false;

}



int getSeqMain(int argc, char** argv) {
    string line;
    parseGetSeqOptions(argc, argv);
    string vcfFileRoot = stripExtension(opt::vcfFile);
    
    std::cerr << "Generating full genome sequences using variants from: " << opt::vcfFile << std::endl;
    std::cerr << "and the reference genome: " << opt::genomeFile << std::endl;
    if (opt::splitNum > 0)
        std::cerr << "with splits at every " << opt::splitNum << " variants" << std::endl;
    
    //std::cerr << "Bootstrap sequences will be output to " << opt::bootSVDnameRoot << "_i_boot.txt" << std::endl;
    if (!opt::bootSVDnameRoot.empty())
        std::cerr << "Bootstrap sequences will be output to " << opt::bootSVDnameRoot << "_i_boot.txt" << std::endl;
    
    // Open connections to read from the vcf and reference genome files
    std::istream* vcfFile = createReader(opt::vcfFile.c_str());
    std::ifstream* genomeFile = new std::ifstream(opt::genomeFile.c_str());
    std::ifstream* accessibleGenomeBed;
    
    std::map<string, string> outgroupSeqs;
    std::map<string, int> fullScaffoldLengths;
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
    unsigned int usedVariantCounter = 0;
    std::vector<string::size_type> splits;
    
    AccessibleGenome* ag;
    if (!opt::accesibleGenBedFile.empty()) {
        accessibleGenomeBed = new std::ifstream(opt::accesibleGenBedFile);
        std::cerr << "Loading the accessible genome annotation" << std::endl;
        ag = new AccessibleGenome(accessibleGenomeBed);
        std::cerr << "Done" << std::endl;
    }
    
    
    while (getline(*vcfFile, line)) {
        if (line[0] == '#' && line[1] == '#') 
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            std::vector<std::string> fields = split(line, '\t');
            numSamples = fields.size()-NUM_NON_GENOTYPE_COLUMNS;
            std::cerr << "numSamples: " << numSamples << std::endl;
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
                if (!opt::bSVD) {
                    if (opt::bWholeGenome) {
                        wgFiles[i-NUM_NON_GENOTYPE_COLUMNS].open(sampleNames[i-NUM_NON_GENOTYPE_COLUMNS].c_str());
                        //wgFiles[i-NUM_NON_GENOTYPE_COLUMNS] << ">" << sampleNames[i-NUM_NON_GENOTYPE_COLUMNS] << std::endl;
                    }
                }
            }
        } else {
            processedVariantCounter++;
            std::vector<std::string> fields = split(line, '\t');
            std::vector<std::string> info = split(fields[7], ';');
            if (fields[0] != currentScaffoldNum) {
                if (currentScaffoldNum != "") {
                    if (opt::genomeFile != "") {
                        for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                            scaffoldStrings[i].append(currentScaffoldReference.substr(inStrPos, string::npos));
                        }
    #ifdef DEBUG
                        if (scaffoldStrings[0].length() != currentScaffoldReference.length()) {
                            std::cerr << "Error!!! Reference scaffold/LG length: " << currentScaffoldReference.length() << " vcf scaffold length: " << scaffoldStrings[0].length() << std::endl;
                        }
#endif
                    }
                    std::ofstream* scaffoldFile;
                    if (opt::splitNum == 0 && !opt::bWholeGenome) scaffoldFile = new std::ofstream(currentScaffoldNum.c_str());
                    
                    std::cerr << currentScaffoldNum << " processed. Total variants: " << processedVariantCounter << std::endl;
                    if (!opt::bSVD) {
                        std::cerr << " Writing output files..." << std::endl;
                    }
                    fullScaffoldLengths[currentScaffoldNum] = (int)scaffoldStrings[0].length();
                    
                    // if requested, reduce the strings to accessible sequence only - this should still include the "right" number of variants because variants can appear only in the accessible sequence
                    if (!opt::accesibleGenBedFile.empty()) {
                        std::cerr << "Reducing scaffoldStrings to accesible genome only.." << scaffoldStrings[0].length() << std::endl;
                        for (int i = 0; i < scaffoldStrings.size(); i++) {
                            scaffoldStrings[i] = ag->getAccessibleSeqForScaffold(currentScaffoldNum,scaffoldStrings[i]);
                        }
                        std::cerr << "after reduction -> scaffoldStrings[0] length: " << scaffoldStrings[0].length() << std::endl;
                        // also needs to be done for the outgroup, if present
                        if (!opt::outgroupFile.empty()) {
                            outgroupSeqs[currentScaffoldNum] = ag->getAccessibleSeqForScaffold(currentScaffoldNum,outgroupSeqs[currentScaffoldNum]);
                        }
                    }
                    if (opt::splitNum > 0) {
                        std::vector<string::size_type> scaledSplits = splits;
                        std::cerr << "Splits" << std::endl;
                        print_vector(splits, std::cerr);
                        if (!opt::accesibleGenBedFile.empty()) { // Need to rescale the splits
                            for (int i = 0; i < splits.size(); i++) {
                                scaledSplits[i] = ag->getAccessibleBPinRegion(currentScaffoldNum, 0, (int)splits[i]);
                            }
                            std::cerr << "Scaled splits:" << std::endl;
                            print_vector(scaledSplits, std::cerr);
                        }
                        
                        if (!opt::outgroupFile.empty()) {
                            print_split_incl_outgroup(currentScaffoldNum, splits, sampleNames, numSamples, scaffoldStrings, processedVariantCounter, outgroupSeqs, "Outgroup",scaledSplits,fullScaffoldLengths[currentScaffoldNum]);
                        } else {
                            print_split(currentScaffoldNum, splits, sampleNames, numSamples, scaffoldStrings, processedVariantCounter,scaledSplits,fullScaffoldLengths[currentScaffoldNum]);
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
                                if (!opt::bSVD) {
                                    print80bpPerLine(wgFiles, i, scaffoldStrings[i]);
                                    scaffoldStrings[i] = "";
                                }
                            }
                        }
                    }
                    splits.clear();
                    if (!opt::bSVD) {
                        processedVariantCounter = 1;
                    }
                    currentScaffoldNum = fields[0];
                    
                    if (opt::genomeFile != "") {
                        while (currentScaffoldNum != thisScaffoldName) {
                            std::cerr << "Starting to read " << thisScaffoldName << std::endl;
                            std::cerr << "No variants in " << thisScaffoldName << std::endl;
                            std::cerr << "currentScaffoldNum " << currentScaffoldNum << std::endl;
                            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                                if (!opt::bSVD) {
                                    wgFiles[i] << ">" << thisScaffoldName << std::endl;
                                }
                            }
                            currentScaffoldReference = readScaffold(genomeFile, thisScaffoldName);
                            std::cerr << "Finished reading" << std::endl;
                            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                                if (!opt::bSVD) {
                                    print80bpPerLine(wgFiles, i, currentScaffoldReference);
                                }
                            }
                            thisScaffoldName.erase(0,1);
                            
                        }
                    }
                } else {
                    std::cerr << "currentScaffoldNum: " << currentScaffoldNum << std::endl;
                    getline(*genomeFile, thisScaffoldName);
                    thisScaffoldName.erase(0,1);
                    currentScaffoldNum = fields[0];
                }
                
                // if (opt::bWholeGenome) printInAllOutputs(wgFiles, numSamples, thisScaffoldName);
                inStrPos = 0;
                             
                std::cerr << "Starting to read " << thisScaffoldName << std::endl;
                for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                    if (!opt::bSVD) {
                        wgFiles[i] << ">" << thisScaffoldName << std::endl;
                    }
                }
                if (opt::genomeFile != "") {
                    currentScaffoldReference = readScaffold(genomeFile, thisScaffoldName);
                }
                thisScaffoldName.erase(0,1);
                std::cerr << "Finished reading" << std::endl;
                std::cerr << "Generating sequences with variants from the VCF file..." << std::endl;
            }
            if (info[0] != "INDEL") {
                int lengthToAppend = (atoi(fields[1].c_str()) - 1) - (int)inStrPos;
                // make sure the length is non-negative (can happen
                // if two consecutive variants have the same coordinate)
                // for now we just ignore the additional variant
                if (lengthToAppend >= 0) {
                    std::vector<int> appendVectorInt(numSamples,0);
                    std::vector<std::string> appendVector(numSamples,"0");
                    for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                        //std::cerr << "Going through genotypes1:" << i << std::endl;
                        //std::cerr << scaffoldStrings.size() << " " << inStrPos << " " << fields[1] << " " << currentScaffoldReference.size() << std::endl;
                        std::vector<string> genotypeFields = split(fields[i], ':');
                        std::vector<char> genotype; genotype.push_back(genotypeFields[0][0]); genotype.push_back(genotypeFields[0][2]);
                        if (opt::bLDhat) {
                            for (int j = 0; j != ((atoi(fields[1].c_str()) - 1)-inStrPos); j++) {
                                scaffoldStrings[i- NUM_NON_GENOTYPE_COLUMNS].append("0");
                            }
                            if (genotype[0] == '0' && genotype[1] == '0')
                                scaffoldStrings[i- NUM_NON_GENOTYPE_COLUMNS].append("0");
                            else if (genotype[0] == '1' && genotype[1] == '1')
                                scaffoldStrings[i- NUM_NON_GENOTYPE_COLUMNS].append("1");
                            else {
                                string ambiguityBase = getAmbiguityCode(fields[3], fields[4]);
                                scaffoldStrings[i- NUM_NON_GENOTYPE_COLUMNS].append("2");
                            }
                        } else {
                            //std::cerr << "lengthToAppend: " << lengthToAppend << std::endl;
                            //std::cerr << "opt::genomeFile: " << opt::genomeFile << std::endl;
                            if (opt::genomeFile != "") {
                                scaffoldStrings[i- NUM_NON_GENOTYPE_COLUMNS].append(currentScaffoldReference.substr(inStrPos, lengthToAppend));
                                
                            }
                            if (opt::bSVD) {
                                std::vector<std::string> genotypeAndZeroOne = returnGenotypeBaseAndZeroOne(fields[3], fields[4], genotype, opt::hetTreatment);
                                appendVector[i- NUM_NON_GENOTYPE_COLUMNS] = genotypeAndZeroOne[0];
                                appendVectorInt[i- NUM_NON_GENOTYPE_COLUMNS] = (int)stringToDouble(genotypeAndZeroOne[1].c_str());
                            } else {
                                if (opt::methylome) {
                                    std::cerr << "Here: " << std::endl;
                                    char currentFastaBase = currentScaffoldReference[atoi(fields[1].c_str())-1];
                                    std::string currentFastaBaseStr(1, currentFastaBase);
                                    string VCFref = fields[3];
                                    string VCFalt = fields[4];
                                    if (currentFastaBaseStr == "C" && currentFastaBaseStr == "G") {
                                        fields[3][0] = complementIUPAC(VCFref[0]);
                                        fields[4][0] = complementIUPAC(VCFalt[0]);
                                    } else if (currentFastaBaseStr == "G" && VCFref == "C") {
                                        fields[3][0] = complementIUPAC(VCFref[0]);
                                        fields[4][0] = complementIUPAC(VCFalt[0]);
                                    }
                                }
                                appendGenotypeBaseToString(scaffoldStrings[i- NUM_NON_GENOTYPE_COLUMNS], fields[3], fields[4], genotype, opt::hetTreatment);
                            }
                        }
                    }
                    if (opt::bSVD) {
                        if(vector_sum(appendVectorInt) > 0) {
                            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                                scaffoldStrings[i].append(appendVector[i]);
                            }
                            usedVariantCounter++;
                        }
                    }
                }
                if (opt::bSVD) {
                    if((int)scaffoldStrings[0].length() != usedVariantCounter) {
                        std::cerr << "usedVariantCounter: " << usedVariantCounter << std::endl;
                        std::cerr << "scaffoldStrings[0].length(): " << scaffoldStrings[0].length() << std::endl;
                        std::cerr << "scaffoldStrings[1].length(): " << scaffoldStrings[1].length() << std::endl;
                        std::cerr << fields[3] << " " << fields[4] << std::endl;
                        std::cerr << "scaffoldStrings[0]: " << scaffoldStrings[0] << std::endl;
                        std::cerr << "scaffoldStrings[1]: " << scaffoldStrings[1] << std::endl;
                    }
                    assert((int)scaffoldStrings[0].length() == usedVariantCounter);
                }
                inStrPos = atoi(fields[1].c_str());

        #ifdef DEBUG
                if (opt::genomeFile != "") {
                if (currentScaffoldReference[inStrPos-1] != fields[3][0]) {
                //    std::cerr << "Error!!! Sequence: " << currentScaffoldReference[inStrPos-1] << " vcf-ref: " << fields[3][0] << std::endl;
                }
                }
        #endif
            }
            if (opt::splitNum > 0) {
                if (processedVariantCounter % opt::splitNum == 0) {
                    splits.push_back(inStrPos);
                    std::cerr << processedVariantCounter << " variants processed..." << std::endl;
                    std::cerr << "Split at bp: " << inStrPos << std::endl;
                    std::cerr << "scaffoldStrings[0].length(): " << scaffoldStrings[0].length() << std::endl;
                }
            } else {
                if (processedVariantCounter % 10000 == 0)
                    std::cerr << processedVariantCounter << " variants processed..." << std::endl;
            }
        }
    }
    
    // Also the final scaffold
    if (opt::genomeFile != "") {
        for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
            scaffoldStrings[i].append(currentScaffoldReference.substr(inStrPos, string::npos));
        }
    }
    std::ofstream* scaffoldFile;
    if (opt::splitNum == 0 && !opt::bWholeGenome) scaffoldFile = new std::ofstream(currentScaffoldNum.c_str());
    
    
    std::cerr << currentScaffoldNum << " processed. Total variants: " << processedVariantCounter << " Writing output files..." << std::endl;
    fullScaffoldLengths[currentScaffoldNum] = (int)scaffoldStrings[0].length();
    
    // if requested, reduce the strings to accessible sequence only - this should still include the "right" number of variants because variants can appear only in the accessible sequence
    if (!opt::accesibleGenBedFile.empty()) {
        std::cerr << "Reducing scaffoldStrings to accesible genome only.." << scaffoldStrings[0].length() << std::endl;
        for (int i = 0; i < scaffoldStrings.size(); i++) {
            scaffoldStrings[i] = ag->getAccessibleSeqForScaffold(currentScaffoldNum,scaffoldStrings[i]);
        }
        std::cerr << "after reduction -> scaffoldStrings[0] length: " << scaffoldStrings[0].length() << std::endl;
        // also needs to be done for the outgroup, if present
        if (!opt::outgroupFile.empty()) {
            outgroupSeqs[currentScaffoldNum] = ag->getAccessibleSeqForScaffold(currentScaffoldNum,outgroupSeqs[currentScaffoldNum]);
        }
    }
    
    
    if (opt::splitNum > 0) {
        std::vector<string::size_type> scaledSplits = splits;
        if (!opt::accesibleGenBedFile.empty()) { // Need to rescale the splits
            for (int i = 0; i < splits.size(); i++) {
                scaledSplits[i] = ag->getAccessibleBPinRegion(currentScaffoldNum, 0, (int)splits[i]);
            }
        }
        if (!opt::outgroupFile.empty()) {
            print_split_incl_outgroup(currentScaffoldNum, splits, sampleNames, numSamples, scaffoldStrings, processedVariantCounter, outgroupSeqs, "Outgroup",scaledSplits,fullScaffoldLengths[currentScaffoldNum]);
        } else {
            print_split(currentScaffoldNum, splits, sampleNames, numSamples, scaffoldStrings, processedVariantCounter,scaledSplits,fullScaffoldLengths[currentScaffoldNum]);
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
            if (opt::bSVD) {
                std::cout << "#NEXUS" << std::endl;
                std::cout << "begin data;" << std::endl;
                std::cout << "dimensions ntax=" << numSamples << " nchar=" << scaffoldStrings[0].length() << ";" << std::endl;
                std::cout << "format datatype=dna missing=." << ";" << std::endl;
                std::cout << "matrix" << std::endl;
            }
            std::vector<string> editedSnVector;
            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                
                string editedSn = sampleNames[i];
                int snLength = (int)sampleNames[i].length();
                if (snLength < 32) {
                    for (int p = snLength; p <= 32; p++) {
                        editedSn = editedSn + " ";
                    }
                }
                editedSnVector.push_back(editedSn.substr(0,32));
                if (!opt::bSVD) {
                    print80bpPerLine(wgFiles, i, scaffoldStrings[i]);
                } else {
                    std::cout << editedSnVector[i] << "\t" << scaffoldStrings[i] << std::endl;
                }
            }
            if (opt::bSVD) {
                std::cout << ";" << std::endl;
                std::cout << "end;" << std::endl;
            }
                
            if (!opt::bootSVDnameRoot.empty()) { // Output bootstrap sequences
                int totalLength = (int)scaffoldStrings[0].length();
                std::cerr << totalLength << "=totalLength;" << std::endl;
                std::random_device rd; // obtain a random number from hardware
                std::mt19937 eng(rd()); // seed the generator
                std::uniform_int_distribution<> randomPos(0, totalLength-1); // define the range
                for (std::vector<std::string>::size_type i = 0; i != 100; i++) {
                    std::string bootFileName = opt::bootSVDnameRoot + "_" + numToString(i) + + "_boot.txt";
                    std::ofstream* bootFile = new std::ofstream(bootFileName.c_str());
                    int numSamples = (int)sampleNames.size();
                    std::vector<string> thisSeqs(numSamples,"");
                    for (int j = 0; j < totalLength; j++) {
                        int pos = randomPos(eng);
                        for (int k = 0; k < numSamples; k++) {
                         //   std::cerr << pos << "=pos;" << std::endl;
                            thisSeqs[k] += scaffoldStrings[k].substr(pos,1);
                        }
                    }
                    *bootFile << "#NEXUS" << std::endl;
                    *bootFile << "begin data;" << std::endl;
                    *bootFile << "dimensions ntax=" << numSamples << " nchar=" << thisSeqs[0].length() << ";" << std::endl;
                    *bootFile << "format datatype=dna missing=." << ";" << std::endl;
                    *bootFile << "matrix" << std::endl;
                    for (int k = 0; k < numSamples; k++) {
                        *bootFile << editedSnVector[k] << "\t" << thisSeqs[k] << std::endl;
                        thisSeqs[k] = "";
                    }
                    *bootFile << ";" << std::endl;
                    *bootFile << "end;" << std::endl;
                    bootFile->close();
                }
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

void print_split(const std::string& currentScaffoldNum, const std::vector<string::size_type>& splits, const std::vector<std::string>& sampleNames, const size_t numSamples, std::vector<std::string>& scaffoldStrings, const unsigned int totalProcessedVariants,const std::vector<string::size_type>& scaledSplits, const int fullScLength) {
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
            print80bpPerLineFile(scaffoldFile, scaffoldStrings[i].substr(0,scaledSplits[0]));
        }
        scaffoldFile->close();
        
        for (std::vector<int>::size_type j = 1; j != splits.size(); j++) {
            string splitFileName = currentScaffoldNum + "_" + numToString(splits[j-1]+1) + "_" + numToString(splits[j]);
            scaffoldFile = new std::ofstream(splitFileName.c_str());
            if (opt::bLDhat)
                *scaffoldFile << numSamples << "\t" << (splits[j] - splits[j-1]) << "\t" << "2" << std::endl;
            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                *scaffoldFile << ">" << sampleNames[i] << std::endl;
                print80bpPerLineFile(scaffoldFile, scaffoldStrings[i].substr(scaledSplits[j-1], scaledSplits[j] - scaledSplits[j-1]));
            }
            scaffoldFile->close();
        }
        
        // Print the last part of the scaffold only if there are a reasonble number of variants
        if (totalProcessedVariants % opt::splitNum > (opt::splitNum * 0.8)) {
            string lastSplitFileName = currentScaffoldNum + "_" + numToString(splits[splits.size()-1]+1) + "_" + numToString(fullScLength);
            
            scaffoldFile = new std::ofstream(lastSplitFileName.c_str());
            if (opt::bLDhat)
                *scaffoldFile << numSamples << "\t" << (scaffoldStrings[0].length() - splits[splits.size()-1]) << "\t" << "2" << std::endl;
            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                *scaffoldFile << ">" << sampleNames[i] << std::endl;
                print80bpPerLineFile(scaffoldFile, scaffoldStrings[i].substr(scaledSplits[scaledSplits.size()-1],std::string::npos));
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

void print_split_incl_outgroup(const std::string& currentScaffoldNum, const std::vector<string::size_type>& splits, const std::vector<std::string>& sampleNames, const size_t numSamples, std::vector<std::string>& scaffoldStrings, const unsigned int totalProcessedVariants, std::map<string, string>& outgroupSeqs, const std::string& outgroupName,const std::vector<string::size_type>& scaledSplits, const int fullScLength) {
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
            print80bpPerLineFile(scaffoldFile, scaffoldStrings[i].substr(0,scaledSplits[0]));
        }
        *scaffoldFile << ">" << outgroupName << std::endl;
        print80bpPerLineFile(scaffoldFile, outgroupSeqs[currentScaffoldNum].substr(0,scaledSplits[0]));
        scaffoldFile->close();
        
        for (std::vector<int>::size_type j = 1; j != splits.size(); j++) {
            string splitFileName = currentScaffoldNum + "_" + numToString(splits[j-1]+1) + "_" + numToString(splits[j]);
            scaffoldFile = new std::ofstream(splitFileName.c_str());
            if (opt::bLDhat)
                *scaffoldFile << numSamples << "\t" << (splits[j] - splits[j-1]) << "\t" << "2" << std::endl;
            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                *scaffoldFile << ">" << sampleNames[i] << std::endl;
                print80bpPerLineFile(scaffoldFile, scaffoldStrings[i].substr(scaledSplits[j-1], scaledSplits[j] - scaledSplits[j-1]));
            }
            *scaffoldFile << ">" << outgroupName << std::endl;
            print80bpPerLineFile(scaffoldFile, outgroupSeqs[currentScaffoldNum].substr(scaledSplits[j-1], scaledSplits[j] - scaledSplits[j-1]));
            scaffoldFile->close();
        }
        
        // Print the last part of the scaffold only if there are a reasonble number of variants
        if (totalProcessedVariants % opt::splitNum > (opt::splitNum * 0.8)) {
            string lastSplitFileName = currentScaffoldNum + "_" + numToString(splits[splits.size()-1]+1) + "_" + numToString(fullScLength);
            
            scaffoldFile = new std::ofstream(lastSplitFileName.c_str());
            if (opt::bLDhat)
                *scaffoldFile << numSamples << "\t" << (scaffoldStrings[0].length() - splits[splits.size()-1]) << "\t" << "2" << std::endl;
            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                *scaffoldFile << ">" << sampleNames[i] << std::endl;
                print80bpPerLineFile(scaffoldFile, scaffoldStrings[i].substr(scaledSplits[splits.size()-1],std::string::npos));
                scaffoldStrings[i] = "";
            }
            *scaffoldFile << ">" << outgroupName << std::endl;
            print80bpPerLineFile(scaffoldFile, outgroupSeqs[currentScaffoldNum].substr(scaledSplits[scaledSplits.size()-1],std::string::npos));
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
            case 'H': arg >> opt::hetTreatment; break;
            case OPT_LDHAT: opt::bLDhat = true; break;
            case OPT_BY_SCAFFOLD: opt::bByScaffold = true; break;
            case OPT_SPLIT: arg >> opt::splitNum; break;
            case OPT_WG: opt::bWholeGenome = true; break;
            case OPT_PN: arg >> opt::outgroupFile; break;
            case OPT_ACC_GEN_BED: arg >> opt::accesibleGenBedFile; break;
            case OPT_SVD: opt::bSVD = true; break;
            case OPT_SVD_BOOT: arg >> opt::bootSVDnameRoot; break;
            case OPT_METH: opt::methylome = true; break;
            case 'h':
                std::cout << GETSEQ_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    if (argc - optind < 1) {
        std::cerr << "missing arguments\n";
        die = true;
    } 
    else if (argc - optind > 2)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (opt::hetTreatment != 'r' && opt::hetTreatment != 'p' && opt::hetTreatment != 'b' && opt::hetTreatment != 'i') {
        std::cerr << "The -H (--het-treatment) option can only have the values 'r', 'p', 'b', or 'i'\n";
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
    //std::cerr << "argc - optind " << argc - optind << std::endl;
    if( argc - optind == 1) {
        opt::genomeFile = argv[optind++];
    }
}

