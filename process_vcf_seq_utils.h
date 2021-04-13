//
//  process_vcf_seq_utils.h
//  process_vcf
//
//  Created by Milan Malinsky on 17/10/2013.
//  Copyright (c) 2013 Milan Malinsky. All rights reserved.
//

#ifndef process_vcf_process_vcf_seq_utils_h
#define process_vcf_process_vcf_seq_utils_h
#include <iostream>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include "process_vcf_IUPAC.h"
#include "process_vcf_utils.h"


inline void appendGenotypeBaseToString(std::string& toExtend, const std::string& ref, const std::string& alt, const std::vector<char>& genotype, char hetTreatment) {
    if (genotype[0] == '0' && genotype[1] == '0')
        toExtend.append(ref);
    else if (genotype[0] == '1' && genotype[1] == '1')
        toExtend.append(alt);
    else {
        if (hetTreatment == 'r') {
            double rn = ((double) rand() / RAND_MAX);
            if (rn <= 0.5) {
                toExtend.append(ref);
            } else {
                toExtend.append(alt);
            }
        } else if(hetTreatment == 'p') {
            if (genotype[0] == '0')
                toExtend.append(ref);
            if (genotype[0] == '1')
                toExtend.append(alt);
        } else if (hetTreatment == 'i') {
            std::string ambiguityBase = getAmbiguityCode(ref, alt);
            toExtend.append(ambiguityBase);
        } else if (hetTreatment == 'b') {
            if (genotype[1] == '0')
                toExtend.append(ref);
            if (genotype[1] == '1')
                toExtend.append(alt);
        }
    }
}

inline std::vector<std::string> returnGenotypeBaseAndZeroOne(const std::string& ref, const std::string& alt, const std::vector<char>& genotype, char hetTreatment) {
    std::vector<std::string> baseZeroOne;
    if (genotype[0] == '0' && genotype[1] == '0') {
        baseZeroOne.push_back(ref); baseZeroOne.push_back("0");
        return baseZeroOne;
    } else if (genotype[0] == '1' && genotype[1] == '1') {
        baseZeroOne.push_back(alt); baseZeroOne.push_back("1");
        return baseZeroOne;
    } else if (genotype[0] == '.' && genotype[1] == '.') { // Missing data
        baseZeroOne.push_back("."); baseZeroOne.push_back("0");
        return baseZeroOne;
    } else {
        if (hetTreatment == 'r') {
            double rn = ((double) rand() / RAND_MAX);
            if (rn <= 0.5) {
                baseZeroOne.push_back(ref); baseZeroOne.push_back("0");
                return baseZeroOne;
            } else {
                baseZeroOne.push_back(alt); baseZeroOne.push_back("1");
                return baseZeroOne;
            }
        } else if(hetTreatment == 'p') {
            if (genotype[0] == '0') {
                baseZeroOne.push_back(ref); baseZeroOne.push_back("0");
                return baseZeroOne;
            } else if (genotype[0] == '1') {
                baseZeroOne.push_back(alt); baseZeroOne.push_back("1");
                return baseZeroOne;
            }
        } else if (hetTreatment == 'b') {
            if (genotype[1] == '0') {
                baseZeroOne.push_back(ref); baseZeroOne.push_back("0");
                return baseZeroOne;
            } else if (genotype[1] == '1') {
                baseZeroOne.push_back(alt); baseZeroOne.push_back("1");
                return baseZeroOne;
            }
        } else {
            exit(1);
        }
    }
    exit(1);
}




inline std::string returnGenotypeBaseZeroOne(const std::vector<char>& genotype, char hetTreatment) {
    if (genotype[0] == '0' && genotype[1] == '0')
        return "0";
    else if (genotype[0] == '1' && genotype[1] == '1')
        return "1";
    else {
        if (hetTreatment == 'r') {
            double rn = ((double) rand() / RAND_MAX);
            if (rn <= 0.5) {
                return "0";
            } else {
                return "1";
            }
        } else if(hetTreatment == 'p') {
            if (genotype[0] == '0')
                return "0";
            if (genotype[0] == '1')
                return "1";
        } else if (hetTreatment == 'b') {
            if (genotype[1] == '0')
                return "0";
            if (genotype[1] == '1')
                return "1";
        } else {
            exit(1);
        }
    }
    exit(1);
}


// Read a scaffold from a reference genome fasta file into a single string +
// put the name of the next scaffold into the "nextScaffoldName" variable
inline std::string readScaffold(std::ifstream*& genomeFile, std::string& nextScaffoldName) {
    std::string scaffoldString;
    scaffoldString.reserve(100000000);
    std::string line = "";
    while (getline(*genomeFile, line) && line[0] != '>') {
        scaffoldString.append(line);
    }
    if (line[0] == '>')
        nextScaffoldName = split(line,' ')[0];
    else
        nextScaffoldName = "";
    return scaffoldString;
}

// Read the last scaffold from a reference genome fasta file into a single string
/*inline std::string readLastScaffold(std::ifstream*& genomeFile) {
    std::string scaffoldString;
    scaffoldString.reserve(100000000);
    std::string line = "";
    while (getline(*genomeFile, line) && line[0] != '>') {
        scaffoldString.append(line);
    }
    return scaffoldString;
} */

// Deal with the situation that there are scaffolds without any variable sites
// read through the genome file until the next scaffold to be read is the scaffold referred to in the vcf file
inline void forwardGenomeToScaffold(const std::string& thisInVCF, std::ifstream*& genomeFile, std::string& nextInGenome) {
    while (thisInVCF != nextInGenome) {
        std::cerr << "Starting to read " << nextInGenome << std::endl;
        std::cerr << "No variants in " << nextInGenome << std::endl;
        std::string currentScaffoldReference = readScaffold(genomeFile, nextInGenome);
        std::cerr << "Finished reading" << std::endl;
        nextInGenome.erase(0,1);
    }
}

#endif
