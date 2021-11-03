//
//  process_vcf_utils.cpp
//  C++ project
//
//  Created by Milan Malinsky on 08/05/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include <iostream>
#include "process_vcf_utils.h"
#include "process_vcf_stats_utils.h"



void printMissingLikelihoodsWarning(const string& chr, const string& pos) {
    std::cerr << "WARNING: Could not fing genotype likelihoods/probabilities (GP, PL, or GL fields) for variant at " << chr << " " << pos << std::endl;
    std::cerr << "WARNING: Did you really mean to use the -g option? Reverting to using called genotypes." << std::endl;
}


// Works only on biallelic markers
void GeneralSetCounts::getSetVariantCounts(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap) {
    
    getBasicCounts(genotypes, posToSpeciesMap);
    
    // If at least one of the outgroup individuals has non-missing data
    // Find out what is the "ancestral allele" - i.e. the one more common in the outgroup
    int AAint;
    try {
        if (setAlleleCounts.at("Outgroup") > 0) {
            if (setRefCounts.at("Outgroup") > setAltCounts.at("Outgroup")) { AAint = 0; }
            else { AAint = 1; }
        }
    } catch (std::out_of_range& e) { AAint = -1; }
    
    // Now fill in the allele frequencies
    for(std::map<string,int>::iterator it = setAltCounts.begin(); it != setAltCounts.end(); ++it) {
        if (setAlleleCounts.at(it->first) > 0) {
            setAAFs[it->first] = (double)setAltCounts.at(it->first)/setAlleleCounts.at(it->first);
            if (AAint == 0) { // Ancestral allele seems to be the ref, so derived is alt
                setDAFs[it->first] = (double)setAltCounts.at(it->first)/setAlleleCounts.at(it->first);
            } else if (AAint == 1) { // Ancestral allele seems to be alt, so derived is ref
                setDAFs[it->first] = (double)setRefCounts.at(it->first)/setAlleleCounts.at(it->first);
            }
        }
    }
}

// Works only on biallelic markers
void GeneralSetCounts::getSetVariantCountsSimple(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap) {
    // std::cerr << fields[0] << "\t" << fields[1] << std::endl;
    getBasicCounts(genotypes, posToSpeciesMap);
    
    // Now fill in the allele frequencies
    for(std::map<string,int>::iterator it = setAltCounts.begin(); it != setAltCounts.end(); ++it) {
        if (setAlleleCounts.at(it->first) > 0) {
            setAAFs[it->first] = (double)setAltCounts.at(it->first)/setAlleleCounts.at(it->first);
        }
    }
}

void GeneralSetCounts::getBasicCounts(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap) {
    // Go through the genotypes - only biallelic markers are allowed
    for (std::vector<std::string>::size_type i = 0; i != genotypes.size(); i++) {
        std::string species = posToSpeciesMap.at(i);
        // The first allele in this individual
        if (genotypes[i][0] == '1') {
            overall++; individualsWithVariant[i] = 1;
            setAltCounts[species]++; setAlleleCounts[species]++;
        } else if (genotypes[i][0] == '0') {
            setAlleleCounts[species]++; setRefCounts[species]++;
            individualsWithVariant[i] = 0;
        }
        // The second allele in this individual
        if (genotypes[i][2] == '1') {
            overall++;
            setAltCounts[species]++; setAlleleCounts[species]++;
            individualsWithVariant[i]++;
        } else if (genotypes[i][2] == '0') {
            setAlleleCounts[species]++; setRefCounts[species]++;
        }
    }
}

void GeneralSetCountsWithComplements::getComplementCounts(const std::vector<string>& populationsToUse) {
    // Now fill in the allele frequencies
    for(std::map<string,int>::iterator it = setAltCounts.begin(); it != setAltCounts.end(); ++it) {
        int complementAAFCount = 0; int complementTotalAlleleCount = 0;
        for (int i = 0; i < populationsToUse.size(); i++) {
            if (it->first != populationsToUse[i]) {
                complementAAFCount += setAltCounts.at(populationsToUse[i]);
                complementTotalAlleleCount += setAlleleCounts.at(populationsToUse[i]);
            }
        }
        setAAFsComplement[it->first] = (double)complementAAFCount/complementTotalAlleleCount;
        setAlleleCountsComplement[it->first] = complementTotalAlleleCount;
    }
}





void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

void splitToDouble(const std::string &s, char delim, std::vector<double> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(stringToDouble(item));
    }
}

std::vector<double> splitToDouble(const std::string &s, char delim) {
    std::vector<double> elems;
    splitToDouble(s, delim, elems);
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

// Initialize a matrix 
void initialize_matrix_double(std::vector<std::vector<double> >& m, int m_rows, int m_columns) {
    if (m_columns == 0)
        m_columns = m_rows;
    for (int i = 0; i < m_columns; i++) {
        std::vector<double> v(m_rows,0);
        m.push_back(v);
    }       
}

// Initialize a matrix 
void initialize_matrix_int(std::vector<std::vector<int> >& m, int m_rows, int m_columns) {
    if (m_columns == 0)
        m_columns = m_rows;
    for (int i = 0; i < m_columns; i++) {
        std::vector<int> v(m_rows,0);
        m.push_back(v);
    }       
}

// Remove a single file extension from the filename
std::string stripPath(const std::string& filename)
{
    size_t slashPos = filename.find_last_of('/');
    if(slashPos == std::string::npos)
        return filename; // no path
    else
        return filename.substr(slashPos+1);
}


// Remove a single file extension from the filename
std::string stripExtension(const std::string& filename)
{
    size_t suffixPos = filename.find_last_of('.');
    if(suffixPos == std::string::npos)
        return filename; // no suffix
    else
        return filename.substr(0, suffixPos);
}

bool isDPinfo (std::string infoField) {
    if (infoField.find("DP=") == string::npos)
        return false;
    else
        return true;
}

bool isFSinfo (std::string infoField) {
    if (infoField.find("FS=") == string::npos)
        return false;
    else
        return true;
}

double stringToDouble(std::string s) {
    double d;
    std::stringstream ss(s); //turn the string into a stream
    ss >> d; //convert
    return d;
}

// Look for MQSB: Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)
bool isMQSBinfo (std::string infoField) {
    if (infoField.find("MQSB=") == string::npos)
        return false;
    else
        return true;
}

// Look for SGB: Segregation based metric (bigger is better)
bool isSGBinfo (std::string infoField) {
    if (infoField.find("SGB=") == string::npos)
        return false;
    else
        return true;
}


double calculateInbreedingCoefficient(std::vector<int>& individualsWithVariant) {
    int naa = 0; int nAa = 0; int nAA = 0;
    for (std::vector<int>::size_type i = 0; i != individualsWithVariant.size(); i++) {
        if (individualsWithVariant[i] == 0) naa++;
        if (individualsWithVariant[i] == 1) nAa++;
        if (individualsWithVariant[i] == 2) nAA++;
    }
    
    // Get the proportions of alt-hom and hets
    double pAA = (double)nAA / individualsWithVariant.size();
    double pAa = (double)nAa / individualsWithVariant.size();
    
    // Allele frequencies
    double p = pAA + (0.5 * pAa);
    double q = 1 - p;
    
    // Get the Hardy-Weinberg prediction for expected number of heterozygotes:
    double HWAa = 2*p*q;

    
    // Get the inbreeding coefficient
    double F = (HWAa - pAa) / HWAa;
    return F;
}

double calculateChiSqPvalForInbreeding(std::vector<int>& individualsWithVariant) {
    // Get the observed numbers of genotypes
    int naa = 0; int nAa = 0; int nAA = 0;
    for (std::vector<int>::size_type i = 0; i != individualsWithVariant.size(); i++) {
        if (individualsWithVariant[i] == 0) naa++;
        if (individualsWithVariant[i] == 1) nAa++;
        if (individualsWithVariant[i] == 2) nAA++;
    }
    
    // Get the observed proportions of alt-hom and hets
    double pAA = (double)nAA / individualsWithVariant.size();
    double pAa = (double)nAa / individualsWithVariant.size();
    // Allele frequencies
    double p = pAA + (0.5 * pAa);
    double q = 1 - p;
    // Get the Hardy-Weinberg prediction for expected proportion of heterozygotes:
    double HWAa = 2*p*q;
    
    // Get the expected numbers of genotypes:
    double expAA = pow(p, 2) * individualsWithVariant.size();
    double expAa = HWAa * individualsWithVariant.size();
    double expaa = pow(q, 2) * individualsWithVariant.size();
    
    // Perform the chi-Sqared test
    std::vector<double> observed{static_cast<double>(nAA),static_cast<double>(nAa),static_cast<double>(naa)};
    std::vector<double> expected{expAA,expAa,expaa};
    double pVal = pearson_chi_sq_goodness_of_fit(observed, expected, 1);
    return pVal;
}



void randomisePhase(char& v1, char& v2) {
    if (v1 != v2) {
        char temp;
        double r = ((double) rand() / (RAND_MAX));
        if (r > 0.5) {
            temp = v1;
            v1 = v2; v2 = temp;
        }
    }
}



void getThisVariantCountsSimple(const std::vector<std::string>& fields, Counts* thisVariantCounts) {
    int numSamples = (int)fields.size()-NUM_NON_GENOTYPE_COLUMNS;
    thisVariantCounts->individualsWithVariant.assign(numSamples,0);
    thisVariantCounts->missingGenotypesPerIndividual.assign(numSamples,false);
    thisVariantCounts->haplotypesWithVariant.assign(numSamples*2,0);
    //std::cerr << "Fields: " << (fields.size()-NUM_NON_GENOTYPE_COLUMNS) << std::endl;
    std::vector<std::string> format = split(fields[8], ':');
    std::vector<std::string> altAlleles = split(fields[4], ',');
    const int numAltAlleles = (int)altAlleles.size();
    
    int alleleAsMissing = -1;
    if (numAltAlleles > 1) {
        string lastAlt = altAlleles.back();
        if (lastAlt == "*") {
            alleleAsMissing = numAltAlleles;
            if ((numAltAlleles - 1) > 1) {
                thisVariantCounts->bIsMultiallelic = true;
            }
        } else {
            if (numAltAlleles > 1) {
                thisVariantCounts->bIsMultiallelic = true;
            }
        }
    }
    
    // Find the position of GQ (genotype quality) in the genotypeData vector below
    if (fields[NUM_NON_GENOTYPE_COLUMNS][1] == '|') { thisVariantCounts->bPhased = true; }
    for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
        char v1 = fields[i][0]; char v2 = fields[i][2];
        if (thisVariantCounts->bPhased == false) randomisePhase(v1,v2);
        int v1int = v1 - '0'; int v2int = v2 - '0';
        if (v1 == '.' || v2 == '.' || v1int == alleleAsMissing || v2int == alleleAsMissing) {
            thisVariantCounts->bAnyMissingGenotypes = true;
            thisVariantCounts->missingGenotypesPerIndividual[i- NUM_NON_GENOTYPE_COLUMNS] = true;
        } else {
            thisVariantCounts->haplotypesWithVariant[2*(i-NUM_NON_GENOTYPE_COLUMNS)] = v1int;
            thisVariantCounts->haplotypesWithVariant[2*(i-NUM_NON_GENOTYPE_COLUMNS)+1] = v2int;
        }
        // Multiallelic sites will be ignored in the genotype record
        if (!thisVariantCounts->bIsMultiallelic) {
            if (v1 == '1') {
                thisVariantCounts->overall++;
                thisVariantCounts->individualsWithVariant[i- NUM_NON_GENOTYPE_COLUMNS]++;
            }
            if (v2 == '1') {
                thisVariantCounts->overall++;
                thisVariantCounts->individualsWithVariant[i-NUM_NON_GENOTYPE_COLUMNS]++;
            }
        }
    }
}



Counts getThisVariantCounts(const std::vector<std::string>& fields) {
    Counts thisVariantCounts;
    bool hasGQ = false; bool hasDP = false;
    thisVariantCounts.individualsWithVariant.assign((fields.size()-NUM_NON_GENOTYPE_COLUMNS),0);
    thisVariantCounts.haplotypesWithVariant.assign((fields.size()-NUM_NON_GENOTYPE_COLUMNS)*2,0);
    //std::cerr << "Fields: " << (fields.size()-NUM_NON_GENOTYPE_COLUMNS) << std::endl;
    // Find the position of DP (per sample read depth) in the genotypeData vector below
    std::vector<std::string> format = split(fields[8], ':');
    std::vector<std::string>::iterator DPit; int DPi = std::numeric_limits<int>::min();
    DPit = find (format.begin(), format.end(), "DP");
    if (DPit == format.end()) {
        // std::cerr << "This variant hasn't got associated per-sample DP info" << std::endl;
    } else {
        DPi = (int)std::distance( format.begin(), DPit );
        hasDP = true;
    }
    // Find the position of GQ (genotype quality) in the genotypeData vector below
    std::vector<std::string>::iterator GQit; int GQi = std::numeric_limits<int>::min();
    GQit = find (format.begin(), format.end(), "GQ");
    if (GQit == format.end()) {
        // std::cerr << "This variant hasn't got associated per-sample GQ info" << std::endl;
    } else {
        GQi = (int)std::distance( format.begin(), GQit );
        hasGQ = true;
    }
    
    if (fields[NUM_NON_GENOTYPE_COLUMNS][1] == '|') { thisVariantCounts.bPhased = true; }
    
    for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
        char v1 = fields[i][0]; char v2 = fields[i][2];
        if (v1 == '.' || v2 == '.') {
            thisVariantCounts.bAnyMissingGenotypes = true;
            thisVariantCounts.missingGenotypesPerIndividual[i- NUM_NON_GENOTYPE_COLUMNS] = true;
        }
        if (thisVariantCounts.bPhased == false) {
            if ((v1 == '0' && v2 == '1') || (v1 == '1' && v2 == '0')) {
            double r = ((double) rand() / (RAND_MAX));
            if (r > 0.5) {
                v1 = '0'; v2 = '1';
            } else {
                v1 = '1'; v2 = '0';
            }
            }
        }
        
        if (v1 == '1') {
            thisVariantCounts.overall++;
            thisVariantCounts.individualsWithVariant[i- NUM_NON_GENOTYPE_COLUMNS]++;
            thisVariantCounts.haplotypesWithVariant[2*(i-NUM_NON_GENOTYPE_COLUMNS)]++;
        }
        if (v2 == '1') {
            thisVariantCounts.overall++;
            thisVariantCounts.individualsWithVariant[i-NUM_NON_GENOTYPE_COLUMNS]++;
            thisVariantCounts.haplotypesWithVariant[2*(i-NUM_NON_GENOTYPE_COLUMNS)+1]++;
        }
            
        std::vector<std::string> genotypeData = split(fields[i], ':');
        
        // read depth at the variant site per individual
        if (hasDP && fields[i][0] != '.') {
            if (atoi(genotypeData[DPi].c_str()) < thisVariantCounts.minimumDepthInAnIndividual) {
                thisVariantCounts.minimumDepthInAnIndividual = atoi(genotypeData[DPi].c_str());
            }
            thisVariantCounts.depthPerIndividual.push_back(atoi(genotypeData[DPi].c_str()));
        }
        // genotype quality at the variant site per individual
        if (hasGQ && fields[i][0] != '.') {
            thisVariantCounts.genotypeQualitiesPerIndividual.push_back(atoi(genotypeData[GQi].c_str()));
        }
    }
    // Also get overall depth for this variant
    std::vector<std::string> info = split(fields[7], ';');
    std::vector<std::string>::iterator overallDPit; int overallDPi = std::numeric_limits<int>::min();
    overallDPit = find_if(info.begin(), info.end(), isDPinfo);
    if (overallDPit == info.end()) {
        // std::cerr << "This variant hasn't got associated overall DP info" << std::endl;
        thisVariantCounts.overallDepth = 0;
    } else {
        overallDPi = (int)std::distance( info.begin(), overallDPit );
        std::vector<std::string> overallDP = split(info[overallDPi], '=');
        thisVariantCounts.overallDepth = atoi((overallDP.back()).c_str());
    }
    
    // Find the position of SGB (Segregation based metric) in the genotypeData vector below
    std::vector<std::string>::iterator SGBit; int SGBi = std::numeric_limits<int>::min();
    SGBit = find_if(format.begin(), format.end(), isSGBinfo);
    if (SGBit != format.end()) {
        // std::cerr << "This variant hasn't got associated per-sample GQ info" << std::endl;
    } else {
        SGBi = (int)std::distance( format.begin(), SGBit );
        std::vector<std::string> overallSGB = split(info[SGBi], '=');
        thisVariantCounts.SGB = stringToDouble(overallSGB.back());
    }
    
    // And get FS (phred-scaled strand-bias p-val) for this variant
    std::vector<std::string>::iterator FSit; int FSi = std::numeric_limits<int>::min();
    FSit = find_if(info.begin(), info.end(), isFSinfo);
    if (FSit == info.end()) {  // Or at least MQSB: Mann-Whitney U test of Mapping Quality vs Strand Bias
        FSit = find_if(info.begin(), info.end(), isMQSBinfo);
        if (FSit != info.end()) {
            FSi = (int)std::distance( info.begin(), FSit );
            std::vector<std::string> overallFS = split(info[FSi], '=');
            thisVariantCounts.MQSBpval = overallFS.back();
        } else {
            // std::cerr << "This variant hasn't got associated FS (strand-bias) info" << std::endl;
        }
    } else {
        FSi = (int)std::distance( info.begin(), FSit );
        std::vector<std::string> overallFS = split(info[FSi], '=');
        thisVariantCounts.FSpval =  overallFS.back();
    }

    // And the inbreeding coefficient
    thisVariantCounts.inbreedingCoefficient = calculateInbreedingCoefficient(thisVariantCounts.individualsWithVariant);
    
    return thisVariantCounts;
}

ThreeSetCounts getThreeSetVariantCounts(const std::vector<std::string>& fields, const std::vector<size_t>& set1_loci, const std::vector<size_t>& set2_loci, const std::vector<size_t>& set3_loci, const std::string& AA) {
    ThreeSetCounts thisVariantCounts;
    thisVariantCounts.individualsWithVariant.assign((fields.size()-NUM_NON_GENOTYPE_COLUMNS),0);
    // std::cerr << fields[0] << "\t" << fields[1] << std::endl;
    for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
        if (fields[i][0] == '1') {
            thisVariantCounts.overall++;
            if (std::find(set1_loci.begin(), set1_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set1_loci.end()) { thisVariantCounts.set1AltCount++; }
            if (std::find(set2_loci.begin(), set2_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set2_loci.end()) { thisVariantCounts.set2AltCount++; }
            if (std::find(set3_loci.begin(), set3_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set3_loci.end()) { thisVariantCounts.set3AltCount++; }
            thisVariantCounts.individualsWithVariant[i- NUM_NON_GENOTYPE_COLUMNS]++;
        }
        if (fields[i][2] == '1') {
            thisVariantCounts.overall++;
            if (std::find(set1_loci.begin(), set1_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set1_loci.end()) { thisVariantCounts.set1AltCount++; }
            if (std::find(set2_loci.begin(), set2_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set2_loci.end()) { thisVariantCounts.set2AltCount++; }
            if (std::find(set3_loci.begin(), set3_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set3_loci.end()) { thisVariantCounts.set3AltCount++; }
            thisVariantCounts.individualsWithVariant[i-NUM_NON_GENOTYPE_COLUMNS]++;
        }
    }
    thisVariantCounts.set1RefCount = (int)(set1_loci.size() * 2) - thisVariantCounts.set1AltCount;
    thisVariantCounts.set2RefCount = (int)(set2_loci.size() * 2) - thisVariantCounts.set2AltCount;
    thisVariantCounts.set3RefCount = (int)(set3_loci.size() * 2) - thisVariantCounts.set3AltCount;
    
    thisVariantCounts.set1AltAF = (double)thisVariantCounts.set1AltCount/(set1_loci.size() * 2);
    thisVariantCounts.set2AltAF = (double)thisVariantCounts.set2AltCount/(set2_loci.size() * 2);
    thisVariantCounts.set3AltAF = (double)thisVariantCounts.set3AltCount/(set3_loci.size() * 2);
    
    // Fill in derived allele frequencies if possible
    if (AA == "ref") {
        thisVariantCounts.set1daAF = (double)thisVariantCounts.set1AltCount/(set1_loci.size() * 2);
        thisVariantCounts.set2daAF = (double)thisVariantCounts.set2AltCount/(set2_loci.size() * 2);
        thisVariantCounts.set3daAF = (double)thisVariantCounts.set3AltCount/(set3_loci.size() * 2);
    } else if (AA == "alt") {
        thisVariantCounts.set1daAF = (double)thisVariantCounts.set1RefCount/(set1_loci.size() * 2);
        thisVariantCounts.set2daAF = (double)thisVariantCounts.set2RefCount/(set2_loci.size() * 2);
        thisVariantCounts.set3daAF = (double)thisVariantCounts.set3RefCount/(set3_loci.size() * 2);
    } else if (AA == "N") {
    } else {
        std::cerr << "Error: Derived allele can only be \"ref\" or \"alt\"" << std::endl;
        exit(1);
    }
    
    return thisVariantCounts;
}


FourSetCounts getFourSetVariantCounts(const std::vector<std::string>& fields, const std::vector<size_t>& set1_loci, const std::vector<size_t>& set2_loci, const std::vector<size_t>& set3_loci, const std::vector<size_t>& set4_loci, const std::string& AA) {
    FourSetCounts thisVariantCounts;
    thisVariantCounts.individualsWithVariant.assign((fields.size()-NUM_NON_GENOTYPE_COLUMNS),0);
    // std::cerr << fields[0] << "\t" << fields[1] << std::endl;
    for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
        if (fields[i][0] == '1') {
            thisVariantCounts.overall++;
            if (std::find(set1_loci.begin(), set1_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set1_loci.end()) { thisVariantCounts.set1AltCount++; }
            if (std::find(set2_loci.begin(), set2_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set2_loci.end()) { thisVariantCounts.set2AltCount++; }
            if (std::find(set3_loci.begin(), set3_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set3_loci.end()) { thisVariantCounts.set3AltCount++; }
            if (std::find(set4_loci.begin(), set4_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set4_loci.end()) { thisVariantCounts.set4AltCount++; }
            thisVariantCounts.individualsWithVariant[i- NUM_NON_GENOTYPE_COLUMNS]++;
        }
        if (fields[i][2] == '1') {
            thisVariantCounts.overall++;
            if (std::find(set1_loci.begin(), set1_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set1_loci.end()) { thisVariantCounts.set1AltCount++; }
            if (std::find(set2_loci.begin(), set2_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set2_loci.end()) { thisVariantCounts.set2AltCount++; }
            if (std::find(set3_loci.begin(), set3_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set3_loci.end()) { thisVariantCounts.set3AltCount++; }
            if (std::find(set4_loci.begin(), set4_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set4_loci.end()) { thisVariantCounts.set4AltCount++; }
            thisVariantCounts.individualsWithVariant[i-NUM_NON_GENOTYPE_COLUMNS]++;
        }
    }
    thisVariantCounts.set1RefCount = (int)(set1_loci.size() * 2) - thisVariantCounts.set1AltCount;
    thisVariantCounts.set2RefCount = (int)(set2_loci.size() * 2) - thisVariantCounts.set2AltCount;
    thisVariantCounts.set3RefCount = (int)(set3_loci.size() * 2) - thisVariantCounts.set3AltCount;
    thisVariantCounts.set4RefCount = (int)(set4_loci.size() * 2) - thisVariantCounts.set4AltCount;
    
    thisVariantCounts.set1AltAF = (double)thisVariantCounts.set1AltCount/(set1_loci.size() * 2);
    thisVariantCounts.set2AltAF = (double)thisVariantCounts.set2AltCount/(set2_loci.size() * 2);
    thisVariantCounts.set3AltAF = (double)thisVariantCounts.set3AltCount/(set3_loci.size() * 2);
    thisVariantCounts.set4AltAF = (double)thisVariantCounts.set4AltCount/(set4_loci.size() * 2);

    // Fill in derived allele frequencies if possible
    if (AA == "ref") {
        thisVariantCounts.set1daAF = (double)thisVariantCounts.set1AltCount/(set1_loci.size() * 2);
        thisVariantCounts.set2daAF = (double)thisVariantCounts.set2AltCount/(set2_loci.size() * 2);
        thisVariantCounts.set3daAF = (double)thisVariantCounts.set3AltCount/(set3_loci.size() * 2);
        thisVariantCounts.set4daAF = (double)thisVariantCounts.set4AltCount/(set4_loci.size() * 2);
    } else if (AA == "alt") {
        thisVariantCounts.set1daAF = (double)thisVariantCounts.set1RefCount/(set1_loci.size() * 2);
        thisVariantCounts.set2daAF = (double)thisVariantCounts.set2RefCount/(set2_loci.size() * 2);
        thisVariantCounts.set3daAF = (double)thisVariantCounts.set3RefCount/(set3_loci.size() * 2);
        thisVariantCounts.set4daAF = (double)thisVariantCounts.set4RefCount/(set4_loci.size() * 2);
    } else if (AA == "N") {
    } else {
        std::cerr << "Error: Derived allele can only be \"ref\" or \"alt\"" << std::endl;
        exit(1);
    }
    
    return thisVariantCounts;
}

ThreeSetCounts getThreeSetVariantCountsAA4(const std::vector<std::string>& fields, const std::vector<size_t>& set1_loci, const std::vector<size_t>& set2_loci, const std::vector<size_t>& set3_loci, const std::vector<size_t>& AA_loci) {
    ThreeSetCounts thisVariantCounts;
    int AAaltCount = 0;
    thisVariantCounts.individualsWithVariant.assign((fields.size()-NUM_NON_GENOTYPE_COLUMNS),0);
    // std::cerr << fields[0] << "\t" << fields[1] << std::endl;
    for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
        if (fields[i][0] == '1') {
            thisVariantCounts.overall++;
            if (std::find(set1_loci.begin(), set1_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set1_loci.end()) { thisVariantCounts.set1AltCount++; }
            if (std::find(set2_loci.begin(), set2_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set2_loci.end()) { thisVariantCounts.set2AltCount++; }
            if (std::find(set3_loci.begin(), set3_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set3_loci.end()) { thisVariantCounts.set3AltCount++; }
            if (std::find(AA_loci.begin(), AA_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != AA_loci.end()) { AAaltCount++; }
            thisVariantCounts.individualsWithVariant[i- NUM_NON_GENOTYPE_COLUMNS]++;
        }
        if (fields[i][2] == '1') {
            thisVariantCounts.overall++;
            if (std::find(set1_loci.begin(), set1_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set1_loci.end()) { thisVariantCounts.set1AltCount++; }
            if (std::find(set2_loci.begin(), set2_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set2_loci.end()) { thisVariantCounts.set2AltCount++; }
            if (std::find(set3_loci.begin(), set3_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set3_loci.end()) { thisVariantCounts.set3AltCount++; }
            if (std::find(AA_loci.begin(), AA_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != AA_loci.end()) { AAaltCount++; }
            thisVariantCounts.individualsWithVariant[i-NUM_NON_GENOTYPE_COLUMNS]++;
        }
    }
    thisVariantCounts.set1RefCount = (int)(set1_loci.size() * 2) - thisVariantCounts.set1AltCount;
    thisVariantCounts.set2RefCount = (int)(set2_loci.size() * 2) - thisVariantCounts.set2AltCount;
    thisVariantCounts.set3RefCount = (int)(set3_loci.size() * 2) - thisVariantCounts.set3AltCount;
    //thisVariantCounts.set4RefCount = (int)(set4_loci.size() * 2) - thisVariantCounts.set4AltCount;
    
    thisVariantCounts.set1AltAF = (double)thisVariantCounts.set1AltCount/(set1_loci.size() * 2);
    thisVariantCounts.set2AltAF = (double)thisVariantCounts.set2AltCount/(set2_loci.size() * 2);
    thisVariantCounts.set3AltAF = (double)thisVariantCounts.set3AltCount/(set3_loci.size() * 2);
    //thisVariantCounts.set4AltAF = (double)thisVariantCounts.set4AltCount/(set4_loci.size() * 2);
    
    // Fill in derived allele frequencies if possible
    if (AAaltCount == 0) {
        thisVariantCounts.set1daAF = (double)thisVariantCounts.set1AltCount/(set1_loci.size() * 2);
        thisVariantCounts.set2daAF = (double)thisVariantCounts.set2AltCount/(set2_loci.size() * 2);
        thisVariantCounts.set3daAF = (double)thisVariantCounts.set3AltCount/(set3_loci.size() * 2);
       // thisVariantCounts.set4daAF = (double)thisVariantCounts.set4AltCount/(set4_loci.size() * 2);
    } else if (AAaltCount == 2) {
        thisVariantCounts.set1daAF = (double)thisVariantCounts.set1RefCount/(set1_loci.size() * 2);
        thisVariantCounts.set2daAF = (double)thisVariantCounts.set2RefCount/(set2_loci.size() * 2);
        thisVariantCounts.set3daAF = (double)thisVariantCounts.set3RefCount/(set3_loci.size() * 2);
        //thisVariantCounts.set4daAF = (double)thisVariantCounts.set4RefCount/(set4_loci.size() * 2);
    } else if (AAaltCount == 1) {
    } else {
        std::cerr << "Error: Outgroup can only be one individual here" << std::endl;
        exit(1);
    }
    
    return thisVariantCounts;
}





bool testBiallelic(const std::string& altField) {
    std::vector<std::string> altVector = split(altField, ',');
    if (altVector.size() == 1) { return true; }
    else { return false; }
}

bool testOverallReadDepth(const int maxReadDepth,const int minReadDepth, const std::string& infoField) {
    std::vector<std::string> info = split(infoField, ';');
    if (info[0] == "INDEL") {
        split(info[1], '=', info);
    } else {
        split(info[0], '=', info);
    }
    int DP = atoi((info.back()).c_str());
    if (DP <= maxReadDepth && DP >= minReadDepth) { return true; }
    else { return false; }
}


// filter out sites where more than maxNumHet individuals are heterozygous 
bool testMaxNumHet(FilterResult& result, std::vector<int>& depthsHetFailed, std::vector<int>& depthsHetPassed, std::vector<int>& numVariantsPerHetCount, int maxNumHet, std::vector<std::vector<int> >& num_indiv_het_vs_depth) {
    // filter out sites where more than MAX_NUM_HET individuals are heterozygous
    int num_hets = 0;
    for (std::vector<std::vector<int> >::size_type i = 0; i < result.counts.individualsWithVariant.size(); i++) {
        if (result.counts.individualsWithVariant[i] == 1)
            num_hets++;
    }
    numVariantsPerHetCount[num_hets]++;
    
    // Randomly sample 1% of sites for num_indiv_het vs depth scatterplot
    double rn = ((double) rand() / RAND_MAX);
    if (rn < 0.01) {
        std::vector<int> this_num_het_depth; this_num_het_depth.push_back(num_hets); this_num_het_depth.push_back(result.counts.overallDepth); 
        num_indiv_het_vs_depth.push_back(this_num_het_depth);
    }
    
    if (num_hets > maxNumHet) {
        depthsHetFailed.push_back(result.counts.overallDepth);
        return false;
    } else {
        depthsHetPassed.push_back(result.counts.overallDepth);
        return true;
    }
}

// Does the same as R function table
std::map<int, int> tabulateVector(std::vector<int>& vec) {
    std::vector<int> vecCopy(vec); 
    std::sort(vecCopy.begin(), vecCopy.end());
    std::vector<int>::iterator it = std::unique(vecCopy.begin(), vecCopy.end());
    vecCopy.resize(std::distance(vecCopy.begin(), it));
    
    std::map<int, int>  table;
    //int pos = 0;
    for (std::vector<int>::size_type i = 0; i != vecCopy.size(); i++) {
        int mycount = std::count(vec.begin(), vec.end(), vecCopy[i]);
        table[vecCopy[i]] = mycount;
        //pos = pos + mycount;
    }
    return table;
}



// Move all doubleton counts to the bottom left corner of the matrix
void rearrange_doubletons(std::vector<std::vector<int> >& doubletons){
    for (std::vector<std::vector<int> >::size_type i = 0; i < doubletons.size(); i++) {
        for (int j = 0; j < i; j++) {
            doubletons[i][j] = doubletons[i][j] + doubletons[j][i];
            doubletons[j][i] = 0;
        }    
    }
}

// For checking if variants fixed in Massoko are also fixed in Malawi
void massokoMalawiSharing(const FilterResult& result, MassokoMalawiResult& sharingResult) {
    if (result.counts.individualsWithVariant[0] == 2 || result.counts.individualsWithVariant[1] == 2 ) {
        for (std::vector<int>::size_type i = 0; i < sharingResult.hetsWithBlue.size(); i++) {
            //std::cerr << ": " << std::endl;
            if (result.counts.individualsWithVariant[i+12] == 0) {
                sharingResult.absentWithBlue[i]++;
            }
            if (result.counts.individualsWithVariant[i+12] == 1) {
                sharingResult.hetsWithBlue[i]++;
            }
            if (result.counts.individualsWithVariant[i+12] == 2) {
                sharingResult.homsWithBlue[i]++;
            }
            if (result.counts.individualsWithVariant[i+12] == 2) {
                sharingResult.absentWithYellow[i]++;
            }
            if (result.counts.individualsWithVariant[i+12] == 1) {
                sharingResult.hetsWithYellow[i]++;
            }
            if (result.counts.individualsWithVariant[i+12] == 0) {
                sharingResult.homsWithYellow[i]++;
            }
        }
    } else if (result.counts.individualsWithVariant[6] == 2 || result.counts.individualsWithVariant[7] == 2) {
        for (std::vector<int>::size_type i = 0; i < sharingResult.hetsWithYellow.size(); i++) {
            if (result.counts.individualsWithVariant[i+12] == 0) {
                sharingResult.absentWithYellow[i]++;
            }
            if (result.counts.individualsWithVariant[i+12] == 1) {
                sharingResult.hetsWithYellow[i]++;
            }
            if (result.counts.individualsWithVariant[i+12] == 2) {
                sharingResult.homsWithYellow[i]++;
            }
            if (result.counts.individualsWithVariant[i+12] == 2) {
                sharingResult.absentWithBlue[i]++;
            }
            if (result.counts.individualsWithVariant[i+12] == 1) {
                sharingResult.hetsWithBlue[i]++;
            }
            if (result.counts.individualsWithVariant[i+12] == 0) {
                sharingResult.homsWithBlue[i]++;
            }
        }
    } else {
        std::cerr << "There is a variant that does not appear fixed in Massoko: " << std::endl;
    }
}

std::vector<string> readSampleNamesFromTextFile(const std::string& sampleNameFile) {
    std::vector<string> sampleNames;
    std::ifstream* sampleFile = new std::ifstream(sampleNameFile.c_str());
    string line;
    while (getline(*sampleFile, line)) {
        sampleNames.push_back(line);
    }
    return sampleNames;
}

std::string suffix(const std::string& seq, size_t len)
{
    assert(seq.length() >= len);
    return seq.substr(seq.length() - len);
}


// Returns true if the filename has an extension indicating it is compressed
bool isGzip(const std::string& filename)
{
    size_t suffix_length = sizeof(GZIP_EXT) - 1;
    
    // Assume files without an extension are not compressed
    if(filename.length() < suffix_length)
        return false;
    
    std::string extension = suffix(filename, suffix_length);
    return extension == GZIP_EXT;
}

// Ensure a filehandle is open
void assertFileOpen(std::ifstream& fh, const std::string& fn)
{
    if(!fh.is_open())
    {
        std::cerr << "Error: could not open " << fn << " for read\n";
        exit(EXIT_FAILURE);
    }
}

// Ensure a filehandle is open
void assertFileOpen(std::ofstream& fh, const std::string& fn)
{
    if(!fh.is_open())
    {
        std::cerr << "Error: could not open " << fn << " for write\n";
        exit(EXIT_FAILURE);
    }
}

//
void assertGZOpen(gzstreambase& gh, const std::string& fn)
{
    if(!gh.good())
    {
        std::cerr << "Error: could not open " << fn << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Open a file that may or may not be gzipped for reading
// The caller is responsible for freeing the handle
std::istream* createReader(const std::string& filename, std::ios_base::openmode mode)
{
    if(isGzip(filename))
    {
        igzstream* pGZ = new igzstream(filename.c_str(), mode);
        assertGZOpen(*pGZ, filename);
        return pGZ;
    }
    else
    {
        std::ifstream* pReader = new std::ifstream(filename.c_str(), mode);
        assertFileOpen(*pReader, filename);
        return pReader;
    }
}

// Open a file that may or may not be gzipped for writing
// The caller is responsible for freeing the handle
std::ostream* createWriter(const std::string& filename,
                           std::ios_base::openmode mode)
{
    if(isGzip(filename))
    {
        ogzstream* pGZ = new ogzstream(filename.c_str(), mode);
        assertGZOpen(*pGZ, filename);
        return pGZ;
    }
    else
    {
        std::ofstream* pWriter = new std::ofstream(filename.c_str(), mode);
        assertFileOpen(*pWriter, filename);
        return pWriter;
    }
}



void print80bpPerLineStdOut(std::ostream& outStream, string toPrint) {
    string::size_type lines = toPrint.length() / 80;
    for (string::size_type j = 0; j <= lines; j++) {
        outStream << toPrint.substr(j*80,80) << std::endl;
    }
}

void print80bpPerLineFile(std::ofstream*& outFile, string toPrint) {
    string::size_type lines = toPrint.length() / 80;
    for (string::size_type j = 0; j <= lines; j++) {
        *outFile << toPrint.substr(j*80,80) << std::endl;
    }
}

std::vector<size_t> locateSet(std::vector<std::string>& sample_names, const std::vector<std::string>& set) {
    std::vector<size_t> setLocs;
    for (std::vector<std::string>::size_type i = 0; i != set.size(); i++) {
        std::vector<std::string>::iterator it = std::find(sample_names.begin(), sample_names.end(), set[i]);
        if (it == sample_names.end()) {
            std::cerr << "Did not find the sample: \"" << set[i] << "\"" << std::endl;
            std::cerr << "Did not find the sample: \"" << sample_names[43] << "\"" << std::endl;
            print_vector(sample_names, std::cerr,',');
        } else {
            size_t loc = std::distance(sample_names.begin(), it);
            setLocs.push_back(loc);
        }
    }
    return setLocs;
}

std::vector<size_t> complementIndices(const size_t fullVectorSize, const std::vector<size_t>& originalIndices) {
    std::vector<size_t> indexComplement;
    std::vector<bool> isInOriginal(fullVectorSize, false);
    for (const size_t j : originalIndices) {
        isInOriginal[j] = true;
    }
    indexComplement.reserve(fullVectorSize - originalIndices.size());
    for (size_t i = 0; i < fullVectorSize; i++) {
        if (!isInOriginal[i])
        indexComplement.push_back(i);
    }
    return indexComplement;
}

size_t locateOneSample(std::vector<std::string>& sample_names, const std::string toFind) {
    size_t pos = 0;
    std::vector<std::string>::iterator it = std::find(sample_names.begin(), sample_names.end(), toFind);
    if (it == sample_names.end()) {
        std::cerr << "Did not find the sample: " << toFind << std::endl;
    } else {
        pos = std::distance(sample_names.begin(), it);
    }
    return pos;
}

std::map<string,string> readMultiFastaToMap(const string& fileName) {
    std::map<string, string> fastaSeqs;
    string line;
    std::ifstream* fastaFile = new std::ifstream(fileName.c_str());
    getline(*fastaFile, line);
    std::vector<std::string> fields = split(line, ' ');
    string currentScaffold = fields[0].substr(1,string::npos);
    fastaSeqs[currentScaffold] = ""; fastaSeqs[currentScaffold].reserve(50000000);
    while (getline(*fastaFile, line)) {
        if (line[0] != '>') {
            fastaSeqs[currentScaffold].append(line);
        } else {
            int l = (int)fastaSeqs[currentScaffold].length();
            std::cerr << currentScaffold << " length: " << l << std::endl;
            fastaSeqs[currentScaffold].shrink_to_fit();
            std::vector<std::string> fields = split(line, ' ');
            currentScaffold = fields[0].substr(1,string::npos);
            fastaSeqs[currentScaffold] = ""; fastaSeqs[currentScaffold].reserve(50000000);
        }
    }
    std::cerr << currentScaffold << " length: " << (int)fastaSeqs[currentScaffold].length() << std::endl;
    return fastaSeqs;
}

std::string readMultiFastaToOneString(const string& fileName, int bytes) {
    std::string fastaSeqs; fastaSeqs.reserve(bytes); fastaSeqs = "";
    string line;
    std::ifstream* fastaFile = new std::ifstream(fileName.c_str());
    getline(*fastaFile, line);
    string currentScaffold = line.substr(1,string::npos);
    while (getline(*fastaFile, line)) {
        if (line[0] != '>') {
            fastaSeqs.append(line);
        } else {
            continue;
        }
    }
    return fastaSeqs;
}


double getExpectedGenotype(const std::vector<double>& thisProbabilities) {
    double Egenotype = thisProbabilities[1] + 2*thisProbabilities[2];
    return Egenotype;
}

void transformFromPhred(std::vector<double>& thisLikelihoods) {

    thisLikelihoods[0] = pow(10,-(thisLikelihoods[0]/10.0));
    thisLikelihoods[1] = pow(10,-(thisLikelihoods[1]/10.0));
    thisLikelihoods[2] = pow(10,-(thisLikelihoods[2]/10.0));
}

void transformFromGL(std::vector<double>& thisLikelihoods) {

    thisLikelihoods[0] = pow(10,(thisLikelihoods[0]/10.0));
    thisLikelihoods[1] = pow(10,(thisLikelihoods[1]/10.0));
    thisLikelihoods[2] = pow(10,(thisLikelihoods[2]/10.0));
}


std::vector<double> GeneralSetCounts::probabilitiesFromLikelihoods(const std::vector<double>& thisLikelihoods, const string& species) {
    std::vector<double> thisProbabilities; thisProbabilities.assign(3, 0.0);
    double multiple0 = thisLikelihoods[0]*setHWEpriorsFromAAFfromGT[species][0];
    double multiple1 = thisLikelihoods[1]*setHWEpriorsFromAAFfromGT[species][1];
    double multiple2 = thisLikelihoods[2]*setHWEpriorsFromAAFfromGT[species][2];
    double sum = multiple0 + multiple1 + multiple2;
    
    thisProbabilities[0] = multiple0/sum;
    thisProbabilities[1] = multiple1/sum;
    thisProbabilities[2] = multiple2/sum;
    
    return thisProbabilities;
}


void GeneralSetCounts::setHWEpriorsFromAFfromGT() {
    double AF;
    // Alternative allele frequencies
    for(std::map<string,double>::iterator it = setAAFs.begin(); it != setAAFs.end(); ++it) {
        if (it->second >= 0) AF = it->second; else AF = averageAAF; // This should be average of AFs across populations where it is known
        setHWEpriorsFromAAFfromGT[it->first][0] = pow((1-AF),2);
        setHWEpriorsFromAAFfromGT[it->first][1] = AF*(1-AF);
        setHWEpriorsFromAAFfromGT[it->first][2] = pow(AF,2);
    }
    // Derived allele frequencies
    for(std::map<string,double>::iterator it = setDAFs.begin(); it != setDAFs.end(); ++it) {
        if (it->second >= 0) AF = it->second; else AF = averageDAF; // This should be average of AFs across populations
        setHWEpriorsFromDAFfromGT[it->first][0] = pow((1-AF),2);
        setHWEpriorsFromDAFfromGT[it->first][1] = AF*(1-AF);
        setHWEpriorsFromDAFfromGT[it->first][2] = pow(AF,2);
    }
}

void GeneralSetCounts::getAFsFromGenotypeLikelihoodsOrProbabilities(const std::vector<std::string>& genotypeFields, const std::map<size_t, string>& posToSpeciesMap, const int likelihoodsOrProbabilitiesTagPosition) {
    if (likelihoodsProbabilitiesType == LikelihoodsProbabilitiesPL || likelihoodsProbabilitiesType == LikelihoodsProbabilitiesGL) {
        setHWEpriorsFromAFfromGT();
    }
    
    for (std::vector<std::string>::size_type i = 0; i < genotypeFields.size(); i++) {
       // std::cerr << "Here posToSpeciesMap.at(i): " << posToSpeciesMap.at(i) << std::endl;
        std::string species; try { species = posToSpeciesMap.at(i); } catch (const std::out_of_range& oor) {
            // std::cerr << "Here the species is not in the map" << std::endl;
            continue;
        }
       // std::cerr << genotypeFields[i] << std::endl;
        std::string thisLikelihoodsOrProbabilitiesString = split(genotypeFields[i], ':')[likelihoodsOrProbabilitiesTagPosition];
        if (thisLikelihoodsOrProbabilitiesString == "." || thisLikelihoodsOrProbabilitiesString == "0,0,0") continue;
        
        else {
            setAlleleProbCounts.at(species) += 2;
            std::vector<double> thisLikelihoodsOrProbabilities = splitToDouble(thisLikelihoodsOrProbabilitiesString,',');
            std::vector<double> thisProbabilities;
            switch (likelihoodsProbabilitiesType)
            {
                case LikelihoodsProbabilitiesPL:
                    transformFromPhred(thisLikelihoodsOrProbabilities);
                   // print_vector(thisLikelihoodsOrProbabilities, std::cerr);
                    thisProbabilities = probabilitiesFromLikelihoods(thisLikelihoodsOrProbabilities,species);
                    break;
                case LikelihoodsProbabilitiesGL: transformFromGL(thisLikelihoodsOrProbabilities);
                    thisProbabilities = probabilitiesFromLikelihoods(thisLikelihoodsOrProbabilities,species);
                    break;
                case LikelihoodsProbabilitiesGP:
                    thisProbabilities = thisLikelihoodsOrProbabilities;
                    break;
            }
            if (setAAFsFromLikelihoods.at(species) == -1) setAAFsFromLikelihoods.at(species) = 0;
            setAAFsFromLikelihoods.at(species) += getExpectedGenotype(thisProbabilities);
        }
    }
    
    for(std::map<string,double>::iterator it = setAAFsFromLikelihoods.begin(); it != setAAFsFromLikelihoods.end(); ++it) {
        if (setAAFsFromLikelihoods.at(it->first) != -1) {
            double AF = it->second/setAlleleProbCounts.at(it->first);
            it->second = AF;
            if (AAint == AncestralAlleleRef) {
                setDAFsFromLikelihoods.at(it->first) = AF;
            } else if (AAint == AncestralAlleleAlt) {
                setDAFsFromLikelihoods.at(it->first) = (1 - AF);
            }
        }
    }
     
}

int GeneralSetCounts::returnFormatTagPosition(std::vector<std::string>& format, const std::string& tag) {
    // Find the position of GQ (genotype quality) in the genotypeData vector below
    std::vector<std::string>::iterator TAGit; int TAGi = std::numeric_limits<int>::min();
    TAGit = find (format.begin(), format.end(), tag);
    if (TAGit == format.end()) {
        // std::cerr << "This variant hasn't got associated per-sample GQ info" << std::endl;
    } else {
        TAGi = (int)std::distance( format.begin(), TAGit );
        //hasGQ = true;
    }
    return TAGi;
}

int GeneralSetCounts::checkForGenotypeLikelihoodsOrProbabilities(const std::vector<std::string>& vcfLineFields) {
    std::vector<std::string> format = split(vcfLineFields[8], ':');
    if (format.size() == 1) return LikelihoodsProbabilitiesAbsent; // The GT tag must be present in the first place
    
    int likelihoodsOrProbabilitiesTagPosition = returnFormatTagPosition(format, "GP");
    if (likelihoodsOrProbabilitiesTagPosition != std::numeric_limits<int>::min()) { likelihoodsProbabilitiesType = LikelihoodsProbabilitiesGP; }
    else {
        likelihoodsOrProbabilitiesTagPosition = returnFormatTagPosition(format, "GL");
        if (likelihoodsOrProbabilitiesTagPosition != std::numeric_limits<int>::min()) { likelihoodsProbabilitiesType = LikelihoodsProbabilitiesGL; }
        else {
            likelihoodsOrProbabilitiesTagPosition = returnFormatTagPosition(format, "PL");
            if (likelihoodsOrProbabilitiesTagPosition != std::numeric_limits<int>::min()) { likelihoodsProbabilitiesType = LikelihoodsProbabilitiesPL; }
        }
    }
    return likelihoodsOrProbabilitiesTagPosition;
}



/*  INPUT MATRIX:
 In the bottom-left part of the matrix (addressed here as diffs_Hets_vs_Homs[i][j]) should be counts of sites where both samples
 are heterozygous
 The top-right half of the matrix contains counts of sites where both samples are homozygous and different from each other (i.e. 1/1::0/0 or 0/0::1/1)
 OUTPUT:
 ratio count(both heterozygous)/count(hom different) in bottom-left half-matrix
 */
void finalize_diffs_Hets_vs_Homs_proportions(std::vector<std::vector<double> >& diffs_Hets_vs_Homs){
    for (std::vector<std::vector<int> >::size_type i = 0; i < diffs_Hets_vs_Homs.size(); i++) {
        for (int j = 0; j < i; j++) {
            diffs_Hets_vs_Homs[i][j] = (diffs_Hets_vs_Homs[i][j]/diffs_Hets_vs_Homs[j][i]);
            diffs_Hets_vs_Homs[j][i] = 0;
        }    
    }
}


