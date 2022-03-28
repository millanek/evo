//
//  process_vcf_utils.h
//  C++ project
//
//  Created by Milan Malinsky on 08/05/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#ifndef C___project_process_vcf_utils_h
#define C___project_process_vcf_utils_h
#include <iostream>
#include <map>
#include <vector>
#include <sstream>
#include <fstream>
#include <limits>
#include <assert.h>
#include <algorithm>
#include <getopt.h>
#include <cstdlib>
#include <stdexcept>
#include <math.h>
#include <random>
#include <unordered_map>
#include <ctime>
#include "gzstream.h"
using std::string;
#define PROGRAM_BIN "evo"
#define PACKAGE_BUGREPORT "mm21@sanger.ac.uk"
#define GZIP_EXT ".gz"

#define LikelihoodsProbabilitiesAbsent 0
#define LikelihoodsProbabilitiesGP 1
#define LikelihoodsProbabilitiesGL 2
#define LikelihoodsProbabilitiesPL 3

#define AncestralAlleleMissing -1
#define AncestralAlleleRef 0
#define AncestralAlleleAlt 1


// VCF format constant
static const int NUM_NON_GENOTYPE_COLUMNS=9;  // 8 mendatory columns + 1 column with definition of the genotype columns

void printMissingLikelihoodsWarning(const string& chr, const string& pos);
void splitToDouble(const std::string &s, char delim, std::vector<double> &elems);
std::vector<double> splitToDouble(const std::string &s, char delim);
std::vector<std::string> split(const std::string &s, char delim);
// Find which fields in the VCF file corresponds to samples listed in a given set
std::vector<size_t> locateSet(std::vector<std::string>& sample_names, const std::vector<std::string>& set);
size_t locateOneSample(std::vector<std::string>& sample_names, const std::string toFind);

// Converting numbers (int, double, size_t, and char) to string
template <typename T> std::string numToString(T i) {
    std::string ret;
    std::stringstream out;
    out << i;
    ret = out.str();
    return ret;
}


template <typename T> void reset_matrix_to_zero(std::vector<std::vector<T> >& m) {
    for (int i = 0; i < m.size(); i++) {
        std::fill(m[i].begin(), m[i].end(), 0);
    }
}

class GeneralSetCounts {
public:
    GeneralSetCounts(const std::map<string, std::vector<size_t>>& setsToPosMap, const int nSamples) : overall(0), likelihoodsProbabilitiesType(LikelihoodsProbabilitiesAbsent) {
        for(std::map<string, std::vector<size_t>>::const_iterator it = setsToPosMap.begin(); it != setsToPosMap.end(); ++it) {
            setRefCounts[it->first] = 0; setAltCounts[it->first] = 0; setAlleleCounts[it->first] = 0;
            setAlleleProbCounts[it->first] = 0;
            setAAFs[it->first] = -1.0; setDAFs[it->first] = -1.0;
            setAAFsFromLikelihoods[it->first] = -1.0; setDAFsFromLikelihoods[it->first] = -1.0;
            setSizes.push_back(it->second.size());
            setHWEpriorsFromAAFfromGT[it->first].assign(3, -1.0);
            setHWEpriorsFromDAFfromGT[it->first].assign(3, -1.0);
        }
        individualsWithVariant.assign(nSamples, -1);
    };
    
    void getSetVariantCountsSimple(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap);
    void getSetVariantCounts(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap);
    void getAFsFromGenotypeLikelihoodsOrProbabilities(const std::vector<std::string>& genotypeFields, const std::map<size_t, string>& posToSpeciesMap, const int likelihoodsOrProbabilitiesTagPosition);
    int checkForGenotypeLikelihoodsOrProbabilities(const std::vector<std::string>& vcfLineFields);
    int returnFormatTagPosition(std::vector<std::string>& format, const std::string& tag);
    void setHWEpriorsFromAFfromGT();
    std::vector<double> probabilitiesFromLikelihoods(const std::vector<double>& thisLikelihoods, const string& species);
    
    int overall; int AAint;
    std::map<string,int> setRefCounts;
    std::map<string,int> setAltCounts;
    std::map<string,int> setAlleleCounts; // The number of non-missing alleles for this set
    std::map<string,int> setAlleleProbCounts; // The number of non-missing alleles for this set in terms of likelihoods/probabilities
    std::vector<size_t> setSizes;
    std::map<string,double> setAAFs; double averageAAF; // Allele frequencies - alternative allele
    std::map<string,double> setDAFs; double averageDAF;// Allele frequencies - derived allele
    std::map<string,double> setAAFsFromLikelihoods; double averageAAFFromLikelihoods; // Allele frequencies - alternative allele
    std::map<string,double> setDAFsFromLikelihoods; double averageDAFFromLikelihoods;// Allele frequencies - derived allele
    std::vector<int> individualsWithVariant; // 0 homRef, 1 het, 2 homAlt
    int likelihoodsProbabilitiesType;
    std::map<string,std::vector<double> > setHWEpriorsFromAAFfromGT;
    std::map<string,std::vector<double> > setHWEpriorsFromDAFfromGT;
    // std::vector<int> set1individualsWithVariant; std::vector<int> set2individualsWithVariant;
    // std::vector<int> set3individualsWithVariant; std::vector<int> set4individualsWithVariant;
    
private:
    void getBasicCounts(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap);
};

class GeneralSetCountsWithComplements : public GeneralSetCounts {
    public:
    GeneralSetCountsWithComplements(const std::map<string, std::vector<size_t>>& setsToPosMap, const int nSamples) : GeneralSetCounts(setsToPosMap,nSamples) {
        for(std::map<string, std::vector<size_t>>::const_iterator it = setsToPosMap.begin(); it != setsToPosMap.end(); ++it) {
            setAAFsComplement[it->first] = -1.0; setDAFsComplement[it->first] = -1.0; setAlleleCountsComplement[it->first] = 0;
        }
    }
    std::map<string,double> setAAFsComplement; // Allele frequencies - alternative allele, in the complement of the set
    std::map<string,double> setDAFsComplement; // Allele frequencies - derived allele, in the complement of the set
    std::map<string,int> setAlleleCountsComplement; // The number of non-missing alleles for the complement of this set
    
    void getComplementCounts(const std::vector<string>& populationsToUse);
};



class Counts {
public:
    Counts() : overall(0), minimumDepthInAnIndividual(std::numeric_limits<int>::max()), overallDepth(0), inbreedingCoefficient(0), chiSqPvalForInbreeding(1), bPhased(false), bIsMultiallelic(false), bAnyMissingGenotypes(false), SGB(std::numeric_limits<double>::max()) {};
    
    int overall;
    int minimumDepthInAnIndividual;
    int overallDepth;
    std::string FSpval; // Phred-scaled pvalue
    std::string MQSBpval; // Phred-scaled pvalue
    double SGB; // 
    double inbreedingCoefficient;
    double chiSqPvalForInbreeding;
    bool bPhased;
    std::vector<int> individualsWithVariant;
    std::vector<int> depthPerIndividual;
    std::vector<int> genotypeQualitiesPerIndividual;
    std::vector<int> haplotypesWithVariant;
    std::vector<bool> missingGenotypesPerIndividual;
    bool bAnyMissingGenotypes;
    bool bIsMultiallelic;
};


class MultiallelicCounts {
public:
    MultiallelicCounts(int nSamples, int starPosInput) : inbreedingCoefficient(0), chiSqPvalForInbreeding(1), bPhased(false) {
        individualsHets.assign(nSamples, -1);
        missingIndividualsDot.assign(nSamples, 0);
        missingIndividualsAny.assign(nSamples, 0);
        missingHaplotypesDot.assign(2*nSamples, 0);
        
        missingHaplotypesStar.assign(2*nSamples, 0);
        haplotypeVariants.assign(2*nSamples, -1);
        starPos = starPosInput;
    };
    
    int starPos;
    double inbreedingCoefficient;
    double chiSqPvalForInbreeding;
    bool bPhased;
    std::vector<int> individualsHets;
    std::vector<int> haplotypeVariants;
    
    std::vector<int> missingIndividualsDot;  // Individuals which are denoted in the VCF as missing data ./.
    std::vector<int> missingIndividualsAny;  // Individuals with missing data (as above) or which have either haplotype denoted in the VCF as * alleles (missing at the site because of an indel nearby)
    
    std::vector<int> missingHaplotypesDot;  // Haplotypes which are denoted in the VCF as missing data ./.
    std::vector<int> missingHaplotypesStar;   // Haplotypes which are denoted in the VCF as * alleles (missing at the site because of an indel nearby)
    
    void getMultiallelicCounts(const std::vector<string>& genotypes);
    double getPiThisVariant();
    double getHeterozygosityThisVariant();
    
};


class SingleSetCounts {
public:
    SingleSetCounts() : overall(0), setCount(0) {};
    
    int overall;
    int setCount;
    std::vector<int> individualsWithVariant;
    std::vector<int> setindividualsWithVariant;
};

class ThreeSetCounts {
public:
    ThreeSetCounts() : overall(0), set1AltCount(0), set2AltCount(0), set3AltCount(0), set1RefCount(0), set2RefCount(0), set3RefCount(0), set1AltAF(-1), set2AltAF(-1), set3AltAF(-1), set1daAF(-1), set2daAF(-1), set3daAF(-1), bAnyMissingGenotypes(false),  bIndel(false), n_alt_alleles(0), set1_n_withoutMissing(0), set2_n_withoutMissing(0), set3_n_withoutMissing(0) {};
    
    
    int overall; bool bAnyMissingGenotypes; bool bIndel;
    int set1AltCount; int set2AltCount; int set3AltCount;
    int set1RefCount; int set2RefCount; int set3RefCount;
    double set1AltAF; double set2AltAF; double set3AltAF;  // Allele frequencies - alternative allele
    double set1daAF; double set2daAF; double set3daAF;  // Allele frequencies - derived allele
    std::vector<int> individualsWithVariant;
    std::vector<int> set1individualsWithVariant; std::vector<int> set2individualsWithVariant;
    std::vector<int> set3individualsWithVariant;
    std::vector<bool> missingGenotypesPerIndividual;
    int n_alt_alleles; // the number of alternative alleles at a locus
    int set1_n_withoutMissing;
    int set2_n_withoutMissing;
    int set3_n_withoutMissing;
    std::vector<int> haplotypesWithVariant;
    std::vector<int> set1HaplotypeVariant; std::vector<int> set2HaplotypeVariant; std::vector<int> set3HaplotypeVariant;
};



class SetCountGeneric {
public:
    SetCountGeneric() : setRefCount(0), setAltCount(0), setAltAF(0), setDaAF(0) {};
    
    int setRefCount; int setAltCount;
    double setAltAF; double setDaAF; 
    std::vector<int> setindividualsWithVariant;
    
    bool isPolymorhic() {
        if (setAltAF != 0 && setAltAF != 1)
            return true;
        else
            return false;
    }
};

class FourSetCounts {
public:
    FourSetCounts() : overall(0), set1AltCount(0), set2AltCount(0), set3AltCount(0), set4AltCount(0), set1RefCount(0), set2RefCount(0), set3RefCount(0), set4RefCount(0), set1AltAF(-1), set2AltAF(-1), set3AltAF(-1), set4AltAF(-1), set1daAF(-1), set2daAF(-1), set3daAF(-1), set4daAF(-1) {};
    
    int overall;
    int set1AltCount; int set2AltCount; int set3AltCount; int set4AltCount;
    int set1RefCount; int set2RefCount; int set3RefCount; int set4RefCount;
    double set1AltAF; double set2AltAF; double set3AltAF; double set4AltAF; // Allele frequencies - alternative allele
    double set1daAF; double set2daAF; double set3daAF; double set4daAF; // Allele frequencies - derived allele
    std::vector<int> individualsWithVariant;
    std::vector<int> set1individualsWithVariant; std::vector<int> set2individualsWithVariant;
    std::vector<int> set3individualsWithVariant; std::vector<int> set4individualsWithVariant;
};

class SetCounts {
public:
    SetCounts() : overall(0), set1Count(0), set2Count(0), fisher_pval(1), chi_sq_pval(1),set1_n_withoutMissing(0), set2_n_withoutMissing(0), bAnyMissingGenotypes(false), n_alt_alleles(0), bIndel(false), p1(0), p2(0) {};
    
    int overall;
    bool bAnyMissingGenotypes;
    std::vector<int> individualsWithVariant;
    std::vector<bool> missingGenotypesPerIndividual;
    std::vector<int> haplotypesWithVariant;
    std::vector<int> set1HaplotypeVariant;
    std::vector<int> set2HaplotypeVariant;
    std::vector<int> set1individualsWithVariant;
    std::vector<int> set2individualsWithVariant;
    int set1Count;
    int set2Count;
    int set1_n_withoutMissing;
    int set2_n_withoutMissing;
    int n_alt_alleles; // the number of alternative alleles at a locus
    bool bIndel;
    double fisher_pval;
    double chi_sq_pval;
    double p1; double p2; // allele frequencies
    
    void calculateAlleleFrequencies() {
        p1 = (double)set1Count/set1_n_withoutMissing;
        p2 = (double)set2Count/set2_n_withoutMissing;
    }
    
    void reset() {
        overall = 0;
        set1Count = 0;
        set2Count = 0;
        fisher_pval = 0;
        chi_sq_pval = 0;
    }
};


// Results object
class FilterResult
{
public:
    FilterResult() : overallQuality(0), overallQualityPassed(true), overallDepthPassed(true), strandBiasPassed(true), inbreedingCoeffPassed(true), counts(), biallelicPassed(true), siteInvariant(true) {};
    
    int overallQuality;  // -10log_10 p(no variant)
    bool overallQualityPassed;
    bool overallDepthPassed;
    bool strandBiasPassed;
    bool inbreedingCoeffPassed;
    Counts counts;
    bool biallelicPassed;
    bool siteInvariant;
};

// For checking if variants fixed in Massoko are also fixed in Malawi
class MassokoMalawiResult {
public:
    void init(std::vector<int>::size_type size) {
        absentWithBlue.assign(size, 0);
        hetsWithBlue.assign(size, 0);
        homsWithBlue.assign(size, 0);
        absentWithYellow.assign(size, 0);
        hetsWithYellow.assign(size, 0);
        homsWithYellow.assign(size, 0);
    };
    
    std::vector<int> absentWithBlue;
    std::vector<int> hetsWithBlue;
    std::vector<int> homsWithBlue;
    std::vector<int> absentWithYellow;
    std::vector<int> hetsWithYellow;
    std::vector<int> homsWithYellow;
};



// Print an arbitrary matrix (vector of vectors)
template <class T> void print_matrix(T matrix, std::ostream& outFile) {
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            if (j == (matrix[i].size()-1))
                outFile << matrix[i][j] << std::endl;
            else 
                outFile << matrix[i][j] << "\t";
        }    
    }    
}

// Print an arbitrary vector to a file
template <class T> void print_vector(T vector, std::ostream& outFile, char delim = '\t', bool endLine = true) {
    for (int i = 0; i < vector.size(); i++) {
        if (i == (vector.size()-1)) {
            if (endLine) outFile << vector[i] << std::endl;
            else outFile << vector[i];
        } else {
            outFile << vector[i] << delim;
        }
    }
}


template <class T> double vector_sum(T vector) {
    double sum = 0;
    for (int i = 0; i < vector.size(); i++) {
        sum += vector[i];
    }
    return sum;
}

template<typename T>
std::vector<T> slice(std::vector<T> const &v, int m, int n)
{
    auto first = v.cbegin() + m;
    auto last = v.cbegin() + n + 1;
    
    std::vector<T> vec(first, last);
    return vec;
}


template <class T> double vector_average(T vector) {
    double sum = 0;
    for (int i = 0; i < vector.size(); i++) {
        sum += vector[i];
    }
    double average = (double)sum / (double)vector.size();
    return average;
}

template <class T> double vector_average_withRegion(T vector, int regionLength) {
    double sum = 0;
    for (int i = 0; i < vector.size(); i++) {
        sum += vector[i];
    }
    double average = (double)sum / (double)regionLength;
    return average;
}


template <typename T> std::map<int, int> tabulateVectorTemplate(T& vec) {
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

class BadConversion : public std::runtime_error {
public:
    BadConversion(std::string const& s)
    : std::runtime_error(s)
    { }
};

inline double convertToDouble(std::string const& s)
{
    std::istringstream i(s);
    double x;
    if (!(i >> x))
        throw BadConversion("convertToDouble(\"" + s + "\")");
    return x;
}

inline bool file_exists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

void getSetVariantCounts(GeneralSetCounts* c, const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap);
void getSetVariantCountsSimple(GeneralSetCounts* c, const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap);

double stringToDouble(std::string s);

double calculateInbreedingCoefficient(std::vector<int>& individualsWithVariant);
double calculateChiSqPvalForInbreeding(std::vector<int>& individualsWithVariant);

// Does the same as R function table
std::map<int, int> tabulateVector(std::vector<int>& vec);

void initialize_matrix_double(std::vector<std::vector<double> >& m, int m_rows, int m_columns = 0);
void initialize_matrix_int(std::vector<std::vector<int> >& m, int m_rows, int m_columns = 0);

std::vector<size_t> complementIndices(const size_t fullVectorSize, const std::vector<size_t>& originalIndices);

// Remove a single file extension from the filename
std::string stripExtension(const std::string& filename);
std::string stripPath(const std::string& filename);

Counts getThisVariantCounts(const std::vector<std::string>& fields);
void getThisVariantCountsSimple(const std::vector<std::string>& fields, Counts* thisVariantCounts);
ThreeSetCounts getThreeSetVariantCounts(const std::vector<std::string>& fields, const std::vector<size_t>& set1_loci, const std::vector<size_t>& set2_loci, const std::vector<size_t>& set3_loci, const std::string& AA);
FourSetCounts getFourSetVariantCounts(const std::vector<std::string>& fields, const std::vector<size_t>& set1_loci, const std::vector<size_t>& set2_loci, const std::vector<size_t>& set3_loci, const std::vector<size_t>& set4_loci, const std::string& AA = "N");
ThreeSetCounts getThreeSetVariantCountsAA4(const std::vector<std::string>& fields, const std::vector<size_t>& set1_loci, const std::vector<size_t>& set2_loci, const std::vector<size_t>& set3_loci, const std::vector<size_t>& AA_loci);


bool testBiallelic(const std::string& altField);

bool testOverallReadDepth(const int maxReadDepth, const int minReadDepth, const std::string& infoField);

// filter out sites where more than maxNumHet individuals are heterozygous 
bool testMaxNumHet(FilterResult& result, std::vector<int>& depthsHetFailed, std::vector<int>& depthsHetPassed, std::vector<int>& numVariantsPerHetCount, int maxNumHet, std::vector<std::vector<int> >& num_indiv_het_vs_depth);

//void doubleton_analysis(std::vector<std::vector<int> >& doubletons, FilterResult& result, int numChromosomes);
// Move all doubleton counts to the bottom left corner of the matrix
void rearrange_doubletons(std::vector<std::vector<int> >& doubletons);

// For checking if variants fixed in Massoko are also fixed in Malawi
void massokoMalawiSharing(const FilterResult& result, MassokoMalawiResult& sharingResult);

// Read the sample names file (one sample name per line)
std::vector<std::string> readSampleNamesFromTextFile(const std::string& sampleNameFile);
// Open a file that may or may not be gzipped for reading - The caller is responsible for freeing the handle
void assertGZOpen(gzstreambase& gh, const std::string& fn);
void assertFileOpen(std::ifstream& fh, const std::string& fn);
void assertFileOpen(std::ofstream& fh, const std::string& fn);
std::istream* createReader(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in);
std::ostream* createWriter(const std::string& filename, std::ios_base::openmode mode = std::ios_base::out);
void print80bpPerLineStdOut(std::ostream& outFile, std::string toPrint);
void print80bpPerLineFile(std::ofstream*& outFile, string toPrint);


// Reading in genome files
std::map<std::string,std::string> readMultiFastaToMap(const std::string& fileName);
std::string readMultiFastaToOneString(const string& fileName, int bytes = 50000000);

/*  INPUT MATRIX:
    In the bottom-left part of the matrix (addressed here as diffs_Hets_vs_Homs[i][j]) should be counts of sites where both samples
    are heterozygous
    The top-right half of the matrix contains counts of sites where both samples are homozygous and different from each other (i.e. 1/1::0/0 or 0/0::1/1)
    OUTPUT:
    ratio count(both heterozygous)/count(hom different) in bottom-left half-matrix */
void finalize_diffs_Hets_vs_Homs_proportions(std::vector<std::vector<double> >& diffs_Hets_vs_Homs);


#endif
