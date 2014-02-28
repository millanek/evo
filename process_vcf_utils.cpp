//
//  process_vcf_utils.cpp
//  C++ project
//
//  Created by Milan Malinsky on 08/05/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include <iostream>
#include "zreaderDaniel.h"
#include "process_vcf_utils.h"


void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

// Initialize a matrix 
void initialize_matrix_double(std::vector<std::vector<double> >& m, int m_size) {
    for (int i = 0; i < m_size; i++) { 
        std::vector<double> v(m_size,0);
        m.push_back(v);
    }       
}

// Initialize a matrix 
void initialize_matrix_int(std::vector<std::vector<int> >& m, int m_size) {
    for (int i = 0; i < m_size; i++) { 
        std::vector<int> v(m_size,0);
        m.push_back(v);
    }       
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

Counts getThisVariantCounts(const std::vector<std::string>& fields) {
    Counts thisVariantCounts;
    thisVariantCounts.individualsWithVariant.assign((fields.size()-NUM_NON_GENOTYPE_COLUMNS),0);
    //std::cerr << "Fields: " << (fields.size()-NUM_NON_GENOTYPE_COLUMNS) << std::endl;
    for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
        if (fields[i][0] == '1') { 
            thisVariantCounts.overall++;
            thisVariantCounts.individualsWithVariant[i- NUM_NON_GENOTYPE_COLUMNS]++;
        }
        if (fields[i][2] == '1') { 
            thisVariantCounts.overall++;
            thisVariantCounts.individualsWithVariant[i-NUM_NON_GENOTYPE_COLUMNS]++;
        }
        std::vector<std::string> genotypeData = split(fields[i], ':');
        
        if (atoi(genotypeData[2].c_str()) < thisVariantCounts.minimumDepthInAnIndividual) {
            thisVariantCounts.minimumDepthInAnIndividual = atoi(genotypeData[2].c_str());
        }
        
    }
    // Also get overall depth for this variant
    std::vector<std::string> info = split(fields[7], ';');
    if (info[0] == "INDEL") {
        split(info[1], '=', info);
    } else {
        split(info[0], '=', info);
    }
    thisVariantCounts.overallDepth = atoi((info.back()).c_str());
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



bool testBiallelic(const std::string& altField) {
    std::vector<std::string> altVector = split(altField, ',');
    if (altVector.size() == 1) { return true; }
    else { return false; }
}

bool testOverallReadDepth(const int maxReadDepth, const std::string& infoField) {
    std::vector<std::string> info = split(infoField, ';');
    if (info[0] == "INDEL") {
        split(info[1], '=', info);
    } else {
        split(info[0], '=', info);
    }
    int DP = atoi((info.back()).c_str());
    if (DP <= maxReadDepth) { return true; }
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

// Open a file that may or may not be gzipped for reading
// The caller is responsible for freeing the handle
std::istream* createReader(const std::string& filename)
{
    if(isGzip(filename))
    {
        dpj::zifstream* pGZ = new dpj::zifstream(filename);
        // assertGZOpen(*pGZ, filename);
        return pGZ;
    }
    else
    {
        std::ifstream* pReader = new std::ifstream(filename.c_str());
        assertFileOpen(*pReader, filename);
        return pReader;
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
            std::cerr << "Did not find the sample: " << set[i] << std::endl;
        } else {
            size_t loc = std::distance(sample_names.begin(), it);
            setLocs.push_back(loc);
        }
    }
    return setLocs;
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
    string currentScaffold = line.substr(1,string::npos);
    fastaSeqs[currentScaffold] = ""; fastaSeqs[currentScaffold].reserve(50000000);
    while (getline(*fastaFile, line)) {
        if (line[0] != '>') {
            fastaSeqs[currentScaffold].append(line);
        } else {
            // std::cerr << currentScaffold << " length: " << ancSeqs[currentScaffold].length() << std::endl;
            currentScaffold = line.substr(1,string::npos);
            fastaSeqs[currentScaffold] = ""; fastaSeqs[currentScaffold].reserve(50000000);
        }
    }
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


