//
//  process_vcf_coding_sequences.h
//  vcf_process
//
//  Created by Milan Malinsky on 16/09/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#ifndef __vcf_process__process_vcf_coding_sequences__
#define __vcf_process__process_vcf_coding_sequences__

#include "process_vcf_seq_utils.h"
#include "process_vcf_annotation_tools.h"
#include "process_vcf_stats_utils.h"
#include <unordered_map>
#include <unordered_set>


class pNsets {
public:
    pNsets() : set1vsSet2pN(0), sets1and2vsSet3pN(0), withinSet1andSet2pN(0), withinSet1pN(0), withinSet2pN(0), withinSet3pN(0), initialised(false) {};
    
    pNsets(std::ifstream*& setsFile) {
        string line;
        string set1String; string set2String; string set3String;
        getline(*setsFile, set1String);
        getline(*setsFile, set2String);
        getline(*setsFile, set3String);
        std::vector<string> set1 = split(set1String, ',');
        std::vector<string> set2 = split(set2String, ',');
        std::vector<string> set3 = split(set3String, ',');
        for (int i = 0; i < set1.size(); i++) {
            set1Loci.insert(atoi(set1[i].c_str()));
        }
        for (int i = 0; i < set2.size(); i++) {
            set2Loci.insert(atoi(set2[i].c_str()));
        }
        for (int i = 0; i < set3.size(); i++) {
            set3Loci.insert(atoi(set3[i].c_str()));
        }
        std::cerr << "set 1: "; print_vector_stream(set1, std::cerr);
        std::cerr << "set 2: "; print_vector_stream(set2, std::cerr);
        std::cerr << "set 3: "; print_vector_stream(set3, std::cerr);
        initialised = true;
    }
    
    std::unordered_set<size_t> set1Loci;
    std::unordered_set<size_t> set2Loci;
    std::unordered_set<size_t> set3Loci;
    
    double withinSet1pN;
    double withinSet2pN;
    double withinSet3pN;
    double withinSet1andSet2pN;
    double set1vsSet2pN;
    double sets1and2vsSet3pN;
    bool initialised = false;
    
};


class CDSComparisonMatrices {
public:
    CDSComparisonMatrices(int numSamples) {
        initialize_matrix_double(N_d_jk, numSamples); initialize_matrix_double(N_jk, numSamples);
        initialize_matrix_double(S_d_jk, numSamples); initialize_matrix_double(S_jk, numSamples);
        initialize_matrix_double(tS_N_jk, numSamples); initialize_matrix_double(tS_N_jk, numSamples);
        initialize_matrix_double(tS_S_jk, numSamples); initialize_matrix_double(tS_S_jk, numSamples);
        initialize_matrix_double(tV_N_jk, numSamples); initialize_matrix_double(tV_N_jk, numSamples);
        initialize_matrix_double(tV_S_jk, numSamples); initialize_matrix_double(tV_S_jk, numSamples);
    }
    std::vector<std::vector<double> > N_d_jk; std::vector<std::vector<double> > N_jk;
    std::vector<std::vector<double> > S_d_jk; std::vector<std::vector<double> > S_jk;
    // Separate matrices to incorporate unequal tS/tV mutation probabilities:
    std::vector<std::vector<double> > tS_N_jk; std::vector<std::vector<double> > tS_S_jk;
    std::vector<std::vector<double> > tV_N_jk; std::vector<std::vector<double> > tV_S_jk;
};


class CDSH1H2ComparisonMatrices {
public:
    CDSH1H2ComparisonMatrices(int numSamples) {
        H1p = new CDSComparisonMatrices(numSamples);
        H2p = new CDSComparisonMatrices(numSamples);
        H1H2p = new CDSComparisonMatrices(numSamples);
    }
    
    CDSComparisonMatrices* H1p;
    CDSComparisonMatrices* H2p;
    CDSComparisonMatrices* H1H2p;
};


enum t_codon {
    TTT, TTC,
    TTA, TTG,
    CTA, CTC, CTG, CTT,
    ATA, ATC, ATT,
    ATG,
    GTA, GTC, GTG, GTT,
    TCA, TCC, TCG, TCT,
    AGC, AGT,
    CCA, CCC, CCG, CCT,
    ACA, ACC, ACG, ACT,
    GCA, GCC, GCG, GCT,
    TAC, TAT,
    TAA, TAG, TGA,
    CAC, CAT,
    CAA, CAG,
    AAC, AAT,
    AAA, AAG,
    GAC, GAT,
    GAA, GAG,
    TGC, TGT,
    TGG,
    CGA, CGC, CGG, CGT,
    AGA,AGG,
    GGA, GGC, GGG, GGT
};


static std::unordered_map<std::string, t_codon> s_mapStringToTcodon =
{
    { "TTT", TTT }, { "TTC", TTC },
    { "TTA", TTA }, { "TTG", TTG },
    { "CTA", CTA }, { "CTC", CTC }, { "CTG", CTG }, { "CTT", CTT },
    { "ATA", ATA }, { "ATC", ATC }, { "ATT", ATT },
    { "ATG", ATG },
    { "GTA", GTA }, { "GTC", GTC }, { "GTG", GTG }, { "GTT", GTT },
    { "TCA", TCA }, { "TCC", TCC }, { "TCG", TCG }, { "TCT", TCT },
    { "AGC", AGC }, { "AGT", AGT },
    { "CCA", CCA }, { "CCC", CCC }, { "CCG", CCG }, { "CCT", CCT },
    { "ACA", ACA }, { "ACC", ACC }, { "ACG", ACG }, { "ACT", ACT },
    { "GCA", GCA }, { "GCC", GCC }, { "GCG", GCG }, { "GCT", GCT },
    { "TAC", TAC }, { "TAT", TAT },
    { "TAA", TAA }, { "TAG", TAG }, { "TGA", TGA },
    { "CAC", CAC }, { "CAT", CAT },
    { "CAA", CAA }, { "CAG", CAG },
    { "AAT", AAT }, { "AAC", AAC },
    { "AAA", AAA }, { "AAG", AAG },
    { "GAC", GAC }, { "GAT", GAT },
    { "GAA", GAA }, { "GAG", GAG },
    { "TGC", TGC }, { "TGT", TGT },
    { "TGG", TGG },
    { "CGA", CGA }, { "CGC", CGC }, { "CGG", CGG }, { "CGT", CGT },
    { "AGA", AGA }, { "AGG", AGG },
    { "GGA", GGA }, { "GGC", GGC }, { "GGG", GGG }, { "GGT", GGT }
};

static std::unordered_map<std::string, double> s_mapCodonToExpDist =
{
    { "TTT", 8.0/3.0 }, { "TTC", 8.0/3.0 },
    { "TTA", 7.0/3.0 }, { "TTG", 7.0/3.0 },                                             // Leu
    { "CTA", 5.0/3.0 }, { "CTG", 5.0/3.0 }, { "CTC", 2.0 }, { "CTT", 2.0 },             // Leu
    { "ATA", 7.0/3.0 }, { "ATC", 7.0/3.0 }, { "ATT", 7.0/3.0 },
    { "ATG", 3.0 },
    { "GTA", 2.0 }, { "GTC", 2.0 }, { "GTG", 2.0 }, { "GTT", 2.0 },
    { "TCA", 2.0 }, { "TCC", 2.0 }, { "TCG", 2.0 }, { "TCT", 2.0 },                     // Ser
    { "AGC", 8.0/3.0 }, { "AGT", 8.0/3.0 },                                             // Ser
    { "CCA", 2.0 }, { "CCC", 2.0 }, { "CCG", 2.0 }, { "CCT", 2.0 },
    { "ACA", 2.0 }, { "ACC", 2.0 }, { "ACG", 2.0 }, { "ACT", 2.0 },
    { "GCA", 2.0 }, { "GCC", 2.0 }, { "GCG", 2.0 }, { "GCT", 2.0 },
    { "TAC", 8.0/3.0 }, { "TAT", 8.0/3.0 },
    { "TAA", 7.0/3.0 }, { "TAG", 8.0/3.0 }, { "TGA", 8.0/3.0 },                         // Stop
    { "CAC", 8.0/3.0 }, { "CAT", 8.0/3.0 },
    { "CAA", 8.0/3.0 }, { "CAG", 8.0/3.0 },
    { "AAT", 8.0/3.0 }, { "AAC", 8.0/3.0 },
    { "AAA", 8.0/3.0 }, { "AAG", 8.0/3.0 },
    { "GAC", 8.0/3.0 }, { "GAT", 8.0/3.0 },
    { "GAA", 8.0/3.0 }, { "GAG", 8.0/3.0 },
    { "TGC", 8.0/3.0 }, { "TGT", 8.0/3.0 },
    { "TGG", 3.0 },
    { "CGA", 5.0/3.0 }, { "CGG", 5.0/3.0 }, { "CGC", 2.0 }, { "CGT", 2.0 },             // Arg
    { "AGA", 7.0/3.0 }, { "AGG", 7.0/3.0 },                                             // Arg
    { "GGA", 2.0 }, { "GGC", 2.0 }, { "GGG", 2.0 }, { "GGT", 2.0 }
};





// Process the fifth column from the annatation file to get this gene name
// and information about whether this gene is 'partial' (passed on through the bool parameter)
std::string getGeneName(const std::string& geneColumn, bool& bPartial);

// Various QC checks to verify that we really have a protein coding sequence (length, start/stop codons)
// Return false if the length of the coding sequence is not divisible by three
bool codingSequenceErrorChecks(const std::string& geneSeq, const std::string& transcriptName, const std::vector<std::vector<std::string> >& annotation, const int k, std::ofstream*& badStartStopCodonFile);
// Calculate some statistics about the sequences
void getStatsHaploidSeq(const std::vector<std::string>& allSeqs, std::vector<std::string>& statsThisGene, double tStVratio = 0.5);
void getStatsBothPhasedHaps(const std::vector<std::string>& allSeqs, const std::vector<std::string>& allSeqsH2, std::vector<string>& statsThisGene, std::vector<std::vector<double> >& combinedVectorForPCA, pNsets* sets, double tStVratio = 0.5, bool nonCodingNull = false);
void getStatsIUPAC(const std::vector<std::string>& allSeqs, const std::string& refSeq, const std::string& transcriptName, std::vector<std::string>& statsThisGene, std::ofstream*& prematureStopCodonFile, const std::vector<std::string>& sampleNames);

void parseGetCodingSeqOptions(int argc, char** argv);

int getCodingSeqMain(int argc, char** argv);

inline std::string getAminoAcid(const std::string& nuc)
{
    if (nuc == "TTT" || nuc == "TTC") 
        return "Phe";
    if (nuc == "TTA" || nuc == "TTG" || nuc == "CTA" || nuc == "CTC" || nuc == "CTG" || nuc == "CTT")
        return "Leu";
    if (nuc == "ATA" || nuc == "ATC" || nuc == "ATT")
        return "Ile";
    if (nuc == "ATG")
        return "Met";
    if (nuc == "GTA" || nuc == "GTC" || nuc == "GTG" || nuc == "GTT")
        return "Val";
    if (nuc == "TCA" || nuc == "TCC" || nuc == "TCG" || nuc == "TCT")
        return "Ser";
    if (nuc == "CCA" || nuc == "CCC" || nuc == "CCG" || nuc == "CCT")
        return "Pro";
    if (nuc == "ACA" || nuc == "ACC" || nuc == "ACG" || nuc == "ACT")
        return "Thr";
    if (nuc == "GCA" || nuc == "GCC" || nuc == "GCG" || nuc == "GCT")
        return "Ala";
    if (nuc == "TAC" || nuc == "TAT")
        return "Tyr";
    if (nuc == "TAA" || nuc == "TAG" || nuc == "TGA")
        return "Stop";
    if (nuc == "CAC" || nuc == "CAT")
        return "His";
    if (nuc == "CAA" || nuc == "CAG")
        return "Gln";
    if (nuc == "AAC" || nuc == "AAT")
        return "Asn";
    if (nuc == "AAA" || nuc == "AAG")
        return "Lys";
    if (nuc == "GAC" || nuc == "GAT")
        return "Asp";
    if (nuc == "GAA" || nuc == "GAG")
        return "Glu";
    if (nuc == "TGC" || nuc == "TGT")
        return "Cys";
    if (nuc == "TGG")
        return "Trp";
    if (nuc == "CGA" || nuc == "CGC" || nuc == "CGG" || nuc == "CGT")
        return "Arg";
    if (nuc == "AGC" || nuc == "AGT")
        return "Ser";
    if (nuc == "AGA" || nuc == "AGG")
        return "Arg";
    if (nuc == "GGA" || nuc == "GGC" || nuc == "GGG" || nuc == "GGT")
        return "Gly";
    else
        return "Uncrecognised codon....";
    
}


// Stop codons are 'X'
// And unrecodnised sequences are 'Z'
inline char getAminoAcidOneLetterCode(const std::string& nuc)
{
    assert(nuc.length() == 3);
    if (nuc == "TTT" || nuc == "TTC")
        return 'F';
    if (nuc == "TTA" || nuc == "TTG" || nuc == "CTA" || nuc == "CTC" || nuc == "CTG" || nuc == "CTT")
        return 'L';
    if (nuc == "ATA" || nuc == "ATC" || nuc == "ATT")
        return 'I';
    if (nuc == "ATG")
        return 'M';
    if (nuc == "GTA" || nuc == "GTC" || nuc == "GTG" || nuc == "GTT")
        return 'V';
    if (nuc == "TCA" || nuc == "TCC" || nuc == "TCG" || nuc == "TCT")
        return 'S';
    if (nuc == "CCA" || nuc == "CCC" || nuc == "CCG" || nuc == "CCT")
        return 'P';
    if (nuc == "ACA" || nuc == "ACC" || nuc == "ACG" || nuc == "ACT")
        return 'T';
    if (nuc == "GCA" || nuc == "GCC" || nuc == "GCG" || nuc == "GCT")
        return 'A';
    if (nuc == "TAC" || nuc == "TAT")
        return 'Y';
    if (nuc == "TAA" || nuc == "TAG" || nuc == "TGA")
        return 'X';
    if (nuc == "CAC" || nuc == "CAT")
        return 'H';
    if (nuc == "CAA" || nuc == "CAG")
        return 'Q';
    if (nuc == "AAC" || nuc == "AAT")
        return 'N';
    if (nuc == "AAA" || nuc == "AAG")
        return 'K';
    if (nuc == "GAC" || nuc == "GAT")
        return 'D';
    if (nuc == "GAA" || nuc == "GAG")
        return 'E';
    if (nuc == "TGC" || nuc == "TGT")
        return 'C';
    if (nuc == "TGG")
        return 'W';
    if (nuc == "CGA" || nuc == "CGC" || nuc == "CGG" || nuc == "CGT")
        return 'R';
    if (nuc == "AGC" || nuc == "AGT")
        return 'S';
    if (nuc == "AGA" || nuc == "AGG")
        return 'R';
    if (nuc == "GGA" || nuc == "GGC" || nuc == "GGG" || nuc == "GGT")
        return 'G';
    else
        return 'Z';
    
}


inline double getExpectedNumberOfNonsynonymousSites(const std::string& nuc)
{
    switch (s_mapStringToTcodon[nuc])
    {
        case TTT: return 8.0/3.0; break;    case TTC: return 8.0/3.0; break;    // Phe
        case TTA: return 7.0/3.0; break;    case TTG: return 7.0/3.0; break;    // Leu
        case CTA: return 5.0/3.0; break;    case CTG: return 5.0/3.0; break;    // Leu
        case CTC: return 2.0; break;        case CTT: return 2.0; break;        // Leu
        case ATA: return 7.0/3.0; break;    case ATC: return 7.0/3.0; break;    // Ile
        case ATT: return 7.0/3.0; break;                                        // Ile
        case ATG: return 3.0; break;                                            // Met
        case GTA: return 2.0; break;        case GTC: return 2.0; break;        // Val
        case GTG: return 2.0; break;        case GTT: return 2.0; break;        // Val
        case TCA: return 2.0; break;        case TCC: return 2.0; break;        // Ser
        case TCG: return 2.0; break;        case TCT: return 2.0; break;        // Ser
        case AGC: return 8.0/3.0; break;    case AGT: return 8.0/3.0; break;    // Ser
        case CCA: return 2.0; break;        case CCC: return 2.0; break;        // Pro
        case CCG: return 2.0; break;        case CCT: return 2.0; break;        // Pro
        case ACA: return 2.0; break;        case ACC: return 2.0; break;        // Thr
        case ACG: return 2.0; break;        case ACT: return 2.0; break;        // Thr
        case GCA: return 2.0; break;        case GCC: return 2.0; break;        // Ala
        case GCG: return 2.0; break;        case GCT: return 2.0; break;        // Ala
        case GGA: return 2.0; break;        case GGC: return 2.0; break;        // Gly
        case GGG: return 2.0; break;        case GGT: return 2.0; break;        // Gly
        case TAC: return 8.0/3.0; break;    case TAT: return 8.0/3.0; break;    // Tyr
        case TAA: return 7.0/3.0; break;    case TAG: return 8.0/3.0; break;    // Stop
        case TGA: return 8.0/3.0; break;                                        // Stop
        case CAC: return 8.0/3.0; break;    case CAT: return 8.0/3.0; break;    // His
        case CAA: return 8.0/3.0; break;    case CAG: return 8.0/3.0; break;    // Gln
        case AAC: return 8.0/3.0; break;    case AAT: return 8.0/3.0; break;    // Asn
        case AAA: return 8.0/3.0; break;    case AAG: return 8.0/3.0; break;    // Lys
        case GAC: return 8.0/3.0; break;    case GAT: return 8.0/3.0; break;    // Asp
        case GAA: return 8.0/3.0; break;    case GAG: return 8.0/3.0; break;    // Glu
        case TGC: return 8.0/3.0; break;    case TGT: return 8.0/3.0; break;    // Cys
        case TGG: return 3.0; break;                                            // Trp
        case CGC: return 2.0; break;        case CGT: return 2.0; break;        // Arg
        case CGA: return 5.0/3.0; break;    case CGG: return 5.0/3.0; break;    // Arg
        case AGA: return 7.0/3.0; break;    case AGG: return 7.0/3.0; break;    // Arg
        }
}


inline double getExpectedNumberOfNonsynonymousSitesIf(const std::string& nuc)
{
    // Phe
    if (nuc == "TTT" || nuc == "TTC")
        return 8.0/3.0;
    // Leu
    if (nuc == "TTA" || nuc == "TTG")
        return 7.0/3.0;
    if (nuc == "CTA" || nuc == "CTG")
        return 5.0/3.0;
    if (nuc == "CTC" || nuc == "CTT")
        return 2.0;
    // Ile
    if (nuc == "ATA" || nuc == "ATC" || nuc == "ATT")
        return 7.0/3.0;
    // Met
    if (nuc == "ATG")
        return 3.0;
    // Val
    if (nuc == "GTA" || nuc == "GTC" || nuc == "GTG" || nuc == "GTT")
        return 2.0;
    // Ser
    if (nuc == "TCA" || nuc == "TCC" || nuc == "TCG" || nuc == "TCT")
        return 2.0;
    if (nuc == "AGC" || nuc == "AGT")
        return 8.0/3.0;
    // Pro
    if (nuc == "CCA" || nuc == "CCC" || nuc == "CCG" || nuc == "CCT")
        return 2.0;
    // Thr
    if (nuc == "ACA" || nuc == "ACC" || nuc == "ACG" || nuc == "ACT")
        return 2.0;
    // Ala
    if (nuc == "GCA" || nuc == "GCC" || nuc == "GCG" || nuc == "GCT")
        return 2.0;
    // Tyr
    if (nuc == "TAC" || nuc == "TAT")
        return 8.0/3.0;
    // Stop
    if (nuc == "TAA")
        return 7.0/3.0;
    if (nuc == "TGA" || nuc == "TAG")
        return 8.0/3.0;
    // His
    if (nuc == "CAC" || nuc == "CAT")
        return 8.0/3.0;
    // Gln
    if (nuc == "CAA" || nuc == "CAG")
        return 8.0/3.0;
    // Asn
    if (nuc == "AAC" || nuc == "AAT")
        return 8.0/3.0;
    // Lys
    if (nuc == "AAA" || nuc == "AAG")
        return 8.0/3.0;
    // Asp
    if (nuc == "GAC" || nuc == "GAT")
        return 8.0/3.0;
    // Glu
    if (nuc == "GAA" || nuc == "GAG")
        return 8.0/3.0;
    // Cys
    if (nuc == "TGC" || nuc == "TGT")
        return 8.0/3.0;
    // Trp
    if (nuc == "TGG")
        return 3.0;
    // Arg
    if (nuc == "CGA" || nuc == "CGG")
        return 5.0/3.0;
    if (nuc == "CGC" || nuc == "CGT")
        return 2.0;
    if (nuc == "AGA" || nuc == "AGG")
        return 7.0/3.0;
    // Gly
    if (nuc == "GGA" || nuc == "GGC" || nuc == "GGG" || nuc == "GGT")
        return 2.0;
    else
        return 0;
    
}


static std::unordered_map<std::string, double> s_mapCodonToExpTs =
{
    { "TTT", 2.0/3.0 }, { "TTC", 2.0/3.0 },
    { "TTA", 1.0/3.0 }, { "TTG", 1.0/3.0 },                                             // Leu
    { "CTA", 1.0/3.0 }, { "CTG", 1.0/3.0 }, { "CTC", 2.0/3.0 }, { "CTT", 2.0/3.0 },     // Leu
    { "ATA", 1.0 }, { "ATC", 2.0/3.0 }, { "ATT", 2.0/3.0 },                             // Ile
    { "ATG", 1.0 },                                                                     // Met
    { "GTA", 2.0/3.0 }, { "GTC", 2.0/3.0 }, { "GTG", 2.0/3.0 }, { "GTT", 2.0/3.0 },     // Val
    { "TCA", 2.0/3.0 }, { "TCC", 2.0/3.0 }, { "TCG", 2.0/3.0 }, { "TCT", 2.0/3.0 },     // Ser
    { "AGC", 2.0/3.0 }, { "AGT", 2.0/3.0 },                                             // Ser
    { "CCA", 2.0/3.0 }, { "CCC", 2.0/3.0 }, { "CCG", 2.0/3.0 }, { "CCT", 2.0/3.0 },     // Pro
    { "ACA", 2.0/3.0 }, { "ACC", 2.0/3.0 }, { "ACG", 2.0/3.0 }, { "ACT", 2.0/3.0 },     // Thr
    { "GCA", 2.0/3.0 }, { "GCC", 2.0/3.0 }, { "GCG", 2.0/3.0 }, { "GCT", 2.0/3.0 },     // Ala
    { "TAC", 2.0/3.0 }, { "TAT", 2.0/3.0 },                                             // Tyr
    { "TAA", 1.0/3.0 }, { "TAG", 2.0/3.0 }, { "TGA", 2.0/3.0 },                         // Stop
    { "CAC", 2.0/3.0 }, { "CAT", 2.0/3.0 },                                             // His
    { "CAA", 2.0/3.0 }, { "CAG", 2.0/3.0 },                                             // Gln
    { "AAT", 2.0/3.0 }, { "AAC", 2.0/3.0 },                                             // Asn
    { "AAA", 2.0/3.0 }, { "AAG", 2.0/3.0 },                                             // Lys
    { "GAC", 2.0/3.0 }, { "GAT", 2.0/3.0 },                                             // Asp
    { "GAA", 2.0/3.0 }, { "GAG", 2.0/3.0 },                                             // Glu
    { "TGC", 2.0/3.0 }, { "TGT", 2.0/3.0 },                                             // Cys
    { "TGG", 1.0 },                                                                     // Trp
    { "CGA", 2.0/3.0 }, { "CGG", 2.0/3.0 }, { "CGC", 2.0/3.0 }, { "CGT", 2.0/3.0 },     // Arg
    { "AGA", 2.0/3.0 }, { "AGG", 2.0/3.0 },                                             // Arg
    { "GGA", 2.0/3.0 }, { "GGC", 2.0/3.0 }, { "GGG", 2.0/3.0 }, { "GGT", 2.0/3.0 }      // Gly
};


// Expected number of synonymous transitions is always 1 - getExpectedNumberOfNonsynonymousTransitions()
inline double getExpectedNumberOfNonsynonymousTransitions(const std::string& nuc)
{
    // Phe
    if (nuc == "TTT" || nuc == "TTC")
        return 2.0/3.0;
    // Leu
    if (nuc == "TTA" || nuc == "TTG")
        return 1.0/3.0;
    if (nuc == "CTA" || nuc == "CTG")
        return 1.0/3.0;
    if (nuc == "CTC" || nuc == "CTT")
        return 2.0/3.0;
    // Ile
    if (nuc == "ATA")
        return 1.0;
    if (nuc == "ATC" || nuc == "ATT")
        return 2.0/3.0;
    // Met
    if (nuc == "ATG")
        return 1.0;
    // Val
    if (nuc == "GTA" || nuc == "GTC" || nuc == "GTG" || nuc == "GTT")
        return 2.0/3.0;
    // Ser
    if (nuc == "TCA" || nuc == "TCC" || nuc == "TCG" || nuc == "TCT")
        return 2.0/3.0;
    if (nuc == "AGC" || nuc == "AGT")
        return 2.0/3.0;
    // Pro
    if (nuc == "CCA" || nuc == "CCC" || nuc == "CCG" || nuc == "CCT")
        return 2.0/3.0;
    if (nuc == "ACA" || nuc == "ACC" || nuc == "ACG" || nuc == "ACT")
        return 2.0/3.0;
    if (nuc == "GCA" || nuc == "GCC" || nuc == "GCG" || nuc == "GCT")
        return 2.0/3.0;
    if (nuc == "TAC" || nuc == "TAT")
        return 2.0/3.0;
    if (nuc == "TAA" || nuc == "TAG" || nuc == "TGA")
        return 1.0/3.0;
    if (nuc == "CAC" || nuc == "CAT")
        return 2.0/3.0;
    if (nuc == "CAA" || nuc == "CAG")
        return 2.0/3.0;
    if (nuc == "AAC" || nuc == "AAT")
        return 2.0/3.0;
    if (nuc == "AAA" || nuc == "AAG")
        return 2.0/3.0;
    if (nuc == "GAC" || nuc == "GAT")
        return 2.0/3.0;
    if (nuc == "GAA" || nuc == "GAG")
        return 2.0/3.0;
    if (nuc == "TGC" || nuc == "TGT")
        return 2.0/3.0;
    if (nuc == "TGG")
        return 1.0;
    // Arg
    if (nuc == "CGA" || nuc == "CGC" || nuc == "CGG" || nuc == "CGT" || nuc == "AGA" || nuc == "AGG")
        return 2.0/3.0;
    // Gly
    if (nuc == "GGA" || nuc == "GGC" || nuc == "GGG" || nuc == "GGT")
        return 2.0/3.0;
    else
        return 0;
    
}





inline int getCodonDistance(const std::string& refCdn, const std::string& altCdn) {
    int numDiffs = 0;
    assert(refCdn.length() == altCdn.length());
    for (std::string::size_type i = 0; i != refCdn.length(); i++) {
        if (refCdn[i] != altCdn[i]) {
            numDiffs = numDiffs + 1;
        }
    }
    return numDiffs;
}

inline bool isSingleChangeTransition(const std::string& refCdn, const std::string& altCdn, int diffPos = -1) {
    assert(refCdn.length() == altCdn.length());
    assert(diffPos <= 2);
    if (diffPos == -1) {
        for (std::string::size_type i = 0; i != refCdn.length(); i++) {
            if (refCdn[i] != altCdn[i])
                diffPos = (int)i;
        }
    }

    if (refCdn[diffPos] == 'C' && altCdn[diffPos] == 'T')
        return true;
    else if (refCdn[diffPos] == 'T' && altCdn[diffPos] == 'C')
        return true;
    else if (refCdn[diffPos] == 'A' && altCdn[diffPos] == 'G')
        return true;
    else if (refCdn[diffPos] == 'G' && altCdn[diffPos] == 'A')
        return true;
    else
        return false;
}

static std::unordered_map<std::string, bool> s_mapCodonPairToSynonymous =
{
    // Phe
    { "TTTTTC", true }, { "TTCTTT", true },
    // Leu
    { "TTATTG", true }, { "TTACTA", true }, { "TTGTTA", true }, { "TTGCTG", true },
    { "CTATTA", true }, { "CTACTC", true }, { "CTACTG", true }, { "CTACTT", true },
    { "CTCCTA", true }, { "CTCCTG", true }, { "CTCCTT", true },
    { "CTGTTG", true }, { "CTGCTC", true }, { "CTGCTA", true }, { "CTGCTT", true },
    { "CTTCTA", true }, { "CTTCTC", true }, { "CTTCTG", true },
    // Ile
    { "ATAATC", true }, { "ATAATT", true },
    { "ATCATA", true }, { "ATCATT", true },
    { "ATTATA", true }, { "ATTATC", true },
    // Val
    { "GTAGTC", true }, { "GTAGTG", true }, { "GTAGTT", true },
    { "GTCGTA", true }, { "GTCGTG", true }, { "GTCGTT", true },
    { "GTGGTA", true }, { "GTGGTC", true }, { "GTGGTT", true },
    { "GTTGTA", true }, { "GTTGTC", true }, { "GTTGTG", true },
    // Ser
    { "TCATCC", true }, { "TCATCG", true }, { "TCATCT", true },
    { "TCCTCA", true }, { "TCCTCG", true }, { "TCCTCT", true },
    { "TCGTCA", true }, { "TCGTCC", true }, { "TCGTCT", true },
    { "TCTTCA", true }, { "TCTTCC", true }, { "TCTTCG", true },
    { "AGCAGT", true }, { "AGTAGC", true },
    // Pro
    { "CCACCC", true }, { "CCACCG", true }, { "CCACCT", true },
    { "CCCCCA", true }, { "CCCCCG", true }, { "CCCCCT", true },
    { "CCGCCA", true }, { "CCGCCC", true }, { "CCGCCT", true },
    { "CCTCCA", true }, { "CCTCCC", true }, { "CCTCCG", true },
    // Thr
    { "ACAACC", true }, { "ACAACG", true }, { "ACAACT", true },
    { "ACCACA", true }, { "ACCACG", true }, { "ACCACT", true },
    { "ACGACA", true }, { "ACGACC", true }, { "ACGACT", true },
    { "ACTACA", true }, { "ACTACC", true }, { "ACTACG", true },
    // Ala
    { "GCAGCC", true }, { "GCAGCG", true }, { "GCAGCT", true },
    { "GCCGCA", true }, { "GCCGCG", true }, { "GCCGCT", true },
    { "GCGGCA", true }, { "GCGGCC", true }, { "GCGGCT", true },
    { "GCTGCA", true }, { "GCTGCC", true }, { "GCTGCG", true },
    // Tyr
    { "TACTAT", true }, { "TATTAC", true },
    // Stop
    { "TAATAG", true }, { "TAATGA", true },
    { "TAGTAA", true }, { "TGATAA", true },
    // His
    { "CACCAT", true }, { "CATCAC", true },
    //
    { "CAACAG", true }, { "CAGCAA", true },
    //
    { "AACAAT", true }, { "AATAAC", true },
    //
    { "AAAAAG", true }, { "AAGAAA", true },
    //
    { "GACGAT", true }, { "GATGAC", true },
    //
    { "GAAGAG", true }, { "GAGGAA", true },
    // Cys
    { "TGCTGT", true }, { "TGTTGC", true },
    // Arg
    { "CGACGC", true }, { "CGACGG", true }, { "CGACGT", true }, { "CGAAGA", true },
    { "CGCCGA", true }, { "CGCCGG", true }, { "CGCCGT", true },
    { "CGGCGC", true }, { "CGGCGA", true }, { "CGGCGT", true }, { "CGGAGG", true },
    { "CGTCGA", true }, { "CGTCGC", true }, { "CGTCGG", true },
    { "AGACGA", true }, { "AGAAGG", true }, { "AGGAGA", true }, { "AGGCGG", true },
    // Gly
    { "GGAGGC", true }, { "GGAGGG", true }, { "GGAGGT", true },
    { "GGCGGA", true }, { "GGCGGG", true }, { "GGCGGT", true },
    { "GGGGGA", true }, { "GGGGGC", true }, { "GGGGGT", true },
    { "GGTGGA", true }, { "GGTGGC", true }, { "GGTGGG", true },
};

inline double calculateN(const std::string& refCdn, const std::string& altCdn, int diffNum, bool refAncestral) {
    assert(refCdn.length() == altCdn.length());
    assert(diffNum >= 0 && diffNum <= 3);
    double N = 0;
    
    if (diffNum == 0) {
        return s_mapCodonToExpDist[refCdn];
    }
    
    if (diffNum == 1) {
        if (refAncestral) {
            return s_mapCodonToExpDist[refCdn];
        } else { // Also consider the reverse order of mutations if we don't know the ancestral allele
            return (s_mapCodonToExpDist[refCdn] + s_mapCodonToExpDist[altCdn])/2;
        }
    } else {
        if (diffNum == 2) {
            // Find where the diffs are
            std::vector<std::string::size_type> diffPos;
            for (std::string::size_type i = 0; i != refCdn.length(); i++) {
                if (refCdn[i] != altCdn[i]) {
                    diffPos.push_back(i);
                }
            }
            
            // Then calculate N for the possible mutation paths:
            double Nsum = 0;
            // e.g. TAA -> TGA -> TGG
            std::string stepCdn = refCdn;
            stepCdn[diffPos[0]] = altCdn[diffPos[0]];
            Nsum = (s_mapCodonToExpDist[refCdn] + s_mapCodonToExpDist[stepCdn])/2;
            
            // e.g. TAA -> TAG -> TGG
            stepCdn = refCdn;
            stepCdn[diffPos[1]] = altCdn[diffPos[1]];
            Nsum = Nsum + ((s_mapCodonToExpDist[refCdn] + s_mapCodonToExpDist[stepCdn])/2);
            
            if (refAncestral) {
                return Nsum / 2;
            } else { // Also consider the reverse order of mutations if we don't know the ancestral allele
                // e.g. TAA <- TGA <- TGG
                stepCdn = altCdn;
                stepCdn[diffPos[0]] = refCdn[diffPos[0]];
                Nsum = Nsum + ((s_mapCodonToExpDist[altCdn] + s_mapCodonToExpDist[stepCdn])/2);
                
                // e.g. TAA <- TAG <- TGG
                stepCdn = altCdn; // the reverse order of mutations
                stepCdn[diffPos[1]] = refCdn[diffPos[1]];
                Nsum = Nsum + ((s_mapCodonToExpDist[altCdn] + s_mapCodonToExpDist[stepCdn])/2);
                // Get the average N for the four mutation paths:
                return Nsum / 4;
            }
        }
        
        // This could surely be simplified but I write out all the possible mutation paths explicitly
        // one could even have a lookup table for all the three letter pairs and what N scores they give
        if (diffNum == 3) {
            double Nsum = 0;
            
            // e.g. AAA -> TAA -> TGA -> TGG
            std::string stepCdn = refCdn; stepCdn[0] = altCdn[0];
            std::string step2Cdn = stepCdn; step2Cdn[1] = altCdn[1];
            Nsum = (s_mapCodonToExpDist[refCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3;
            
            // e.g. AAA -> TAA -> TAG -> TGG
            stepCdn = refCdn; stepCdn[0] = altCdn[0];
            step2Cdn = stepCdn; step2Cdn[2] = altCdn[2];
            Nsum = Nsum + ((s_mapCodonToExpDist[refCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
            
            // e.g. AAA -> AGA -> TGA -> TGG
            stepCdn = refCdn; stepCdn[1] = altCdn[1];
            step2Cdn = stepCdn; step2Cdn[0] = altCdn[0];
            Nsum = Nsum + ((s_mapCodonToExpDist[refCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
            
            // e.g. AAA -> AGA -> AGG -> TGG
            stepCdn = refCdn; stepCdn[1] = altCdn[1];
            step2Cdn = stepCdn; step2Cdn[2] = altCdn[2];
            Nsum = Nsum + ((s_mapCodonToExpDist[refCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
            
            // e.g. AAA -> AAG -> AGG -> TGG
            stepCdn = refCdn; stepCdn[2] = altCdn[2];
            step2Cdn = stepCdn; step2Cdn[1] = altCdn[1];
            Nsum = Nsum + ((s_mapCodonToExpDist[refCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
            
            // e.g. AAA -> AAG -> TAG -> TGG
            stepCdn = refCdn; stepCdn[2] = altCdn[2];
            step2Cdn = stepCdn; step2Cdn[0] = altCdn[0];
            Nsum = Nsum + ((s_mapCodonToExpDist[refCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
            
            if (refAncestral) {
                return Nsum / 6;
            } else { // Also consider the reverse order of mutations if we don't know the ancestral allele
                // e.g. AAA <- TAA <- TGA <- TGG
                stepCdn = altCdn; stepCdn[2] = refCdn[2];
                step2Cdn = stepCdn; step2Cdn[1] = refCdn[1];
                Nsum = Nsum + ((s_mapCodonToExpDist[altCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
                
                // e.g. AAA <- TAA <- TAG <- TGG
                stepCdn = altCdn; stepCdn[1] = refCdn[1];
                step2Cdn = stepCdn; step2Cdn[2] = refCdn[2];
                Nsum = Nsum + ((s_mapCodonToExpDist[altCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
                
                // e.g. AAA -> AGA -> TGA -> TGG
                stepCdn = altCdn; stepCdn[2] = refCdn[2];
                step2Cdn = stepCdn; step2Cdn[0] = refCdn[0];
                Nsum = Nsum + ((s_mapCodonToExpDist[altCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
                
                // e.g. AAA <- AGA <- AGG <- TGG
                stepCdn = altCdn; stepCdn[0] = refCdn[0];
                step2Cdn = stepCdn; step2Cdn[2] = refCdn[2];
                Nsum = Nsum + ((s_mapCodonToExpDist[altCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
                
                // e.g. AAA <- AAG <- AGG <- TGG
                stepCdn = altCdn; stepCdn[0] = refCdn[0];
                step2Cdn = stepCdn; step2Cdn[1] = refCdn[1];
                Nsum = Nsum + ((s_mapCodonToExpDist[altCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
                
                // e.g. AAA <- AAG <- TAG <- TGG
                stepCdn = altCdn; stepCdn[1] = refCdn[1];
                step2Cdn = stepCdn; step2Cdn[0] = refCdn[0];
                Nsum = Nsum + ((s_mapCodonToExpDist[altCdn] + s_mapCodonToExpDist[stepCdn] + s_mapCodonToExpDist[step2Cdn])/3);
                
                // Get the average N for the twelve mutation paths:
                return Nsum / 12;
            }
        }
        
    }
    return N;
}


inline double calculateNtS(const std::string& refCdn, const std::string& altCdn, int diffNum, bool refAncestral) {
    assert(refCdn.length() == altCdn.length());
    assert(diffNum >= 0 && diffNum <= 3);
    double NtS = 0;
    
    if (diffNum == 0) {
        return s_mapCodonToExpTs[refCdn];
    }
    
    if (diffNum == 1) {
        if (refAncestral) {
            return s_mapCodonToExpTs[refCdn];
        } else { // Also consider the reverse order of mutations if we don't know the ancestral allele
            return (s_mapCodonToExpTs[refCdn] + s_mapCodonToExpTs[altCdn])/2;
        }
    } else {
        if (diffNum == 2) {
            // Find where the diffs are
            std::vector<std::string::size_type> diffPos;
            for (std::string::size_type i = 0; i != refCdn.length(); i++) {
                if (refCdn[i] != altCdn[i]) {
                    diffPos.push_back(i);
                }
            }
            
            // Then calculate N for the possible mutation paths:
            double Nsum = 0;
            // e.g. TAA -> TGA -> TGG
            std::string stepCdn = refCdn;
            stepCdn[diffPos[0]] = altCdn[diffPos[0]];
            Nsum = (s_mapCodonToExpTs[refCdn] + s_mapCodonToExpTs[stepCdn])/2;
            
            // e.g. TAA -> TAG -> TGG
            stepCdn = refCdn;
            stepCdn[diffPos[1]] = altCdn[diffPos[1]];
            Nsum = Nsum + ((s_mapCodonToExpTs[refCdn] + s_mapCodonToExpTs[stepCdn])/2);
            
            if (refAncestral) {
                return Nsum / 2;
            } else { // Also consider the reverse order of mutations if we don't know the ancestral allele
                // e.g. TAA <- TGA <- TGG
                stepCdn = altCdn;
                stepCdn[diffPos[0]] = refCdn[diffPos[0]];
                Nsum = Nsum + ((s_mapCodonToExpTs[altCdn] + s_mapCodonToExpTs[stepCdn])/2);
                
                // e.g. TAA <- TAG <- TGG
                stepCdn = altCdn; // the reverse order of mutations
                stepCdn[diffPos[1]] = refCdn[diffPos[1]];
                Nsum = Nsum + ((s_mapCodonToExpTs[altCdn] + s_mapCodonToExpTs[stepCdn])/2);
                // Get the average N for the four mutation paths:
                return Nsum / 4;
            }
        }
        
        // This could surely be simplified but I write out all the possible mutation paths explicitly
        // one could even have a lookup table for all the three letter pairs and what N scores they give
        if (diffNum == 3) {
            double Nsum = 0;
            
            // e.g. AAA -> TAA -> TGA -> TGG
            std::string stepCdn = refCdn; stepCdn[0] = altCdn[0];
            std::string step2Cdn = stepCdn; step2Cdn[1] = altCdn[1];
            Nsum = (s_mapCodonToExpTs[refCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3;
            
            // e.g. AAA -> TAA -> TAG -> TGG
            stepCdn = refCdn; stepCdn[0] = altCdn[0];
            step2Cdn = stepCdn; step2Cdn[2] = altCdn[2];
            Nsum = Nsum + ((s_mapCodonToExpTs[refCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
            
            // e.g. AAA -> AGA -> TGA -> TGG
            stepCdn = refCdn; stepCdn[1] = altCdn[1];
            step2Cdn = stepCdn; step2Cdn[0] = altCdn[0];
            Nsum = Nsum + ((s_mapCodonToExpTs[refCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
            
            // e.g. AAA -> AGA -> AGG -> TGG
            stepCdn = refCdn; stepCdn[1] = altCdn[1];
            step2Cdn = stepCdn; step2Cdn[2] = altCdn[2];
            Nsum = Nsum + ((s_mapCodonToExpTs[refCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
            
            // e.g. AAA -> AAG -> AGG -> TGG
            stepCdn = refCdn; stepCdn[2] = altCdn[2];
            step2Cdn = stepCdn; step2Cdn[1] = altCdn[1];
            Nsum = Nsum + ((s_mapCodonToExpTs[refCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
            
            // e.g. AAA -> AAG -> TAG -> TGG
            stepCdn = refCdn; stepCdn[2] = altCdn[2];
            step2Cdn = stepCdn; step2Cdn[0] = altCdn[0];
            Nsum = Nsum + ((s_mapCodonToExpTs[refCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
            
            if (refAncestral) {
                return Nsum / 6;
            } else { // Also consider the reverse order of mutations if we don't know the ancestral allele
                // e.g. AAA <- TAA <- TGA <- TGG
                stepCdn = altCdn; stepCdn[2] = refCdn[2];
                step2Cdn = stepCdn; step2Cdn[1] = refCdn[1];
                Nsum = Nsum + ((s_mapCodonToExpTs[altCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
                
                // e.g. AAA <- TAA <- TAG <- TGG
                stepCdn = altCdn; stepCdn[1] = refCdn[1];
                step2Cdn = stepCdn; step2Cdn[2] = refCdn[2];
                Nsum = Nsum + ((s_mapCodonToExpTs[altCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
                
                // e.g. AAA -> AGA -> TGA -> TGG
                stepCdn = altCdn; stepCdn[2] = refCdn[2];
                step2Cdn = stepCdn; step2Cdn[0] = refCdn[0];
                Nsum = Nsum + ((s_mapCodonToExpTs[altCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
                
                // e.g. AAA <- AGA <- AGG <- TGG
                stepCdn = altCdn; stepCdn[0] = refCdn[0];
                step2Cdn = stepCdn; step2Cdn[2] = refCdn[2];
                Nsum = Nsum + ((s_mapCodonToExpTs[altCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
                
                // e.g. AAA <- AAG <- AGG <- TGG
                stepCdn = altCdn; stepCdn[0] = refCdn[0];
                step2Cdn = stepCdn; step2Cdn[1] = refCdn[1];
                Nsum = Nsum + ((s_mapCodonToExpTs[altCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
                
                // e.g. AAA <- AAG <- TAG <- TGG
                stepCdn = altCdn; stepCdn[1] = refCdn[1];
                step2Cdn = stepCdn; step2Cdn[0] = refCdn[0];
                Nsum = Nsum + ((s_mapCodonToExpTs[altCdn] + s_mapCodonToExpTs[stepCdn] + s_mapCodonToExpTs[step2Cdn])/3);
                
                // Get the average N for the twelve mutation paths:
                return Nsum / 12;
            }
        }
        
    }
    return NtS;
}


inline double calculateNd(const std::string& refCdn, const std::string& altCdn, int diffNum) {
    assert(refCdn.length() == altCdn.length());
    std::string refStep = refCdn;
    int countNS = 0;
    double Nd = 0;
    
    if (diffNum == 1) {
        if (s_mapCodonPairToSynonymous.count(refCdn+altCdn) == 0)
            Nd = 1.0;
    }
    
    if (diffNum == 2) {
        std::vector<std::string::size_type> diffPos;
        for (std::string::size_type i = 0; i != refCdn.length(); i++) {
            if (refCdn[i] != altCdn[i]) {
                diffPos.push_back(i);
            }
        }
        refStep[diffPos[0]] = altCdn[diffPos[0]];
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0)
            countNS++;
        if (s_mapCodonPairToSynonymous.count(refStep+altCdn) == 0)
            countNS++;
        refStep = refCdn; // the reverse order of mutations
        refStep[diffPos[1]] = altCdn[diffPos[1]];
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0)
            countNS++;
        if (s_mapCodonPairToSynonymous.count(refStep+altCdn) == 0)
            countNS++;
        Nd = countNS/2.0;
    }
    
    if (diffNum == 3) { // six different mutation pathways (orders of mutations)
        refStep[0] = altCdn[0]; // 1
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0)
            countNS++;
        string refStep2 = refStep;
        refStep2[1] = altCdn[1];
        if (s_mapCodonPairToSynonymous.count(refStep+refStep2) == 0)
            countNS++;
        if (s_mapCodonPairToSynonymous.count(refStep2+altCdn) == 0)
            countNS++;
        refStep = refCdn;                           // 2
        refStep[2] = altCdn[2];
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0)
            countNS++;
        refStep2 = refStep;
        refStep2[1] = altCdn[1];
        if (s_mapCodonPairToSynonymous.count(refStep+refStep2) == 0)
            countNS++;
        if (s_mapCodonPairToSynonymous.count(refStep2+altCdn) == 0)
            countNS++;
        refStep = refCdn;                           // 3
        refStep[0] = altCdn[0];
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0)
            countNS++;
        refStep2 = refStep;
        refStep2[2] = altCdn[2];
        if (s_mapCodonPairToSynonymous.count(refStep+refStep2) == 0)
            countNS++;
        if (s_mapCodonPairToSynonymous.count(refStep2+altCdn) == 0)
            countNS++;
        refStep = refCdn;                           // 4
        refStep[2] = altCdn[2];
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0)
            countNS++;
        refStep2 = refStep;
        refStep2[0] = altCdn[0];
        if (s_mapCodonPairToSynonymous.count(refStep+refStep2) == 0)
            countNS++;
        if (s_mapCodonPairToSynonymous.count(refStep2+altCdn) == 0)
            countNS++;
        refStep = refCdn;                           // 5
        refStep[1] = altCdn[1];
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0)
            countNS++;
        refStep2 = refStep;
        refStep2[0] = altCdn[0];
        if (s_mapCodonPairToSynonymous.count(refStep+refStep2) == 0)
            countNS++;
        if (s_mapCodonPairToSynonymous.count(refStep2+altCdn) == 0)
            countNS++;
        refStep = refCdn;                           // 6
        refStep[1] = altCdn[1];
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0)
            countNS++;
        refStep2 = refStep;
        refStep2[2] = altCdn[2];
        if (s_mapCodonPairToSynonymous.count(refStep+refStep2) == 0)
            countNS++;
        if (s_mapCodonPairToSynonymous.count(refStep2+altCdn) == 0)
            countNS++;
        Nd = countNS/6.0;
    }
    assert(Nd <= diffNum);
    return Nd;
}

inline std::vector<double> calculateNd_tS_tV(const std::string& refCdn, const std::string& altCdn, int diffNum) {
    assert(refCdn.length() == altCdn.length());
    std::string refStep = refCdn;
    int countNS = 0;
    int countNS_tS = 0;
    int countNS_tV = 0;
    double Nd = 0;
    double tS_Nd = 0;
    double tV_Nd = 0;
    
    if (diffNum == 1) {
        if (s_mapCodonPairToSynonymous.count(refCdn+altCdn) == 0) {
            Nd = 1.0;
            if (isSingleChangeTransition(refCdn,altCdn)) {
                tS_Nd = 1.0; tV_Nd = 0.0;
            } else {
                tS_Nd = 0.0; tV_Nd = 1.0;
            }
        }
    }
    
    if (diffNum == 2) {
        // std::cerr << "diffNum = " << diffNum << std::endl;
        std::vector<std::string::size_type> diffPos;
        for (std::string::size_type i = 0; i != refCdn.length(); i++) {
            if (refCdn[i] != altCdn[i]) {
                diffPos.push_back(i);
            }
        }
        refStep[diffPos[0]] = altCdn[diffPos[0]];
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0) {
            countNS++;
            if (isSingleChangeTransition(refCdn,refStep,(int)diffPos[0]))
                countNS_tS++;
            else
                countNS_tV++;
        }
        if (s_mapCodonPairToSynonymous.count(refStep+altCdn) == 0) {
            countNS++;
            if (isSingleChangeTransition(refStep,altCdn,(int)diffPos[1]))
                countNS_tS++;
            else
                countNS_tV++;
        }
        refStep = refCdn; // the reverse order of mutations
        refStep[diffPos[1]] = altCdn[diffPos[1]];
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0) {
            countNS++;
            if (isSingleChangeTransition(refCdn,refStep,(int)diffPos[1]))
                countNS_tS++;
            else
                countNS_tV++;
        }
        if (s_mapCodonPairToSynonymous.count(refStep+altCdn) == 0) {
            countNS++;
            if (isSingleChangeTransition(refStep,altCdn,(int)diffPos[0]))
                countNS_tS++;
            else
                countNS_tV++;
        }
        Nd = countNS/2.0;
        tS_Nd = countNS_tS/2.0;
        tV_Nd = countNS_tV/2.0;
    }
    
    if (diffNum == 3) { // six different mutation pathways (orders of mutations)
        refStep[0] = altCdn[0]; // 1
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0) {
            countNS++;
            if (isSingleChangeTransition(refCdn,refStep,0))
                countNS_tS++;
            else
                countNS_tV++;
        }
        string refStep2 = refStep;
        refStep2[1] = altCdn[1];
        if (s_mapCodonPairToSynonymous.count(refStep+refStep2) == 0) {
            countNS++;
            if (isSingleChangeTransition(refCdn,refStep,1))
                countNS_tS++;
            else
                countNS_tV++;
        }
        if (s_mapCodonPairToSynonymous.count(refStep2+altCdn) == 0) {
            countNS++;
            if (isSingleChangeTransition(refCdn,refStep,2))
                countNS_tS++;
            else
                countNS_tV++;
        }
        refStep = refCdn;                           // 2
        refStep[2] = altCdn[2];
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0) {
            countNS++;
            if (isSingleChangeTransition(refCdn,refStep,2))
                countNS_tS++;
            else
                countNS_tV++;
        }
        refStep2 = refStep;
        refStep2[1] = altCdn[1];
        if (s_mapCodonPairToSynonymous.count(refStep+refStep2) == 0) {
            countNS++;
            if (isSingleChangeTransition(refCdn,refStep,1))
                countNS_tS++;
            else
                countNS_tV++;
        }
        if (s_mapCodonPairToSynonymous.count(refStep2+altCdn) == 0) {
            countNS++;
            if (isSingleChangeTransition(refCdn,refStep,0))
                countNS_tS++;
            else
                countNS_tV++;
        }
        refStep = refCdn;                           // 3
        refStep[0] = altCdn[0];
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0) {
            countNS++;
            if (isSingleChangeTransition(refCdn,refStep,0))
                countNS_tS++;
            else
                countNS_tV++;
        }
        refStep2 = refStep;
        refStep2[2] = altCdn[2];
        if (s_mapCodonPairToSynonymous.count(refStep+refStep2) == 0) {
            countNS++;
            if (isSingleChangeTransition(refCdn,refStep,2))
                countNS_tS++;
            else
                countNS_tV++;
        }
        if (s_mapCodonPairToSynonymous.count(refStep2+altCdn) == 0) {
            countNS++;
            if (isSingleChangeTransition(refCdn,refStep,1))
                countNS_tS++;
            else
                countNS_tV++;
        }
        refStep = refCdn;                           // 4
        refStep[2] = altCdn[2];
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0) {
            countNS++;
            if (isSingleChangeTransition(refCdn,refStep,2))
                countNS_tS++;
            else
                countNS_tV++;
        }
        refStep2 = refStep;
        refStep2[0] = altCdn[0];
        if (s_mapCodonPairToSynonymous.count(refStep+refStep2) == 0) {
            countNS++;
            if (isSingleChangeTransition(refCdn,refStep,0))
                countNS_tS++;
            else
                countNS_tV++;
        }
        if (s_mapCodonPairToSynonymous.count(refStep2+altCdn) == 0) {
            countNS++;
            if (isSingleChangeTransition(refCdn,refStep,1))
                countNS_tS++;
            else
                countNS_tV++;
        }
        refStep = refCdn;                           // 5
        refStep[1] = altCdn[1];
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0) {
            countNS++;
            if (isSingleChangeTransition(refCdn,refStep,1))
                countNS_tS++;
            else
                countNS_tV++;
        }
        refStep2 = refStep;
        refStep2[0] = altCdn[0];
        if (s_mapCodonPairToSynonymous.count(refStep+refStep2) == 0) {
            countNS++;
            if (isSingleChangeTransition(refCdn,refStep,0))
                countNS_tS++;
            else
                countNS_tV++;
        }
        if (s_mapCodonPairToSynonymous.count(refStep2+altCdn) == 0) {
            countNS++;
            if (isSingleChangeTransition(refCdn,refStep,2))
                countNS_tS++;
            else
                countNS_tV++;
        }
        refStep = refCdn;                           // 6
        refStep[1] = altCdn[1];
        if (s_mapCodonPairToSynonymous.count(refCdn+refStep) == 0) {
            countNS++;
            if (isSingleChangeTransition(refCdn,refStep,1))
                countNS_tS++;
            else
                countNS_tV++;
        }
        refStep2 = refStep;
        refStep2[2] = altCdn[2];
        if (s_mapCodonPairToSynonymous.count(refStep+refStep2) == 0) {
            countNS++;
            if (isSingleChangeTransition(refCdn,refStep,2))
                countNS_tS++;
            else
                countNS_tV++;
        }
        if (s_mapCodonPairToSynonymous.count(refStep2+altCdn) == 0) {
            countNS++;
            if (isSingleChangeTransition(refCdn,refStep,0))
                countNS_tS++;
            else
                countNS_tV++;
        }
        Nd = countNS/6.0;
        tS_Nd = countNS_tS/6.0;
        tV_Nd = countNS_tV/6.0;
    }
    assert(Nd <= diffNum);
    std::vector<double> Nd_tS_Nd; Nd_tS_Nd.push_back(Nd); Nd_tS_Nd.push_back(tS_Nd); Nd_tS_Nd.push_back(tV_Nd);
    return Nd_tS_Nd;
}




inline void addAllPairwiseN_S_Nd_Sd_DifferentIndividuals(const std::vector<string>& altCodons, std::map<std::vector<string>::size_type, int>& haveStop, CDSComparisonMatrices& p) {

    for (std::vector<std::string>::size_type j = 0; j != altCodons.size() - 1; j++) {
        if (haveStop[j] == 1 || getAminoAcid(altCodons[j]) == "Uncrecognised codon....")
            continue; // only consider this individual if it did not have a premature stop codon and there are no unrecognised letters
        for (std::vector<std::string>::size_type k = j+1; k != altCodons.size(); k++) {
            if (haveStop[k] == 1 || getAminoAcid(altCodons[k]) == "Uncrecognised codon....")
                continue;
            int d = getCodonDistance(altCodons[j],altCodons[k]);
            // std::cerr << "Got codon distance: d = " << d << std::endl;
//            double n_d_ijk = calculateNd(altCodons[j],altCodons[k], d);
//            double s_d_ijk = d - n_d_ijk;
            //std::cerr << "N_d_jk[j][k] = " << p.N_d_jk[j][k] << std::endl;
            //print_matrix(N_d_jk, std::cout);
            double n_d_ijk = calculateNd(altCodons[j],altCodons[k], d);
          //  std::cerr << "Calculated Nd; n_d_ijk[0] = " << n_d_ijk[0] << std::endl;
          //  std::cerr << "Calculated Nd; n_d_ijk[1] = " << n_d_ijk[1] << std::endl;
          //  std::cerr << "Calculated Nd; n_d_ijk[2] = " << n_d_ijk[2] << std::endl;
            double s_d_ijk = d - n_d_ijk;
            p.N_d_jk[j][k] = p.N_d_jk[j][k] + n_d_ijk;
            p.S_d_jk[j][k] = p.S_d_jk[j][k] + s_d_ijk;

//          p.N_d_jk[j][k] = p.N_d_jk[j][k] + n_d_ijk;
//          p.S_d_jk[j][k] = p.S_d_jk[j][k] + s_d_ijk;
            
            //  n_di = n_di + n_d_ijk; s_di = s_di + s_d_ijk;
            double N_ijk = calculateN(altCodons[j],altCodons[k], d, false);
            //std::cerr << "Calculated N; N_ijk = " <<N_ijk << std::endl;
            double S_ijk = (3 - N_ijk);
            p.N_jk[j][k] = p.N_jk[j][k] + N_ijk; p.S_jk[j][k] = p.S_jk[j][k] + S_ijk;
            //std::cerr << "altCodons[j] = " << altCodons[j] << "; altCodons[k] = " << altCodons[k] << std::endl;
            //std::cerr << "j = " << j << "; k = " << k << std::endl;
            //std::cerr << "d = " << d << "; n_d_ijk = " << n_d_ijk << "; N_ijk = " << N_ijk << std::endl;
            //std::cerr << "N_d_jk[j][k] = " << p.N_d_jk[j][k] << "; N_jk[j][k] = " << p.N_jk[j][k] << std::endl;
            double N_tS = calculateNtS(altCodons[j],altCodons[k], d, false);
            p.tS_N_jk[j][k] = p.tS_N_jk[j][k] + N_tS;
            p.tS_S_jk[j][k] = p.tS_S_jk[j][k] + (1 - N_tS);
            p.tV_N_jk[j][k] = p.tV_N_jk[j][k] + (N_ijk - N_tS);
            p.tV_S_jk[j][k] = p.tV_S_jk[j][k] + (2 - (N_ijk - N_tS));
            

        }
    }
}

inline void addN_S_Nd_Sd_DifferentIndividualsH1againstH2(const std::vector<string>& altCodons, const std::vector<string>& altCodonsH2, std::map<std::vector<string>::size_type, int>& haveStop, std::map<std::vector<string>::size_type, int>& haveStopH2, CDSComparisonMatrices& p) {
    int numSamples = (int)altCodons.size();
    for (std::vector<std::string>::size_type j = 0; j != numSamples; j++) {
        if (haveStop[j] == 1 || getAminoAcid(altCodons[j]) == "Uncrecognised codon....")
            continue; // only consider this individual if it did not have a premature stop codon
        for (std::vector<std::string>::size_type k = 0; k != numSamples; k++) {
            if (haveStopH2[k] == 1 || getAminoAcid(altCodonsH2[k]) == "Uncrecognised codon....")
                continue;
            if (j != k) {
                int d = getCodonDistance(altCodons[j],altCodonsH2[k]);
                // Nd
                double n_d_ijk = calculateNd(altCodons[j],altCodonsH2[k], d);
                double s_d_ijk = d - n_d_ijk;
                p.N_d_jk[j][k] = p.N_d_jk[j][k] + n_d_ijk;
                p.S_d_jk[j][k] = p.S_d_jk[j][k] + s_d_ijk;
                // N
                double N_ijk = calculateN(altCodons[j],altCodonsH2[k], d, false);
                double S_ijk = (3 - N_ijk);
                p.N_jk[j][k] = p.N_jk[j][k] + N_ijk; p.S_jk[j][k] = p.S_jk[j][k] + S_ijk;
                // NtS
                double N_tS = calculateNtS(altCodons[j],altCodonsH2[k], d, false);
                p.tS_N_jk[j][k] = p.tS_N_jk[j][k] + N_tS;
                p.tS_S_jk[j][k] = p.tS_S_jk[j][k] + (1 - N_tS);
                p.tV_N_jk[j][k] = p.tV_N_jk[j][k] + (N_ijk - N_tS);
                p.tV_S_jk[j][k] = p.tV_S_jk[j][k] + (2 - (N_ijk - N_tS));
            }
        }
    }
}

inline void addN_S_Nd_Sd_SameIndividualsH1againstH2(const std::vector<string>& altCodons, const std::vector<string>& altCodonsH2, std::map<std::vector<string>::size_type, int>& haveStop, std::map<std::vector<string>::size_type, int>& haveStopH2, CDSComparisonMatrices& p) {
    int numSamples = (int)altCodons.size();
    for (std::vector<std::string>::size_type j = 0; j != numSamples; j++) {
        if (haveStop[j] == 1 || getAminoAcid(altCodons[j]) == "Uncrecognised codon....")
            continue; // only consider this individual if it did not have a premature stop codon
        if (haveStopH2[j] == 1 || getAminoAcid(altCodonsH2[j]) == "Uncrecognised codon....")
            continue;
        
        int d = getCodonDistance(altCodons[j],altCodonsH2[j]);
        // Nd
        double n_d_ijk = calculateNd(altCodons[j],altCodonsH2[j], d);
        double s_d_ijk = d - n_d_ijk;
        p.N_d_jk[j][j] = p.N_d_jk[j][j] + n_d_ijk;
        p.S_d_jk[j][j] = p.S_d_jk[j][j] + s_d_ijk;
        // N
        double N_ijk = calculateN(altCodons[j],altCodonsH2[j], d, false);
        double S_ijk = (3 - N_ijk);
        p.N_jk[j][j] = p.N_jk[j][j] + N_ijk; p.S_jk[j][j] = p.S_jk[j][j] + S_ijk;
        // NtS -- for incorporating tS/tV mutation probabilities
        double N_tS = calculateNtS(altCodons[j],altCodonsH2[j], d, false);
        p.tS_N_jk[j][j] = p.tS_N_jk[j][j] + N_tS;
        p.tS_S_jk[j][j] = p.tS_S_jk[j][j] + (1 - N_tS);
        p.tV_N_jk[j][j] = p.tV_N_jk[j][j] + (N_ijk - N_tS);
        p.tV_S_jk[j][j] = p.tV_S_jk[j][j] + (2 - (N_ijk - N_tS));
    }
}




#endif /* defined(__vcf_process__process_vcf_coding_sequences__) */
