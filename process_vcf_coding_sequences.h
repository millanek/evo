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

// Process the fifth column from the annatation file to get this gene name
// and information about whether this gene is 'partial' (passed on through the bool parameter)
std::string getGeneName(const std::string& geneColumn, bool& bPartial);

// Various QC checks to verify that we really have a protein coding sequence (length, start/stop codons)
// Return false if the length of the coding sequence is not divisible by three
bool codingSequenceErrorChecks(const std::string& geneSeq, const std::string& transcriptName, const std::vector<std::vector<std::string> >& annotation, const int k, std::ofstream*& badStartStopCodonFile);
// Calculate some statistics about the sequences
void getCodingSequenceStatsPhasedSeq(const std::vector<std::string>& allSeqs, const std::string& refSeq, const std::string& transcriptName, std::vector<std::string>& statsThisGene, std::ofstream*& prematureStopCodonFile, const std::vector<std::string>& sampleNames,Annotation& wgAnnotation);
void getCodingSequenceStatsIUPAC(const std::vector<std::string>& allSeqs, const std::string& refSeq, const std::string& transcriptName, std::vector<std::string>& statsThisGene, std::ofstream*& prematureStopCodonFile, const std::vector<std::string>& sampleNames,Annotation& wgAnnotation);

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


inline double getExpectedNumberOfNonsynonymousSites(const std::string& nuc)
{
    if (nuc == "TTT" || nuc == "TTC")
        return 2.6666666;
    if (nuc == "TTA" || nuc == "TTG" || nuc == "CTA" || nuc == "CTC" || nuc == "CTG" || nuc == "CTT")
        return 1.6666666;
    if (nuc == "ATA" || nuc == "ATC" || nuc == "ATT")
        return 2.3333333;
    if (nuc == "ATG")
        return 3.0;
    if (nuc == "GTA" || nuc == "GTC" || nuc == "GTG" || nuc == "GTT")
        return 2.0;
    if (nuc == "TCA" || nuc == "TCC" || nuc == "TCG" || nuc == "TCT" || nuc == "AGC" || nuc == "AGT")
        return 2.2222222;
    if (nuc == "CCA" || nuc == "CCC" || nuc == "CCG" || nuc == "CCT")
        return 2.0;
    if (nuc == "ACA" || nuc == "ACC" || nuc == "ACG" || nuc == "ACT")
        return 2.0;
    if (nuc == "GCA" || nuc == "GCC" || nuc == "GCG" || nuc == "GCT")
        return 2.0;
    if (nuc == "TAC" || nuc == "TAT")
        return 2.6666666;
    if (nuc == "TAA" || nuc == "TAG" || nuc == "TGA")
        return 2.3333333;
    if (nuc == "CAC" || nuc == "CAT")
        return 2.6666666;
    if (nuc == "CAA" || nuc == "CAG")
        return 2.6666666;
    if (nuc == "AAC" || nuc == "AAT")
        return 2.6666666;
    if (nuc == "AAA" || nuc == "AAG")
        return 2.6666666;
    if (nuc == "GAC" || nuc == "GAT")
        return 2.6666666;
    if (nuc == "GAA" || nuc == "GAG")
        return 2.6666666;
    if (nuc == "TGC" || nuc == "TGT")
        return 2.6666666;
    if (nuc == "TGG")
        return 3.0;
    if (nuc == "CGA" || nuc == "CGC" || nuc == "CGG" || nuc == "CGT" || nuc == "AGA" || nuc == "AGG")
        return 1.6666666;
    if (nuc == "GGA" || nuc == "GGC" || nuc == "GGG" || nuc == "GGT")
        return 2.0;
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

inline bool isSingleChangeSynonymous(const std::string& refCdn, const std::string& altCdn) {
    // Phe
    if((refCdn == "TTT" && altCdn == "TTC") || (refCdn == "TTC" && altCdn == "TTT"))
        return true;
    // Leu
    else if((refCdn == "TTA" && altCdn == "TTG") || (refCdn == "TTA" && altCdn == "CTA"))
        return true;
    else if((refCdn == "TTG" && altCdn == "TTA") || (refCdn == "TTG" && altCdn == "CTG"))
        return true;
    else if((refCdn == "CTA" && altCdn == "TTA") || (refCdn == "CTA" && altCdn == "CTC") || (refCdn == "CTA" && altCdn == "CTG") || (refCdn == "CTA" && altCdn == "CTT"))
        return true;
    else if((refCdn == "CTC" && altCdn == "CTA") || (refCdn == "CTC" && altCdn == "CTG") || (refCdn == "CTC" && altCdn == "CTT"))
        return true;
    else if((refCdn == "CTG" && altCdn == "TTG") || (refCdn == "CTG" && altCdn == "CTA") || (refCdn == "CTG" && altCdn == "CTC") || (refCdn == "CTG" && altCdn == "CTT"))
        return true;
    else if((refCdn == "CTT" && altCdn == "CTA") || (refCdn == "CTT" && altCdn == "CTC") || (refCdn == "CTT" && altCdn == "CTG"))
        return true;
    // Ile
    else if((refCdn == "ATA" && altCdn == "ATC") || (refCdn == "ATA" && altCdn == "ATT"))
        return true;
    else if((refCdn == "ATC" && altCdn == "ATA") || (refCdn == "ATC" && altCdn == "ATT"))
        return true;
    else if((refCdn == "ATT" && altCdn == "ATA") || (refCdn == "ATT" && altCdn == "ATC"))
        return true;
    // Val
    else if((refCdn == "GTA" && altCdn == "GTC") || (refCdn == "GTA" && altCdn == "GTG") || (refCdn == "GTA" && altCdn == "GTT"))
        return true;
    else if((refCdn == "GTC" && altCdn == "GTA") || (refCdn == "GTC" && altCdn == "GTG") || (refCdn == "GTC" && altCdn == "GTT"))
        return true;
    else if((refCdn == "GTG" && altCdn == "GTC") || (refCdn == "GTG" && altCdn == "GTA") || (refCdn == "GTG" && altCdn == "GTT"))
        return true;
    else if((refCdn == "GTT" && altCdn == "GTC") || (refCdn == "GTT" && altCdn == "GTG") || (refCdn == "GTT" && altCdn == "GTA"))
        return true;
    // Pro
    else if((refCdn == "CCA" && altCdn == "CCC") || (refCdn == "CCA" && altCdn == "CCG") || (refCdn == "CCA" && altCdn == "CCT"))
        return true;
    else if((refCdn == "CCC" && altCdn == "CCA") || (refCdn == "CCC" && altCdn == "CCG") || (refCdn == "CCC" && altCdn == "CCT"))
        return true;
    else if((refCdn == "CCG" && altCdn == "CCC") || (refCdn == "CCG" && altCdn == "CCA") || (refCdn == "CCG" && altCdn == "CCT"))
        return true;
    else if((refCdn == "CCT" && altCdn == "CCC") || (refCdn == "CCT" && altCdn == "CCG") || (refCdn == "CCT" && altCdn == "CCA"))
        return true;
    // Thr
    else if((refCdn == "ACA" && altCdn == "ACC") || (refCdn == "ACA" && altCdn == "ACG") || (refCdn == "ACA" && altCdn == "ACT"))
        return true;
    else if((refCdn == "ACC" && altCdn == "ACA") || (refCdn == "ACC" && altCdn == "ACG") || (refCdn == "ACC" && altCdn == "ACT"))
        return true;
    else if((refCdn == "ACG" && altCdn == "ACC") || (refCdn == "ACG" && altCdn == "ACA") || (refCdn == "ACG" && altCdn == "ACT"))
        return true;
    else if((refCdn == "ACT" && altCdn == "ACC") || (refCdn == "ACT" && altCdn == "ACG") || (refCdn == "ACT" && altCdn == "ACA"))
        return true;
    // Ala
    else if((refCdn == "GCA" && altCdn == "GCC") || (refCdn == "GCA" && altCdn == "GCG") || (refCdn == "GCA" && altCdn == "GCT"))
        return true;
    else if((refCdn == "GCC" && altCdn == "GCA") || (refCdn == "GCC" && altCdn == "GCG") || (refCdn == "GCC" && altCdn == "GCT"))
        return true;
    else if((refCdn == "GCG" && altCdn == "GCC") || (refCdn == "GCG" && altCdn == "GCA") || (refCdn == "GCG" && altCdn == "GCT"))
        return true;
    else if((refCdn == "GCT" && altCdn == "GCC") || (refCdn == "GCT" && altCdn == "GCG") || (refCdn == "GCT" && altCdn == "GCA"))
        return true;
    // Tyr
    if((refCdn == "TAC" && altCdn == "TAT") || (refCdn == "TAT" && altCdn == "TAC"))
        return true;
    // His
    if((refCdn == "CAC" && altCdn == "CAT") || (refCdn == "CAT" && altCdn == "CAC"))
        return true;
    // Gln
    if((refCdn == "CAA" && altCdn == "CAG") || (refCdn == "CAG" && altCdn == "CAA"))
        return true;
    // Asn
    if((refCdn == "AAC" && altCdn == "AAT") || (refCdn == "AAT" && altCdn == "AAC"))
        return true;
    // Lys
    if((refCdn == "AAA" && altCdn == "AAG") || (refCdn == "AAG" && altCdn == "AAA"))
        return true;
    // Asp
    if((refCdn == "GAC" && altCdn == "GAT") || (refCdn == "GAT" && altCdn == "GAC"))
        return true;
    // Glu
    if((refCdn == "GAA" && altCdn == "GAG") || (refCdn == "GAG" && altCdn == "GAA"))
        return true;
    // Cys
    if((refCdn == "TGC" && altCdn == "TGT") || (refCdn == "TGT" && altCdn == "TGC"))
        return true;
    // Arg
    else if((refCdn == "CGA" && altCdn == "CGC") || (refCdn == "CGA" && altCdn == "CGG") || (refCdn == "CGA" && altCdn == "CGT") || (refCdn == "CGA" && altCdn == "AGA"))
        return true;
    else if((refCdn == "CGC" && altCdn == "CGA") || (refCdn == "CGC" && altCdn == "CGG") || (refCdn == "CGC" && altCdn == "CGT"))
        return true;
    else if((refCdn == "CGG" && altCdn == "CGC") || (refCdn == "CGG" && altCdn == "CGA") || (refCdn == "CGG" && altCdn == "CGT") || (refCdn == "CGG" && altCdn == "AGG"))
        return true;
    else if((refCdn == "CGT" && altCdn == "CGC") || (refCdn == "CGT" && altCdn == "CGG") || (refCdn == "CGT" && altCdn == "CGA"))
        return true;
    else if((refCdn == "AGA" && altCdn == "AGG") || (refCdn == "AGG" && altCdn == "AGA") || (refCdn == "AGA" && altCdn == "CGA") || (refCdn == "AGG" && altCdn == "CGG"))
        return true;
    // Gly
    else if((refCdn == "GGA" && altCdn == "GGC") || (refCdn == "GGA" && altCdn == "GGG") || (refCdn == "GGA" && altCdn == "GGT"))
        return true;
    else if((refCdn == "GGC" && altCdn == "GGA") || (refCdn == "GGC" && altCdn == "GGG") || (refCdn == "GGC" && altCdn == "GGT"))
        return true;
    else if((refCdn == "GGG" && altCdn == "GGC") || (refCdn == "GGG" && altCdn == "GGA") || (refCdn == "GGG" && altCdn == "GGT"))
        return true;
    else if((refCdn == "GGT" && altCdn == "GGC") || (refCdn == "GGT" && altCdn == "GGG") || (refCdn == "GGT" && altCdn == "GGA"))
        return true;
    // Ser
    else if((refCdn == "TCA" && altCdn == "TCC") || (refCdn == "TCA" && altCdn == "TCG") || (refCdn == "TCA" && altCdn == "TCT"))
        return true;
    else if((refCdn == "TCC" && altCdn == "TCA") || (refCdn == "TCC" && altCdn == "TCG") || (refCdn == "TCC" && altCdn == "TCT"))
        return true;
    else if((refCdn == "TCG" && altCdn == "TCC") || (refCdn == "TCG" && altCdn == "TCA") || (refCdn == "TCG" && altCdn == "TCT"))
        return true;
    else if((refCdn == "TCT" && altCdn == "TCC") || (refCdn == "TCT" && altCdn == "TCG") || (refCdn == "TCT" && altCdn == "TCA"))
        return true;
    else if((refCdn == "AGC" && altCdn == "AGT") || (refCdn == "AGT" && altCdn == "AGC"))
        return true;
    // Stop
    else if((refCdn == "TAA" && altCdn == "TAG") || (refCdn == "TAA" && altCdn == "TGA"))
        return true;
    else if((refCdn == "TAG" && altCdn == "TAA") || (refCdn == "TAG" && altCdn == "TGA"))
        return true;
    else if((refCdn == "TGA" && altCdn == "TAA") || (refCdn == "TGA" && altCdn == "TAG"))
        return true;
    else
        return false;
}



/*inline void findCodonPaths(const std::string& refCdn, const std::string& altCdn) {
    assert(refCdn.lenGCh() == altCdn.length());
    for (std::string::size_type i = 0; i != refCdn.length(); i++) {
        if (refCdn[i] != altCdn[i]) {
            numDiffs = numDiffs + 1;
        }
    }
} */




#endif /* defined(__vcf_process__process_vcf_coding_sequences__) */
