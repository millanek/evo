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
void getStatsPhasedSeq(const std::vector<std::string>& allSeqs, const std::string& refSeq, std::vector<std::string>& statsThisGene, std::ofstream*& prematureStopCodonFile);
void getStatsBothPhasedHaps(const std::vector<std::string>& allSeqs, const std::vector<std::string>& allSeqsH2, std::vector<string>& statsThisGene);
std::vector<double> getPhasedPnPs(const std::vector<std::string>& allSeqs);
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

inline bool isSingleChangeTransition(const std::string& refCdn, const std::string& altCdn) {
    assert(refCdn.length() == altCdn.length());
    std::vector<std::string::size_type> diffPos;
    for (std::string::size_type i = 0; i != refCdn.length(); i++) {
        if (refCdn[i] != altCdn[i])
            diffPos.push_back(i);
    }
    assert(diffPos.size() == 1);

    if (refCdn[diffPos[0]] == 'C' && altCdn[diffPos[0]] == 'T')
        return true;
    else if (refCdn[diffPos[0]] == 'T' && altCdn[diffPos[0]] == 'C')
        return true;
    else if (refCdn[diffPos[0]] == 'A' && altCdn[diffPos[0]] == 'G')
        return true;
    else if (refCdn[diffPos[0]] == 'G' && altCdn[diffPos[0]] == 'A')
        return true;
    else
        return false;
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

inline double calculateN(const std::string& refCdn, const std::string& altCdn, int diffNum, bool refAncestral) {
    assert(refCdn.length() == altCdn.length());
    assert(diffNum >= 0 && diffNum <= 3);
    double N = 0;
    
    if (diffNum == 0) {
        return getExpectedNumberOfNonsynonymousSites(refCdn);
    }
    
    if (diffNum == 1) {
        if (refAncestral) {
            return getExpectedNumberOfNonsynonymousSites(refCdn);
        } else { // Also consider the reverse order of mutations if we don't know the ancestral allele
            return (getExpectedNumberOfNonsynonymousSites(refCdn) + getExpectedNumberOfNonsynonymousSites(altCdn))/2;
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
            Nsum = (getExpectedNumberOfNonsynonymousSites(refCdn) + getExpectedNumberOfNonsynonymousSites(stepCdn))/2;
            
            // e.g. TAA -> TAG -> TGG
            stepCdn = refCdn;
            stepCdn[diffPos[1]] = altCdn[diffPos[1]];
            Nsum = Nsum + ((getExpectedNumberOfNonsynonymousSites(refCdn) + getExpectedNumberOfNonsynonymousSites(stepCdn))/2);
            
            if (refAncestral) {
                return Nsum / 2;
            } else { // Also consider the reverse order of mutations if we don't know the ancestral allele
                // e.g. TAA <- TGA <- TGG
                stepCdn = altCdn;
                stepCdn[diffPos[0]] = refCdn[diffPos[0]];
                Nsum = Nsum + ((getExpectedNumberOfNonsynonymousSites(altCdn) + getExpectedNumberOfNonsynonymousSites(stepCdn))/2);
                
                // e.g. TAA <- TAG <- TGG
                stepCdn = altCdn; // the reverse order of mutations
                stepCdn[diffPos[1]] = refCdn[diffPos[1]];
                Nsum = Nsum + ((getExpectedNumberOfNonsynonymousSites(altCdn) + getExpectedNumberOfNonsynonymousSites(stepCdn))/2);
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
            Nsum = (getExpectedNumberOfNonsynonymousSites(refCdn) + getExpectedNumberOfNonsynonymousSites(stepCdn) + getExpectedNumberOfNonsynonymousSites(step2Cdn))/3;
            
            // e.g. AAA -> TAA -> TAG -> TGG
            stepCdn = refCdn; stepCdn[0] = altCdn[0];
            step2Cdn = stepCdn; step2Cdn[2] = altCdn[2];
            Nsum = Nsum + ((getExpectedNumberOfNonsynonymousSites(refCdn) + getExpectedNumberOfNonsynonymousSites(stepCdn) + getExpectedNumberOfNonsynonymousSites(step2Cdn))/3);
            
            // e.g. AAA -> AGA -> TGA -> TGG
            stepCdn = refCdn; stepCdn[1] = altCdn[1];
            step2Cdn = stepCdn; step2Cdn[0] = altCdn[0];
            Nsum = Nsum + ((getExpectedNumberOfNonsynonymousSites(refCdn) + getExpectedNumberOfNonsynonymousSites(stepCdn) + getExpectedNumberOfNonsynonymousSites(step2Cdn))/3);
            
            // e.g. AAA -> AGA -> AGG -> TGG
            stepCdn = refCdn; stepCdn[1] = altCdn[1];
            step2Cdn = stepCdn; step2Cdn[2] = altCdn[2];
            Nsum = Nsum + ((getExpectedNumberOfNonsynonymousSites(refCdn) + getExpectedNumberOfNonsynonymousSites(stepCdn) + getExpectedNumberOfNonsynonymousSites(step2Cdn))/3);
            
            // e.g. AAA -> AAG -> AGG -> TGG
            stepCdn = refCdn; stepCdn[2] = altCdn[2];
            step2Cdn = stepCdn; step2Cdn[1] = altCdn[1];
            Nsum = Nsum + ((getExpectedNumberOfNonsynonymousSites(refCdn) + getExpectedNumberOfNonsynonymousSites(stepCdn) + getExpectedNumberOfNonsynonymousSites(step2Cdn))/3);
            
            // e.g. AAA -> AAG -> TAG -> TGG
            stepCdn = refCdn; stepCdn[2] = altCdn[2];
            step2Cdn = stepCdn; step2Cdn[0] = altCdn[0];
            Nsum = Nsum + ((getExpectedNumberOfNonsynonymousSites(refCdn) + getExpectedNumberOfNonsynonymousSites(stepCdn) + getExpectedNumberOfNonsynonymousSites(step2Cdn))/3);
            
            if (refAncestral) {
                return Nsum / 6;
            } else { // Also consider the reverse order of mutations if we don't know the ancestral allele
                // e.g. AAA <- TAA <- TGA <- TGG
                stepCdn = altCdn; stepCdn[2] = refCdn[2];
                step2Cdn = stepCdn; step2Cdn[1] = refCdn[1];
                Nsum = Nsum + ((getExpectedNumberOfNonsynonymousSites(altCdn) + getExpectedNumberOfNonsynonymousSites(stepCdn) + getExpectedNumberOfNonsynonymousSites(step2Cdn))/3);
                
                // e.g. AAA <- TAA <- TAG <- TGG
                stepCdn = altCdn; stepCdn[1] = refCdn[1];
                step2Cdn = stepCdn; step2Cdn[2] = refCdn[2];
                Nsum = Nsum + ((getExpectedNumberOfNonsynonymousSites(altCdn) + getExpectedNumberOfNonsynonymousSites(stepCdn) + getExpectedNumberOfNonsynonymousSites(step2Cdn))/3);
                
                // e.g. AAA -> AGA -> TGA -> TGG
                stepCdn = altCdn; stepCdn[2] = refCdn[2];
                step2Cdn = stepCdn; step2Cdn[0] = refCdn[0];
                Nsum = Nsum + ((getExpectedNumberOfNonsynonymousSites(altCdn) + getExpectedNumberOfNonsynonymousSites(stepCdn) + getExpectedNumberOfNonsynonymousSites(step2Cdn))/3);
                
                // e.g. AAA <- AGA <- AGG <- TGG
                stepCdn = altCdn; stepCdn[0] = refCdn[0];
                step2Cdn = stepCdn; step2Cdn[2] = refCdn[2];
                Nsum = Nsum + ((getExpectedNumberOfNonsynonymousSites(altCdn) + getExpectedNumberOfNonsynonymousSites(stepCdn) + getExpectedNumberOfNonsynonymousSites(step2Cdn))/3);
                
                // e.g. AAA <- AAG <- AGG <- TGG
                stepCdn = altCdn; stepCdn[0] = refCdn[0];
                step2Cdn = stepCdn; step2Cdn[1] = refCdn[1];
                Nsum = Nsum + ((getExpectedNumberOfNonsynonymousSites(altCdn) + getExpectedNumberOfNonsynonymousSites(stepCdn) + getExpectedNumberOfNonsynonymousSites(step2Cdn))/3);
                
                // e.g. AAA <- AAG <- TAG <- TGG
                stepCdn = altCdn; stepCdn[1] = refCdn[1];
                step2Cdn = stepCdn; step2Cdn[0] = refCdn[0];
                Nsum = Nsum + ((getExpectedNumberOfNonsynonymousSites(altCdn) + getExpectedNumberOfNonsynonymousSites(stepCdn) + getExpectedNumberOfNonsynonymousSites(step2Cdn))/3);
                
                // Get the average N for the twelve mutation paths:
                return Nsum / 12;
            }
        }
        
    }
    return N;
}


inline double calculateNd(const std::string& refCdn, const std::string& altCdn, int diffNum) {
    assert(refCdn.length() == altCdn.length());
    std::string refStep = refCdn;
    int countNS = 0;
    double Nd = 0;
    
    if (diffNum == 1) {
        if (!isSingleChangeSynonymous(refCdn,altCdn)) {
            Nd = 1.0;
        }
    }
    
    if (diffNum == 2) {
        std::vector<std::string::size_type> diffPos;
        for (std::string::size_type i = 0; i != refCdn.length(); i++) {
            if (refCdn[i] != altCdn[i]) {
                diffPos.push_back(i);
            }
        }
        refStep[diffPos[0]] = altCdn[diffPos[0]];
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        refStep[diffPos[1]] = altCdn[diffPos[1]];
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        refStep = refCdn; // the reverse order of mutations
        refStep[diffPos[1]] = altCdn[diffPos[1]];
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        refStep[diffPos[0]] = altCdn[diffPos[0]];
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        Nd = countNS/2.0;
    }
    
    if (diffNum == 3) { // six different mutation pathways (orders of mutations)
        refStep[0] = altCdn[0]; // 1
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        refStep[1] = altCdn[1];
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        refStep[2] = altCdn[2];
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        refStep = refCdn;                           // 2
        refStep[2] = altCdn[2];
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        refStep[1] = altCdn[1];
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        refStep[0] = altCdn[0];
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        refStep = refCdn;                           // 3
        refStep[0] = altCdn[0];
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        refStep[2] = altCdn[2];
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        refStep[1] = altCdn[1];
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        refStep = refCdn;                           // 4
        refStep[2] = altCdn[2];
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        refStep[0] = altCdn[0];
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        refStep[1] = altCdn[1];
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        refStep = refCdn;                           // 5
        refStep[1] = altCdn[1];
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        refStep[0] = altCdn[0];
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        refStep[2] = altCdn[2];
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        refStep = refCdn;                           // 6
        refStep[1] = altCdn[1];
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        refStep[2] = altCdn[2];
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        refStep[0] = altCdn[0];
        if (!isSingleChangeSynonymous(refCdn, refStep))
            countNS++;
        Nd = countNS/6.0;
    }
    assert(Nd <= diffNum);
    return Nd;
}

inline void addAllPairwiseN_S_Nd_Sd_DifferentIndividuals(const std::vector<string>& altCodons, std::map<std::vector<string>::size_type, int>& haveStop, std::vector<std::vector<double> >& N_d_jk, std::vector<std::vector<double> >& N_jk, std::vector<std::vector<double> >& S_d_jk, std::vector<std::vector<double> >& S_jk) {

    for (std::vector<std::string>::size_type j = 0; j != altCodons.size() - 1; j++) {
        if (haveStop[j] == 1)
            continue; // only consider this individual if it did not have a premature stop codon
        for (std::vector<std::string>::size_type k = j+1; k != altCodons.size(); k++) {
            if (haveStop[k] == 1)
                continue;
            int d = getCodonDistance(altCodons[j],altCodons[k]);
            //std::cerr << "Got codon distance: d = " << d << std::endl;
            double n_d_ijk = calculateNd(altCodons[j],altCodons[k], d);
            //std::cerr << "Calculated Nd; n_d_ijk = " << n_d_ijk << std::endl;
            double s_d_ijk = d - n_d_ijk;
            //std::cerr << "N_d_jk[j][k] = " << N_d_jk[j][k] << std::endl;
            //print_matrix(N_d_jk, std::cout);
            N_d_jk[j][k] = N_d_jk[j][k] + n_d_ijk;
            S_d_jk[j][k] = S_d_jk[j][k] + s_d_ijk;
            //  n_di = n_di + n_d_ijk; s_di = s_di + s_d_ijk;
            double N_ijk = calculateN(altCodons[j],altCodons[k], d, false);
            //std::cerr << "Calculated N; N_ijk = " << N_ijk << std::endl;
            double S_ijk = (3 - N_ijk);
            N_jk[j][k] = N_jk[j][k] + N_ijk; S_jk[j][k] = S_jk[j][k] + S_ijk;
            //std::cerr << "altCodons[j] = " << altCodons[j] << "; altCodons[k] = " << altCodons[k] << std::endl;
            //std::cerr << "j = " << j << "; k = " << k << std::endl;
            //std::cerr << "d = " << d << "; n_d_ijk = " << n_d_ijk << "; N_ijk = " << N_ijk << std::endl;
            //std::cerr << "N_d_jk[j][k] = " << N_d_jk[j][k] << "; N_jk[j][k] = " << N_jk[j][k] << std::endl;

        }
    }
}

inline void addN_S_Nd_Sd_DifferentIndividualsH1againstH2(const std::vector<string>& altCodons, const std::vector<string>& altCodonsH2, std::map<std::vector<string>::size_type, int>& haveStop, std::map<std::vector<string>::size_type, int>& haveStopH2, std::vector<std::vector<double> >& N_d_jk, std::vector<std::vector<double> >& N_jk, std::vector<std::vector<double> >& S_d_jk, std::vector<std::vector<double> >& S_jk) {
    int numSamples = (int)altCodons.size();
    for (std::vector<std::string>::size_type j = 0; j != numSamples; j++) {
        if (haveStop[j] == 1)
            continue; // only consider this individual if it did not have a premature stop codon
        for (std::vector<std::string>::size_type k = 0; k != numSamples; k++) {
            if (haveStopH2[k] == 1)
                continue;
            if (j != k) {
                int d = getCodonDistance(altCodons[j],altCodonsH2[k]);
                double n_d_ijk = calculateNd(altCodons[j],altCodonsH2[k], d);
                double s_d_ijk = d - n_d_ijk;
                N_d_jk[j][k] = N_d_jk[j][k] + n_d_ijk;
                S_d_jk[j][k] = S_d_jk[j][k] + s_d_ijk;
                double N_ijk = calculateN(altCodons[j],altCodonsH2[k], d, false);
                double S_ijk = (3 - N_ijk);
                N_jk[j][k] = N_jk[j][k] + N_ijk; S_jk[j][k] = S_jk[j][k] + S_ijk;
            }
        }
    }
}


#endif /* defined(__vcf_process__process_vcf_coding_sequences__) */
