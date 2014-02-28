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
void getCodingSequenceStats(const std::vector<std::string>& allSeqs, const std::string& refSeq, const std::string& transcriptName, std::vector<std::string>& statsThisGene, std::ofstream*& prematureStopCodonFile, const std::vector<std::string>& sampleNames,Annotation& wgAnnotation);

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

#endif /* defined(__vcf_process__process_vcf_coding_sequences__) */
