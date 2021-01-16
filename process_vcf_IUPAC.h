//
//  process_vcf_IUPAC.h
//  vcf_process
//
//  Created by Milan Malinsky on 16/09/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#ifndef __vcf_process__process_vcf_IUPAC__
#define __vcf_process__process_vcf_IUPAC__

#include <iostream>
#include <assert.h>

// IUPAC aplphabet tools
// Complement a base using the full IUPAC alphabet
// Also allows for lowercase bases
char complementIUPAC(char c);
// Reverse complement a sequence using the full iupac alphabet
std::string reverseComplementIUPAC(const std::string& seq);
// Return IUPAC nucleotide codes for heterozygous calls
std::string getAmbiguityCode(std::string b1, std::string b2);
// Given reference DNA base and IUPAC encoded het - returns DNA base for the alteranative allele
inline char disambiguateIUPAC(char refBase, char IUPACbase) {
    char altBase;
    switch (IUPACbase)
    {
        case 'K': altBase = (refBase == 'G') ? 'T' : 'G'; break;
        case 'M': altBase = (refBase == 'A') ? 'C' : 'A'; break;
        case 'R': altBase = (refBase == 'A') ? 'G' : 'A'; break;
        case 'S': altBase = (refBase == 'C') ? 'G' : 'C'; break;
        case 'W': altBase = (refBase == 'A') ? 'T' : 'A'; break;
        case 'Y': altBase = (refBase == 'C') ? 'T' : 'C'; break;
        default: assert(false);
    }
    return altBase;
}

// Returns both DNA bases for het ambiguous codes
inline std::string returnHetIUPAC(char IUPACbase) {
    char firstBase; char secondBase;
    switch (std::toupper(IUPACbase))
    {
        case 'K': firstBase = 'T'; secondBase = 'G'; break;
        case 'M': firstBase = 'C'; secondBase = 'A'; break;
        case 'R': firstBase = 'G'; secondBase = 'A'; break;
        case 'S': firstBase = 'G'; secondBase = 'C'; break;
        case 'W': firstBase = 'T'; secondBase = 'A'; break;
        case 'Y': firstBase = 'T'; secondBase = 'C'; break;
        default: assert(false);
    }
    std::string bases = ""; bases = firstBase + secondBase;
    return bases;
}



inline bool isDNAonly(char base) {
    bool DNAonly;
    switch (base)
    {
        case 'A': DNAonly = true; break;
        case 'C': DNAonly = true; break;
        case 'G': DNAonly = true; break;
        case 'T': DNAonly = true; break;
        default: DNAonly = false;
    }
    return DNAonly;
}

inline bool isDNAonlySeq(const std::string& seq) {
    bool DNAonly;
    for (std::string::size_type i = 0; i != seq.length(); i++) {
        DNAonly = isDNAonly(seq[i]);
        if (!DNAonly)
            return false;
    }
    return true;
}

#endif /* defined(__vcf_process__process_vcf_IUPAC__) */
