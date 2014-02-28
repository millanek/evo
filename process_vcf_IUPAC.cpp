//
//  process_vcf_IUPAC.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 16/09/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include "process_vcf_IUPAC.h"



// Reverse complement a sequence using the full iupac alphabet
std::string reverseComplementIUPAC(const std::string& seq)
{
    std::string out(seq.length(), 'A');
    int last_pos = (int)seq.length() - 1;
    for(int i = last_pos; i >= 0; --i)
    {
        out[last_pos - i] = complementIUPAC(seq[i]);
    }
    return out;
}


// Complement a base using the full IUPAC alphabet
// Also allows for lowercase bases
char complementIUPAC(char c)
{
    char cmp = '\0';
    bool is_lc = std::islower(c);
    
    switch(std::toupper(c)) {
        case 'A': cmp = 'T'; break;
        case 'C': cmp = 'G'; break;
        case 'G': cmp = 'C'; break;
        case 'T': cmp = 'A'; break;
        case 'M': cmp = 'K'; break;
        case 'R': cmp = 'Y'; break;
        case 'W': cmp = 'W'; break;
        case 'S': cmp = 'S'; break;
        case 'Y': cmp = 'R'; break;
        case 'K': cmp = 'M'; break;
        case 'V': cmp = 'B'; break;
        case 'H': cmp = 'D'; break;
        case 'D': cmp = 'H'; break;
        case 'B': cmp = 'V'; break;
        case 'N': cmp = 'N'; break;
        default:
            assert(false);
    }
    
    if(is_lc)
        cmp = std::tolower(cmp);
    return cmp;
}

// Return IUPAC nucleotide codes for heterozygous calls
std::string getAmbiguityCode(std::string b1, std::string b2) {
    if ((b1 == "A" && b2 == "C") || (b1 == "C" && b2 == "A"))
        return "M";
    else if ((b1 == "A" && b2 == "G") || (b1 == "G" && b2 == "A"))
        return "R";
    else if ((b1 == "A" && b2 == "T") || (b1 == "T" && b2 == "A"))
        return "W";
    else if ((b1 == "C" && b2 == "G") || (b1 == "G" && b2 == "C"))
        return "S";
    else if ((b1 == "C" && b2 == "T") || (b1 == "T" && b2 == "C"))
        return "Y";
    else if ((b1 == "G" && b2 == "T") || (b1 == "T" && b2 == "G"))
        return "K";
    else
        return "?";
}



