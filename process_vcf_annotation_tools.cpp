//
//  process_vcf_annotation_tools.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 16/09/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include "process_vcf_annotation_tools.h"


std::string getIndividualSequenceForThisRegion(const std::vector<std::string>& thisRegionAnnotation, const std::string& strand, const std::string& currentIndividualWholeScaffoldSequence) {
    string geneSeq = "";
    for (std::vector<std::string>::size_type j = 0; j != thisRegionAnnotation.size(); j++) {
        std::vector<string> annotLineVec = split(thisRegionAnnotation[j], '\t');
        geneSeq = geneSeq + currentIndividualWholeScaffoldSequence.substr(atoi(annotLineVec[1].c_str())-1,atoi(annotLineVec[2].c_str())-atoi(annotLineVec[1].c_str()) + 1);
    }
    if (strand == "-")
        geneSeq = reverseComplementIUPAC(geneSeq);
    
    return geneSeq;
}


std::string getReferenceForThisRegion(const std::vector<std::string>& thisRegionAnnotation, const std::string& strand, const std::string& currentScaffoldReference) {
    // std::cerr << "generating gene refseq" << std::endl;
    std::string refSeq = "";
    for (std::vector<std::string>::size_type j = 0; j != thisRegionAnnotation.size(); j++) {
        std::vector<string> annotLineVec = split(thisRegionAnnotation[j], '\t');
        refSeq = refSeq + currentScaffoldReference.substr(atoi(annotLineVec[1].c_str())-1,atoi(annotLineVec[2].c_str())-atoi(annotLineVec[1].c_str()) + 1);
    }
    if (strand == "-")
        refSeq = reverseComplementIUPAC(refSeq);
    // std::cerr << "Done" << std::endl;
    return refSeq;
}
