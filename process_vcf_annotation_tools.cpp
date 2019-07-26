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


std::vector<string> Annotation::getSNPgeneDetails(const string& SNPscaffold, const int SNPlocus) {
    std::vector<string> SNPgeneDetails;
    std::vector<string> scaffoldTranscriptStartEnd = transcriptStartEndMap[SNPscaffold];
    string inGene = ""; string SNPcategory = "nonCoding";
    for (std::vector<std::vector<string> >::size_type i = 0; i != scaffoldTranscriptStartEnd.size(); i++) {
        std::vector<string> startEndVec = split(scaffoldTranscriptStartEnd[i], '\t');
        string thisTranscript = startEndVec[0];
        int geneStart = atoi(startEndVec[1].c_str()); int geneEnd = atoi(startEndVec[2].c_str()); string strand = startEndVec[3];
        //if (SNPlocus == 20001) { print_vector_stream(startEndVec, std::cerr);}
        if (strand == "+") {
            if (SNPlocus >= (geneStart-3000) && SNPlocus < geneStart) {
                SNPcategory = "promoter"; inGene = thisTranscript;
                break;
            }
        } else if (strand == "-") {
            if (SNPlocus > geneEnd && SNPlocus <= (geneEnd+3000)) {
                SNPcategory = "promoter"; inGene = thisTranscript;
                break;
            }
        }
        if (SNPlocus >= geneStart && SNPlocus <= geneEnd) {
            int numDots = (int)std::count(thisTranscript.begin(), thisTranscript.end(), '.');
            SNPcategory = "intron";
            if (numDots == 4)  inGene = geneFromTranscript(thisTranscript);
            else inGene = thisTranscript;
            std::vector<string> exons = annotationMapTranscriptMap[SNPscaffold][thisTranscript];
            // std::cerr << exons.size() << std::endl;
            for (std::vector<string>::size_type j = 0; j != exons.size(); j++) {
                std::vector<string> exonVec = split(exons[j], '\t');
                //if (SNPlocus == 20001) { print_vector_stream(exonVec, std::cerr); }
                if (SNPlocus >= atoi(exonVec[1].c_str()) && SNPlocus <= atoi(exonVec[2].c_str())) {
                    SNPcategory = "exon";
                }
                if (SNPcategory == "exon") { break; }
            }
            break;
        }
        //else if (SNPcategory == "intron") { if (inGene != geneFromTranscript(startEndVec[0])) { break; } }
    }
    SNPgeneDetails.push_back(inGene); SNPgeneDetails.push_back(SNPcategory);
    return SNPgeneDetails;
    }



void Annotation::annotateGeneStartsEnds() {
    for (std::map<std::string, std::vector<std::vector<std::string> > >::iterator it = annotationMap.begin(); it != annotationMap.end(); it++) {
        std::vector<std::vector<std::string> > thisScaffoldAnnotation = it->second;
        std::vector<std::string> thisScaffoldTranscriptStartEnd;
        std::map<std::string, std::vector<std::string> > thisScaffoldTranscriptMap;
        for (std::vector<std::vector<std::string> >::size_type i = 0; i != thisScaffoldAnnotation.size(); i++) {
            std::vector<string> annotLineVec = split(thisScaffoldAnnotation[i][0], '\t');
            //if (i == 0) { print_vector_stream(annotLineVec, std::cerr); }
            string transcriptName = annotLineVec[4]; string strand = annotLineVec[3];
            string transcriptStart; string transcriptEnd;
            if (strand == "+") { transcriptStart = annotLineVec[1]; } else { transcriptEnd = annotLineVec[2]; }
            annotLineVec = split(thisScaffoldAnnotation[i].back(), '\t');
            if (strand == "+") { transcriptEnd = annotLineVec[2]; } else { transcriptStart = annotLineVec[1]; }
            // if (i == 0) { print_vector_stream(annotLineVec, std::cerr); }
            string transcriptStartEnd = transcriptName + "\t" + transcriptStart + "\t" + transcriptEnd + "\t" + strand;
            thisScaffoldTranscriptStartEnd.push_back(transcriptStartEnd);
            thisScaffoldTranscriptMap[transcriptName] = thisScaffoldAnnotation[i];
            //std::cerr << "Annotation processed: " << transcriptName << std::endl;
        }
        //std::cerr << "Annotation processed: " << it->first << std::endl;
        annotationMapTranscriptMap[it->first] = thisScaffoldTranscriptMap;
        transcriptStartEndMap[it->first] = thisScaffoldTranscriptStartEnd;
        thisScaffoldTranscriptMap.clear();
        thisScaffoldTranscriptStartEnd.clear();
        }
}
