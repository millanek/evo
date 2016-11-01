//
//  process_vcf_annotation_tools.h
//  vcf_process
//
//  Created by Milan Malinsky on 16/09/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#ifndef __vcf_process__process_vcf_annotation_tools__
#define __vcf_process__process_vcf_annotation_tools__

#include <iostream>
#include "process_vcf_utils.h"
#include "process_vcf_IUPAC.h"


inline std::string geneFromTranscript(const std::string& transcript) {
    std::vector<std::string> transcriptVec = split(transcript,'.');
    std::string geneName = transcriptVec[0] + "." + transcriptVec[1] + "." + transcriptVec[2] + "." + transcriptVec[3];
    return geneName;
}

class Annotation {
public:
    Annotation() {};
    
    Annotation(std::ifstream*& geneFile, bool includePartial) {
        annotationMap = loadAnnotationMap(geneFile, includePartial);
        getWgGeneTranscriptCounts();
        annotateGeneStartsEnds();
    }
    
    std::map<std::string, std::vector<std::vector<std::string> > > annotationMap;
    std::map<string, int> geneTranscriptCounts;  // How many transcripts has each gene
    
    string getCategoryOfSNP(const string& SNPscaffold, const string& SNPlocus) {
        string SNPcategory = "other non-coding";
        std::vector<std::vector<string> > scaffoldAnnotation = annotationMap[SNPscaffold];
        std::vector<string> scaffoldTranscriptStartEnd = transcriptStartEndMap[SNPscaffold];
        string inGene = "";
        for (std::vector<std::vector<string> >::size_type i = 0; i != scaffoldTranscriptStartEnd.size(); i++) {
            std::vector<string> startEndVec = split(scaffoldTranscriptStartEnd[i], '\t');
            if (SNPlocus >= startEndVec[1] && SNPlocus <= startEndVec[2]) {
                string thisTranscript = startEndVec[0];
                SNPcategory = "intron"; inGene = geneFromTranscript(thisTranscript);
                std::vector<string> exons = annotationMapTranscriptMap[SNPscaffold][thisTranscript];
                for (std::vector<string>::size_type j = 0; j != exons.size(); j++) {
                    std::vector<string> exonVec = split(exons[j], '\t');
                    if (SNPlocus >= exonVec[1] && SNPlocus <= exonVec[2]) {
                        SNPcategory = "exon";
                    }
                    break;
                }
                if (SNPcategory == "exon") { break; }
            } //else if (SNPcategory == "intron") { if (inGene != geneFromTranscript(startEndVec[0])) { break; } }
        }
        return SNPcategory;
    }
    
    int getTranscriptCount(const std::string& geneOrTranscript) {
        int numTranscripts = 0;
        int numDots = (int)std::count(geneOrTranscript.begin(), geneOrTranscript.end(), '.');
        if (numDots == 3)
            numTranscripts = geneTranscriptCounts[geneOrTranscript];
        else if (numDots == 4) {
            //std::vector<std::string> transcriptVec = split(geneOrTranscript,'.');
            //std::string geneName = transcriptVec[0] + "." + transcriptVec[1] + "." + transcriptVec[2] + "." + transcriptVec[3];
            std::string geneName = geneFromTranscript(geneOrTranscript);
            numTranscripts = geneTranscriptCounts[geneName];
        } else {
            std::cerr << "NumDots: " << numDots << std::endl;
            std::cerr << "geneOrTranscript: " << geneOrTranscript << std::endl;
            assert(false);
        }
        return numTranscripts;
    }
    
private:
    std::map<std::string, std::vector<std::string> > transcriptStartEndMap; // Start and end positions of each transcript
    std::map<std::string, std::map<std::string, std::vector<std::string> > > annotationMapTranscriptMap; // Quickly find the exons of any transcript
    
    
    std::string getTranscriptName(const string& geneColumn, bool& bPartial) {
        string thisTranscriptName;
        std::vector<string> findIfPartial = split(geneColumn, '-');
        if (findIfPartial.size() == 2) {
            bPartial = true;
            thisTranscriptName = findIfPartial[1];
        } else if (findIfPartial.size() == 1) {
            thisTranscriptName = geneColumn;
        }
        return thisTranscriptName;
    }
    
    // Load up the annotation file
    std::map<string, std::vector<std::vector<string> > > loadAnnotationMap(std::ifstream*& geneFile, const bool bUsePartial) {
        std::map<string, std::vector<std::vector<string> > > annotationMap;
        std::vector<std::vector<string> > annotation;
        string line;
        std::vector<string> currentTranscript;
        getline(*geneFile, line);
        currentTranscript.push_back(line);
        std::vector<string> annotLineVec = split(line, '\t');
        std::vector<string> findIfPartial = split(annotLineVec[4], '-');
        string currentScaffold = annotLineVec[0];
        bool bThisTranscriptPartial; string currentTranscriptName = getTranscriptName(annotLineVec[4], bThisTranscriptPartial);
        while (getline(*geneFile, line)) {
            std::vector<string> annotLineVec = split(line, '\t');
            bool bCurrentLinePartial = false; string currentLineTranscriptName = getTranscriptName(annotLineVec[4], bCurrentLinePartial);
            if (bCurrentLinePartial)
                bThisTranscriptPartial = true;
            if (annotLineVec[0] == currentScaffold) {
                if (currentLineTranscriptName == currentTranscriptName) {
                    currentTranscript.push_back(line);
                } else {
                    if (!bThisTranscriptPartial || bUsePartial)
                        annotation.push_back(currentTranscript);
                    currentTranscript.clear();
                    currentTranscript.push_back(line);
                    currentTranscriptName = currentLineTranscriptName;
                    bThisTranscriptPartial = bCurrentLinePartial;
                }
            } else {
                if (!bThisTranscriptPartial || bUsePartial)
                    annotation.push_back(currentTranscript);
                annotationMap[currentScaffold] = annotation;
                annotation.clear();
                currentTranscript.clear();
                currentTranscript.push_back(line);
                currentScaffold = annotLineVec[0];
                currentTranscriptName = currentLineTranscriptName;
                bThisTranscriptPartial = bCurrentLinePartial;
            }
        }
        return annotationMap;
    }
    
    void getAnnotationPerGeneDetails(const std::vector<std::vector<string> >& scaffoldAnnotation) {
        int thisGeneTranscriptCount;
        string previousGeneName = "";
        for (std::vector<string>::size_type i = 0; i != scaffoldAnnotation.size(); i++) {
            std::vector<string> annotLineVec = split(scaffoldAnnotation[i][0], '\t');
            std::string geneName = geneFromTranscript(annotLineVec[4]);
            //if (geneName == "mz.mrna.s0.284")
                // std::cerr << "geneName: " << geneName << std::endl;
            if (previousGeneName == geneName) {
                thisGeneTranscriptCount++;
            } else {
                geneTranscriptCounts[previousGeneName] = thisGeneTranscriptCount;
                thisGeneTranscriptCount = 1;
                previousGeneName = geneName;
            }
            if (i == (scaffoldAnnotation.size()-1)) {   // Add the final gene in the scaffold
                geneTranscriptCounts[previousGeneName] = thisGeneTranscriptCount;
            }
        }
    }
    
    
    void getWgGeneTranscriptCounts() {
        for (std::map<std::string, std::vector<std::vector<std::string> > >::const_iterator it = annotationMap.begin(); it != annotationMap.end(); it++) {
            getAnnotationPerGeneDetails(it->second);
        }
    }
    
    void annotateGeneStartsEnds() {
        for (std::map<std::string, std::vector<std::vector<std::string> > >::iterator it = annotationMap.begin(); it != annotationMap.end(); it++) {
            std::vector<std::vector<std::string> > thisScaffoldAnnotation = it->second;
            std::vector<std::string> thisScaffoldTranscriptStartEnd;
            std::map<std::string, std::vector<std::string> > thisScaffoldTranscriptMap;
            for (std::vector<std::vector<std::string> >::size_type i = 0; i != thisScaffoldAnnotation.size(); i++) {
                std::vector<string> annotLineVec = split(thisScaffoldAnnotation[i][0], '\t');
                string transcriptName = annotLineVec[4]; string transcriptStart = annotLineVec[1];
                annotLineVec = split(thisScaffoldAnnotation[i][thisScaffoldAnnotation[i].size()-1], '\t');
                string transcriptEnd = annotLineVec[2];
                string transcriptStartEnd = transcriptName + "\t" + transcriptStart + "\t" + transcriptEnd;
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
    
    
};


// Get the sequence for a particular individual, for a given region of the scaffold, as defined by the annotation vector
std::string getIndividualSequenceForThisRegion(const std::vector<std::string>& thisRegionAnnotation, const std::string& strand, const std::string& currentIndividualWholeScaffoldSequence);
// Get the reference sequence for a given region of the scaffold defined by the annotation vector
std::string getReferenceForThisRegion(const std::vector<std::string>& thisRegionAnnotation, const std::string& strand, const std::string& currentScaffoldReference);


class AccessibleGenome {
public:
    AccessibleGenome() {};
    
    AccessibleGenome(std::ifstream*& bedFile) {
        accessibleGenomeMap = loadAccessibleGenomeMap(bedFile);
    }
    std::map<std::string, std::vector<std::vector<int> > > accessibleGenomeMap;
    
    int getAccessibleBPinRegion(const string& scaffold, const int start, const int end) {
        
        std::vector<std::vector<int> > aGThisSc = accessibleGenomeMap[scaffold];
        //std::cerr << "There are " << aGThisSc[0].size() << " accessible intervals" << std::endl;
        /*for (std::map<std::string, std::vector<std::vector<int> > >::iterator i = accessibleGenomeMap.begin(); i != accessibleGenomeMap.end(); i++) {
            std::cerr << "There is scaffold: " << i->first << " in the map" << std::endl;
        } */
        int numBP = 0;
        
        // Binary search to find the first element in the accessible genome annotation whose end coordinate is greater
        // or equal to the start of the region in question
        std::vector<int>::iterator itStart = lower_bound(aGThisSc[1].begin(),aGThisSc[1].end(),start);
        std::vector<int>::size_type index = std::distance(aGThisSc[1].begin(), itStart);
        // Add the first element:
        numBP = numBP + (aGThisSc[1][index] - start) + 1; index++;
        // Sum the lengths
        while (aGThisSc[0][index] <= end) {
            if (aGThisSc[1][index] < end)
                numBP = numBP + (aGThisSc[1][index] - aGThisSc[0][index]);
            else
                numBP = numBP + (end - aGThisSc[0][index]);
            index++;
        }
        return numBP;
    }
    
    bool findIfBPaccessible(const string& scaffold, const int bp) {
        std::vector<std::vector<int> > aGThisSc = accessibleGenomeMap[scaffold];
        
        // Binary search to find the first element in the accessible genome annotation whose end coordinate is greater
        // or equal to the basepair in question
        std::vector<int>::iterator itStart = lower_bound(aGThisSc[1].begin(),aGThisSc[1].end(),bp);
        std::vector<int>::size_type index = std::distance(aGThisSc[1].begin(), itStart);
        
        
        if (aGThisSc[0][index] <= bp) {
            return false;
        } else {
            return true;
        }
    }
    
    string getAccessibleSeqForScaffold(const string& scaffold, const string& fullString) {
        std::vector<std::vector<int> > aGThisSc = accessibleGenomeMap[scaffold];
        
        string accessibleString; accessibleString.reserve(fullString.length()); accessibleString = "";
        for (std::vector<int>::size_type i = 0; i < aGThisSc[0].size(); i++) {
            accessibleString.append(fullString.substr(aGThisSc[0][i],aGThisSc[1][i]- aGThisSc[0][i]));
        }
        return accessibleString;
    }
    
private:
    // Load up the file specifying the accessible genome (needs to be sorted by chromosome)
    std::map<string, std::vector<std::vector<int> > > loadAccessibleGenomeMap(std::ifstream*& bedFile) {
        std::map<string, std::vector<std::vector<int> > > accessibleGenomeMap;
        std::vector<std::vector<int> > accessibleGenomeThisScaffold;
        std::vector<int> featureStarts;
        std::vector<int> featureEnds;
        string line;
        getline(*bedFile, line);
        std::vector<string> currentFeature = split(line, '\t');
        string currentScaffold = currentFeature[0];
        //std::cerr << "Loading scaffold: " << currentScaffold << std::endl;
        while (getline(*bedFile, line)) {
            featureStarts.push_back(atoi(currentFeature[1].c_str()));
            featureEnds.push_back(atoi(currentFeature[2].c_str()));
           // std::cerr << "Loading scaffold: " << currentFeature[0] << "\t" << currentScaffold << std::endl;
            if (currentFeature[0] != currentScaffold) {
                accessibleGenomeThisScaffold.push_back(featureStarts);
                accessibleGenomeThisScaffold.push_back(featureEnds);
                accessibleGenomeMap[currentScaffold] = accessibleGenomeThisScaffold;
                accessibleGenomeThisScaffold.clear(); featureStarts.clear(); featureEnds.clear();
                currentScaffold = currentFeature[0];
          //      std::cerr << "Loading scaffold: " << currentScaffold << std::endl;
            }
            currentFeature = split(line, '\t');
            //std::cerr << "Loading scaffold: " << currentFeature[0] << "\t" << currentScaffold << std::endl;
        }
        // Final line / final scaffold
        accessibleGenomeThisScaffold.push_back(featureStarts);
        accessibleGenomeThisScaffold.push_back(featureEnds);
        accessibleGenomeMap[currentScaffold] = accessibleGenomeThisScaffold;
        return accessibleGenomeMap;
    }

};



#endif /* defined(__vcf_process__process_vcf_annotation_tools__) */
