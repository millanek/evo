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


class SimpleCoordsBed {
public:
    SimpleCoordsBed() {};

    SimpleCoordsBed(std::ifstream*& bedFile, std::map<int, string>& linearToGenomeCoordMap) {
        loadSimpleCoords(bedFile, linearToGenomeCoordMap);
    }
    
private:
    // Require: coordinates are sorted by scaffold
    std::map<string, std::vector<std::vector<string> > > loadSimpleCoords(std::ifstream*& bedFile, std::map<int,string>& linearToGenomeCoordMap) {
        std::map<string, std::vector<std::vector<string> > > coordsScaffoldMap;
        std::vector<std::vector<string> > coordsInScaffold;
        int linearPosition = 0; // Linear positions will be 0-indexed (will start at zero)
        string line;
        // Do the first line to find the name of the first scaffold
        getline(*bedFile, line);
        std::vector<string> bedVector = split(line, '\t');
        string currentScaffold = bedVector[0];
        int left = atoi(bedVector[1].c_str());
        int right = atoi(bedVector[2].c_str());
        int intervalLength = right - left;
        coordsInScaffold.push_back(bedVector);
        
        for (int i = 0; i < intervalLength; i++) {
            // Mapping will be 1-indexed, as in VCF files
            linearToGenomeCoordMap[linearPosition+i] = currentScaffold+"\t"+numToString(left+i+1);
        } linearPosition = linearPosition + intervalLength;
        
        //
        while (getline(*bedFile, line)) {
            std::vector<string> bedVector = split(line, '\t');
            string scaffold = bedVector[0];
            int left = atoi(bedVector[1].c_str());
            int right = atoi(bedVector[2].c_str());
            int intervalLength = right - left;
            
            for (int i = 0; i < intervalLength; i++) {
                // Mapping will be 1-indexed, as in VCF files
                int pos = linearPosition+i;
                if (pos % 1000000 == 0)
                    std::cerr << "Loaded and mapped " << pos/1000000 << "Mb" << std::endl;
                linearToGenomeCoordMap[pos] = currentScaffold+"\t"+numToString(left+i+1);
            } linearPosition = linearPosition + intervalLength;
            
            if (scaffold == currentScaffold) {
                coordsInScaffold.push_back(bedVector);
            } else {
                coordsScaffoldMap[currentScaffold] = coordsInScaffold;
                coordsInScaffold.clear();
                currentScaffold = scaffold;
            }
        }
        return coordsScaffoldMap;
    }
};

class LinkedCoordsBed {
public:
    LinkedCoordsBed() {};
    
    LinkedCoordsBed(std::ifstream*& bedFile) {
        elementCoordsVector = loadLinkedCoords(bedFile);
    }
    
    std::vector<std::vector<std::vector<string> > > elementCoordsVector;
    
    
    std::vector<std::vector<string> > getElementOuterBoundaries() {
        std::vector<std::vector<string> > allOuterBounds;
        for (int i = 0; i < (int)elementCoordsVector.size(); i++) {
            std::vector<std::vector<string> > thisElement = elementCoordsVector[i];
            //if (thisElement[0].size() < 4) {
            //    std::cerr << "thisElement[0].size()" << thisElement.size() << std::endl;
            //}
            std::vector<string> elementOuterBounds;
            elementOuterBounds.push_back(thisElement[0][0]);
            elementOuterBounds.push_back(thisElement[0][1]);
            elementOuterBounds.push_back(thisElement.back()[2]);
            allOuterBounds.push_back(elementOuterBounds);
        }
        return allOuterBounds;
    }
    
    std::vector<string> getElementNames() {
        std::vector<string> elementNames;
        for (int i = 0; i < (int)elementCoordsVector.size(); i++) {
            std::vector<std::vector<string> > thisElement = elementCoordsVector[i];
            //if (thisElement[0].size() < 4) {
            //    std::cerr << "thisElement[0].size()" << thisElement.size() << std::endl;
            //}
            elementNames.push_back(thisElement[0][3]);
        }
        return elementNames;
    }
    
    std::vector<double> getDxyPerElement(std::unordered_map<string, double>& posDxyMap) {
        std::vector<double> elementDxyValues;
        for (int i = 0; i < (int)elementCoordsVector.size(); i++) {
            std::vector<std::vector<string> > thisElement = elementCoordsVector[i];
            double elementDxyTotal = 0; int elementLengthTotal = 0;
            for (int j = 0; j < (int)thisElement.size(); j++) {
                std::vector<string> thisCoordVector = thisElement[j];
                string currentScaffold = thisCoordVector[0];
                int left = atoi(thisCoordVector[1].c_str());
                int right = atoi(thisCoordVector[2].c_str());
                int intervalLength = right - left;
                for (int k = 0; k < intervalLength; k++) {
                    string genomeLocation = currentScaffold+"\t"+numToString(left+k+1);
                    if (posDxyMap.count(genomeLocation) == 1) {
                        elementDxyTotal = elementDxyTotal + posDxyMap[genomeLocation];
                    }
                }
                elementLengthTotal = elementLengthTotal + intervalLength;
            }
            elementDxyValues.push_back(elementDxyTotal/elementLengthTotal);
        }
        return elementDxyValues;
    }
    
private:
    // Require: coordinates are sorted by scaffold
    std::vector<std::vector<std::vector<string> > > loadLinkedCoords(std::ifstream*& bedFile) {
        std::vector<std::vector<std::vector<string> > > elementCoordsVector;
        std::vector<std::vector<string> > thisElement;
        
        // Do the first line to find the name of the first scaffold and of the first element
        string line;
        getline(*bedFile, line);
        std::vector<string> bedVector = split(line, '\t');
        string currentElementName = bedVector[3];
        thisElement.push_back(bedVector);
        
        //
        while (getline(*bedFile, line)) {
            std::vector<string> bedVector = split(line, '\t');
            string element = bedVector[3];
            if (element == currentElementName) {
                thisElement.push_back(bedVector);
            } else {
                elementCoordsVector.push_back(thisElement);
                thisElement.clear();
                currentElementName = element;
                thisElement.push_back(bedVector);
            }
        }
        return elementCoordsVector;
    }
};



    
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
                int numDots = (int)std::count(thisTranscript.begin(), thisTranscript.end(), '.');
                SNPcategory = "intron";
                if (numDots == 4)  inGene = geneFromTranscript(thisTranscript);
                else inGene = thisTranscript;
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
    
    
    std::vector<string> getSNPgeneDetails(const string& SNPscaffold, const string& SNPlocus) {
        std::vector<string> SNPgeneDetails;
        std::vector<string> scaffoldTranscriptStartEnd = transcriptStartEndMap[SNPscaffold];
        string inGene = ""; string SNPcategory = "nonCoding";
        for (std::vector<std::vector<string> >::size_type i = 0; i != scaffoldTranscriptStartEnd.size(); i++) {
            std::vector<string> startEndVec = split(scaffoldTranscriptStartEnd[i], '\t');
            if (SNPlocus >= startEndVec[1] && SNPlocus <= startEndVec[2]) {
                string thisTranscript = startEndVec[0];
                SNPcategory = "intron";
                int numDots = (int)std::count(thisTranscript.begin(), thisTranscript.end(), '.');
                SNPcategory = "intron";
                if (numDots == 4)  inGene = geneFromTranscript(thisTranscript);
                else inGene = thisTranscript;
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
        SNPgeneDetails.push_back(inGene); SNPgeneDetails.push_back(SNPcategory);
        return SNPgeneDetails;
    }
    
    
    int getTranscriptCount(const std::string& geneOrTranscript) {
        int numTranscripts = 0;
        int numDots = (int)std::count(geneOrTranscript.begin(), geneOrTranscript.end(), '.');
        if (numDots < 4)
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
            int numDots = (int)std::count(annotLineVec[4].begin(), annotLineVec[4].end(), '.');
            std::string geneName;
            if (numDots == 4)  geneName = geneFromTranscript(annotLineVec[4]);
            else geneName = annotLineVec[4];
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





class BedCoordinateFeatures {
public:
    BedCoordinateFeatures() : initialised(false) {};
    
    BedCoordinateFeatures(std::ifstream*& bedFile, bool debug = false) {
        BedFeatureMap = loadBedFeatureMap(bedFile, debug);
        initialised = true;
    }
    bool initialised;

    
    // Count how much of a genomic interval is covered by the bed features (returns the number of basepairs)
    // Left coordinates of features are 0-indexed (bed file)
    // Query is 1-indexed
    // query:              -----------------------
    // feature: 1)  ---                                   (f[1] <= start)
    //          2)                                  ---   (f[0] >= end)
    //          3)     ------                             (f[0] < start && f[1] <= end)
    //          4)              --------                  (f[0] >=start && f[1] <= end)
    //          5)                             -------    (f[0] >=start && f[1] > end)
    //          6)       -----------------------------    (f[0] < start && f[1] > end)
    int getNumBPinRegion(const string& scaffold, const int start, const int end) {
        assert(start < end);
        try {
            std::vector<std::vector<int> > featuresThisSc = BedFeatureMap.at(scaffold);
            /*for (std::map<std::string, std::vector<std::vector<int> > >::iterator i = accessibleGenomeMap.begin(); i != accessibleGenomeMap.end(); i++) {
             std::cerr << "There is scaffold: " << i->first << " in the map" << std::endl;
             } */
            
            
            // Binary search to find the first feature whose end coordinate is greater
            // or equal to the start of the region in question
            std::vector<int>::iterator itStart = lower_bound(featuresThisSc[1].begin(),featuresThisSc[1].end(),start);
            int numBP = 0;
            if (itStart != featuresThisSc[1].end()) {  // if (start < f[1])    ---  excludind case 1)
                std::vector<int>::size_type index = std::distance(featuresThisSc[1].begin(), itStart);
                
                // Sum the lengths
                while (featuresThisSc[0][index] <= end && index < featuresThisSc[0].size()) { // if (f[0] >= end)   ---    ecluding case 2)
                    //std::cerr << "There are " << featuresThisSc[0].size() << " intervals" << std::endl;
                    //std::cerr << "Starts: " << std::endl;
                    //print_vector(featuresThisSc[0], std::cerr);
                    //std::cerr << "Ends: " << std::endl;
                    //print_vector(featuresThisSc[1], std::cerr);
                    // Now we know the overlap is > 0
                    //std::cerr << "start, end, featuresThisSc[0][index], featuresThisSc[1][index]: " << start << ", " << end << ", " << featuresThisSc[0][index] << ", " << featuresThisSc[1][index] << std::endl;
                    
                    if (featuresThisSc[0][index] < start && featuresThisSc[1][index] <= end)
                        numBP = numBP + (featuresThisSc[1][index] - start) + 1;
                    else if (featuresThisSc[0][index] >= start && featuresThisSc[1][index] <= end)
                        numBP = numBP + (featuresThisSc[1][index] - featuresThisSc[0][index]);
                    else if (featuresThisSc[0][index] >= start && featuresThisSc[1][index] > end)
                        numBP = numBP + (end - featuresThisSc[0][index]);
                    else if (featuresThisSc[0][index] < start && featuresThisSc[1][index] > end)
                        numBP = numBP + (end - start) + 1;
                    index++;
                }
            }
            return numBP;
        } catch (const std::out_of_range& oor) {
            //std::cerr << "No features on scaffold: " << scaffold << std::endl;
            return 0;
        }
    }
    
    std::vector<std::vector <string> > getFeaturesinRegion(const string& scaffold, const int start, const int end) {
        std::vector<std::vector <string> > allFeatures;
        try {
            std::vector<std::vector<int> > featuresThisSc = BedFeatureMap.at(scaffold);
            //std::cerr << "There are " << aGThisSc[0].size() << " accessible intervals" << std::endl;
            /*for (std::map<std::string, std::vector<std::vector<int> > >::iterator i = accessibleGenomeMap.begin(); i != accessibleGenomeMap.end(); i++) {
             std::cerr << "There is scaffold: " << i->first << " in the map" << std::endl;
             } */
            
            
            // Binary search to find the first feature whose end coordinate is greater
            // or equal to the start of the region in question
            std::vector<int>::iterator itStart = lower_bound(featuresThisSc[1].begin(),featuresThisSc[1].end(),start);
            std::vector<int>::size_type index = std::distance(featuresThisSc[1].begin(), itStart);
            // Now check that the start coordinate of the feature is within
            
            // Sum the lengths
            while (featuresThisSc[0][index] <= end) {
                std::vector <string> feature;
                feature.push_back(scaffold); feature.push_back(numToString(featuresThisSc[0][index]));
                feature.push_back(numToString(featuresThisSc[1][index]));
                index++;
            }
            return allFeatures;
        } catch (const std::out_of_range& oor) {
            std::cerr << "No features on scaffold: " << scaffold << std::endl;
            return allFeatures;
        }
    }
    
    
    
    // In the bed file, the left coordinate is zero indexed
    // bp is a 1-indexed coordinate
    bool findIfBPinBedFile(const string& scaffold, const int bp) {
        try {
            std::vector<std::vector<int> > featuresThisSc = BedFeatureMap.at(scaffold);
            
            // Binary search to find the first element in the accessible genome annotation whose end coordinate is greater
            // or equal to the basepair in question
            std::vector<int>::iterator itStart = lower_bound(featuresThisSc[1].begin(),featuresThisSc[1].end(),bp);
            std::vector<int>::size_type index = std::distance(featuresThisSc[1].begin(), itStart);
            
            // Then if the start coordinate is less than the coordinate of the basepair, the bp is covered by the bed file features
            if (featuresThisSc[0][index] < bp) {
                return true;
            } else {
                return false;
            }
        } catch (const std::out_of_range& oor) {
            return false;
        }
    }
    
protected:
    std::map<std::string, std::vector<std::vector<int> > > BedFeatureMap;
    
    // Load up the file specifying the genomic coordinates (the bed file needs to be sorted by chromosome)
    std::map<string, std::vector<std::vector<int> > > loadBedFeatureMap(std::ifstream*& bedFile, bool debug) {
        std::map<string, std::vector<std::vector<int> > > BedFeatureMap;
        std::vector<std::vector<int> > BedFeaturesThisScaffold;
        std::vector<int> featureStarts;
        std::vector<int> featureEnds;
        string line;
        string previousScaffold = "";
        //std::cerr << "Loading scaffold: " << currentScaffold << std::endl;
        while (getline(*bedFile, line)) {
            std::vector<string> currentFeature = split(line, '\t');
            if (currentFeature[0] != previousScaffold && previousScaffold != "") {
                if (debug) {
                    std::cerr << "Loading scaffold: " << previousScaffold << std::endl;
                    std::cerr << "Starts: " << std::endl;
                    print_vector(featureStarts, std::cerr);
                    std::cerr << "Ends: " << std::endl;
                    print_vector(featureEnds, std::cerr);
                }
                BedFeaturesThisScaffold.push_back(featureStarts);
                BedFeaturesThisScaffold.push_back(featureEnds);
                BedFeatureMap[previousScaffold] = BedFeaturesThisScaffold;
                BedFeaturesThisScaffold.clear(); featureStarts.clear(); featureEnds.clear();
                // std::cerr << "Loading scaffold: " << currentScaffold << std::endl;
            }
            featureStarts.push_back(atoi(currentFeature[1].c_str()));
            featureEnds.push_back(atoi(currentFeature[2].c_str()));
            // std::cerr << "Loading scaffold: " << currentFeature[0] << "\t" << currentScaffold << std::endl;
            previousScaffold = currentFeature[0];
            //std::cerr << "Loading scaffold: " << currentFeature[0] << "\t" << currentScaffold << std::endl;
        }
        // Final line / final scaffold
        BedFeaturesThisScaffold.push_back(featureStarts);
        BedFeaturesThisScaffold.push_back(featureEnds);
        BedFeatureMap[previousScaffold] = BedFeaturesThisScaffold;
        return BedFeatureMap;
    }
    
};



class AccessibleGenome : public BedCoordinateFeatures {
public:
    AccessibleGenome() {};
    
    AccessibleGenome(std::ifstream*& bedFile) {
        BedFeatureMap = loadBedFeatureMap(bedFile, false);
    }
    
    int getAccessibleBPinRegion(const string& scaffold, const int start, const int end) {
        return getNumBPinRegion(scaffold, start, end);
    }
    
    // Accessible genome is a bed file (left coordinate is zero indexed)
    // bp is a 1-indexed coordinate
    bool findIfBPaccessible(const string& scaffold, const int bp) {
        return findIfBPinBedFile(scaffold,bp);
    }
    
    string getAccessibleSeqForScaffold(const string& scaffold, const string& fullString) {
        std::vector<std::vector<int> > aGThisSc = BedFeatureMap[scaffold];
        
        string accessibleString; accessibleString.reserve(fullString.length()); accessibleString = "";
        for (std::vector<int>::size_type i = 0; i < aGThisSc[0].size(); i++) {
            accessibleString.append(fullString.substr(aGThisSc[0][i],aGThisSc[1][i]- aGThisSc[0][i]));
        }
        return accessibleString;
    }

};



#endif /* defined(__vcf_process__process_vcf_annotation_tools__) */
