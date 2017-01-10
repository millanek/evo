//
//  process_vcf_coding_sequences.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 16/09/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include "process_vcf_coding_sequences.h"

#define SUBPROGRAM "getCodingSeq"

#define DEBUG 1

static const char *CODINGSEQ_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] VCF_FILE GENOME_SEQUENCE.fa ANNOTATION.gffExtract\n"
"Obtain full genome sequences from a VCF file (e.g. for multiple alignment and phylogenetic analyses), output to STD_OUT\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -n, --non-coding                        the ANNOTATION.gffExtract file does not specify coding sequences\n"
"                                               don't do any stats or analyses intended for coding sequences\n"
"       -p, --coding-partial                    include genes whose annotation includes CDS overlapping with UTRs\n"
"                                               i.e. the precise start and or end of coding sequence is unknown\n"
"                                               ('5_prime_partial=true' or '3_prime_partial=true' in gff3 CDS line)\n"
"                                               (NOT COMPATIBLE WITH THE -n OPTION)\n"
"       --only-stats                            only calculate statistics for coding sequences\n"
"                                               do not output the sequences themselves\n"
"       -H,   --het-treatment <r|p|i>           r: assign het bases randomly (default); p: use the phase information in a VCF outputting haplotype 1 for each individual; i: use IUPAC codes\n"
"       -s SAMPLES.txt, --samples=SAMPLES.txt   supply a file of sample identifiers to be used for the output\n"
"                                               (default: sample ids from the vcf file are used)\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hpws:cH:";
static std::vector<string> stopsTranscriptRecord;

enum { OPT_ONLY_STATS, OPT_USE_PHASE };


static const struct option longopts[] = {
    { "non-coding",   required_argument, NULL, 'n' },
    { "samples",   required_argument, NULL, 's' },
    { "partial",   no_argument, NULL, 'p' },
    { "het-treatment",   required_argument, NULL, 'H' },
    { "only-stats",   no_argument, NULL, OPT_ONLY_STATS },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string genomeFile;
    static string geneFile;
    static bool bIsCoding = true;
    static bool bUsePartial = false;
    static bool bOnlyStats = false;
    static char hetTreatment = 'r';
    static string sampleNameFile;
    
}


void printPerGeneSummaries(std::ofstream*& stopsPerGeneSummaryFile, Annotation& wgAnnotation) {
    std::string previousGene = "";
    int numTranscripsWithStopsThisGene = 0;
    double stopAlleleFrequencySum = 0;
    std::vector<double> stopPosPctVec;
    for (std::vector<string>::size_type i = 0; i != stopsTranscriptRecord.size(); i++) {
         std::cerr << "previousGene: " << previousGene << std::endl;
        std::vector<string> oneTransctipt = split(stopsTranscriptRecord[i], '\t');
        if (geneFromTranscript(oneTransctipt[0]) == previousGene) {
            numTranscripsWithStopsThisGene++;
            stopAlleleFrequencySum += convertToDouble(oneTransctipt[3]);
            double stopPosPct = convertToDouble(oneTransctipt[1]) / convertToDouble(oneTransctipt[2]);
            stopPosPctVec.push_back(stopPosPct);
        } else {
            if (previousGene != "")
                *stopsPerGeneSummaryFile << previousGene << "\t" << numTranscripsWithStopsThisGene << "\t" << wgAnnotation.getTranscriptCount(previousGene) << "\t" << stopAlleleFrequencySum/numTranscripsWithStopsThisGene << "\t" << vector_average(stopPosPctVec) << std::endl;
            previousGene = geneFromTranscript(oneTransctipt[0]);
            numTranscripsWithStopsThisGene = 1; stopAlleleFrequencySum = convertToDouble(oneTransctipt[3]); stopPosPctVec.clear();
            double stopPosPct = convertToDouble(oneTransctipt[1]) / convertToDouble(oneTransctipt[2]);
            stopPosPctVec.push_back(stopPosPct);
        }
        
        // print the final gene
        if (i == (stopsTranscriptRecord.size() - 1)) {
            *stopsPerGeneSummaryFile << geneFromTranscript(oneTransctipt[0]) << "\t" << numTranscripsWithStopsThisGene << "\t" << wgAnnotation.getTranscriptCount(previousGene) << "\t" << stopAlleleFrequencySum/numTranscripsWithStopsThisGene << "\t" << vector_average(stopPosPctVec) << std::endl;
        }
        
        //oneTransctipt[0]
    }
}

int getCodingSeqMain(int argc, char** argv) {
    parseGetCodingSeqOptions(argc, argv);
    string vcfFileRoot = stripExtension(opt::vcfFile);
    string geneFileRoot = stripExtension(opt::geneFile);
    
    if (opt::bIsCoding) {
        std::cerr << "Generating gene coding sequences using variants from: " << opt::vcfFile << std::endl;
        std::cerr << "and the reference genome: " << opt::genomeFile << std::endl;
        std::cerr << "with coding sequence annotation defined in: " << opt::geneFile << std::endl;
    } else {
        std::cerr << "Generating sequences using variants from: " << opt::vcfFile << std::endl;
        std::cerr << "and the reference genome: " << opt::genomeFile << std::endl;
        std::cerr << "extracting regions according to annotation defined in: " << opt::geneFile << std::endl;
    }
        
    // Open connections to read from the vcf and reference genome files
    std::istream* vcfFile = createReader(opt::vcfFile.c_str());
    std::ifstream* genomeFile = new std::ifstream(opt::genomeFile.c_str());
    std::ifstream* geneFile = new std::ifstream(opt::geneFile.c_str());
    
    // Prepare an object for output files
    std::ofstream* badStartStopCodonFile;
    if (opt::bIsCoding) badStartStopCodonFile = new std::ofstream("badStartStopCodonList.txt");
    
    string line;
    string currentScaffoldNum = "";
    string currentScaffoldReference;
    string::size_type inStrPos;
    size_t numSamples;
    std::vector<string> sampleNames;
    string thisScaffoldName;
    unsigned int processedVariantCounter = 0;
    std::vector<string> scaffoldStrings;
    std::vector<std::vector<string> > statsAllGenes;
    string statsFileName = geneFileRoot + "_stats.txt";
    std::ofstream* statsFile = new std::ofstream(statsFileName.c_str());
    string stopsFileName = geneFileRoot + "_prematureStops.txt";
    *statsFile << "transcript" << "\t" << "length_in_nucleotides" << "\t" << "segregating_sites(ss)" << "\t" << "ss_proportion" << "\t" << "length_in_AA" << "\t" << "num_of_AA_with_synonymous_changes(synAAs)" << "\t" << "synAAs_proportion" << "\t" << "non_synonymous_AA_substitutions(nsAAs)" << "\t" << "nsAAs_proportion" << "\t" << "synonymousMAFaverage" << "\t" << "nonsynonymousMAFaverage" << "\t" << "pN" << "\t" << "pS" << std::endl;
    std::cout << "transcript" << "\t" << "length_in_nucleotides" << "\t" << "segregating_sites(ss)" << "\t" << "ss_proportion" << "\t" << "length_in_AA" << "\t" << "num_of_AA_with_synonymous_changes(synAAs)" << "\t" << "synAAs_proportion" << "\t" << "non_synonymous_AA_substitutions(nsAAs)" << "\t" << "nsAAs_proportion" << "\t" << "synonymousMAFaverage" << "\t" << "nonsynonymousMAFaverage" << "\t" << "pN" << "\t" << "pS" << std::endl;
    std::ofstream* stopsFile = new std::ofstream(stopsFileName.c_str());
    *stopsFile << "transcript" << "\t" << "stopAA_position" << "\t" << "transcript_length" << "\t" << "stop_allele_frequency" << "\t" << "individuals_with_stop" << std::endl;
    string stopsPerGeneFileName = geneFileRoot + "_prematureStops_perGene.txt";
    std::ofstream* stopsPerGeneSummaryFile = new std::ofstream(stopsPerGeneFileName.c_str());
    *stopsPerGeneSummaryFile << "gene" << "\t" << "numStops" << "\t" << "numTranscripts" << "\t" << "avg_stop_allele_frequency" << "\t" << "avg_stop_AA_position(%_of_trancript_length)" << std::endl;
    
    // Load up the gene (CDS) annotation file
    Annotation wgAnnotation(geneFile, opt::bUsePartial);
    std::vector<std::vector<string> > annotation;
    
    //print_vector_stream(wgAnnotation.annotationMap["scaffold_0"][0], std::cerr);
    
    while (getline(*vcfFile, line)) {
        if (line[0] == '#' && line[1] == '#')
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            std::vector<std::string> fields = split(line, '\t');
            numSamples = fields.size()-NUM_NON_GENOTYPE_COLUMNS;
            // Initialize vectors
            scaffoldStrings.resize(numSamples);
            for (std::vector<string>::size_type i = 0; i != scaffoldStrings.size(); i++) {
                scaffoldStrings[i].reserve(100000000);
                scaffoldStrings[i] = "";
            }
            
            for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                if (opt::sampleNameFile.empty())
                    sampleNames.push_back(fields[i]);
                else
                    sampleNames = readSampleNamesFromTextFile(opt::sampleNameFile);
            }
            
        } else {
            processedVariantCounter++;
            std::vector<std::string> fields = split(line, '\t');
            // Also get overall depth for this variant
            std::vector<std::string> info = split(fields[7], ';');
            if (fields[0] != currentScaffoldNum) {
                if (currentScaffoldNum != "") {
                    for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                        scaffoldStrings[i].append(currentScaffoldReference.substr(inStrPos, string::npos));
                    }
                    
#ifdef DEBUG
                    if (scaffoldStrings[0].length() != currentScaffoldReference.length()) {
                        std::cerr << "Error!!! Reference scaffold length: " << currentScaffoldReference.length() << " vcf scaffold length: " << scaffoldStrings[0].length() << std::endl;
                    }
#endif
                    
                    std::cerr << currentScaffoldNum << std::endl;
                    annotation = wgAnnotation.annotationMap[currentScaffoldNum]; // Get annotation for this scaffold
                    
                    std::cerr << "Going through the annotation..." << std::endl;
                    for (std::vector<std::vector<string> >::size_type k = 0; k != annotation.size(); k++) {
                        std::vector<string> annotLineVec = split(annotation[k][0], '\t');
                        // std::cerr << "Gene:" << annotLineVec[4] << std::endl;
                        std::ofstream* geneOutFiles;
                        bool geneLengthDivisibleByThree = true;
                        std::vector<string> allSeqs;
                        std::vector<string> statsThisGene; statsThisGene.push_back(annotLineVec[4]);
                        string refSeq = getReferenceForThisRegion(annotation[k], annotLineVec[3], currentScaffoldReference);
                        if (opt::bIsCoding)
                            geneLengthDivisibleByThree = codingSequenceErrorChecks(refSeq, annotLineVec[4], annotation, (int)k, badStartStopCodonFile);
                        if (geneLengthDivisibleByThree && !opt::bOnlyStats)
                            geneOutFiles = new std::ofstream(annotLineVec[4].c_str());
                        for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                            string geneSeq = getIndividualSequenceForThisRegion(annotation[k], annotLineVec[3], scaffoldStrings[i]);
                            if (geneLengthDivisibleByThree) {
                                if (!opt::bOnlyStats) *geneOutFiles << ">" << sampleNames[i] << std::endl;
                                if (!opt::bOnlyStats) *geneOutFiles << geneSeq << std::endl;
                                allSeqs.push_back(geneSeq);
                            }
                        }
                       // std::cerr << "Got all sequences for: " << annotLineVec[4] << std::endl;
                        
                        // Get statistics for the sequences
                        if (opt::bIsCoding && geneLengthDivisibleByThree) {
                            if (opt::hetTreatment == 'p' || opt::hetTreatment == 'r') {
                                getCodingSequenceStatsPhasedSeq(allSeqs, refSeq, annotLineVec[4], statsThisGene,stopsFile, sampleNames, wgAnnotation);
                                getCodingSequencePhasedPnPs(allSeqs, refSeq, annotLineVec[4], statsThisGene,stopsFile, sampleNames, wgAnnotation);
                                print_vector_stream(statsThisGene, std::cout);
                            } else {
                                getCodingSequenceStatsIUPAC(allSeqs, refSeq, annotLineVec[4], statsThisGene,stopsFile, sampleNames, wgAnnotation);
                            }
                            statsAllGenes.push_back(statsThisGene);
                            // std::cerr << "Got statistics for: " << annotLineVec[4] << std::endl;
                        }
                        if (!opt::bOnlyStats) geneOutFiles->close();
                    }
                    for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                        scaffoldStrings[i] = "";
                    }
                    processedVariantCounter = 1;
                    currentScaffoldNum = fields[0];
                    forwardGenomeToScaffold(currentScaffoldNum, genomeFile, thisScaffoldName);
                } else {
                    getline(*genomeFile, thisScaffoldName);
                    thisScaffoldName.erase(0,1);
                    currentScaffoldNum = fields[0];
                }
                
                print_matrix(statsAllGenes, *statsFile);
                statsFile->flush();
                statsAllGenes.clear();
                printPerGeneSummaries(stopsPerGeneSummaryFile, wgAnnotation);
                stopsTranscriptRecord.clear();
                
                inStrPos = 0;
                std::cerr << "Starting to read " << thisScaffoldName << std::endl;
                currentScaffoldReference = readScaffold(genomeFile, thisScaffoldName);
                thisScaffoldName.erase(0,1);
                std::cerr << "Finished reading" << std::endl;
                std::cerr << "Generating sequences with variants from the VCF file..." << std::endl;
                
            }
            if (info[0] != "INDEL") {
                for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                    //std::cerr << "Going through genotypes1:" << i << std::endl;
                    //std::cerr << scaffoldStrings.size() << " " << inStrPos << " " << fields[1] << " " << currentScaffoldReference.size() << std::endl;
                    scaffoldStrings[i- NUM_NON_GENOTYPE_COLUMNS].append(currentScaffoldReference.substr(inStrPos, (atoi(fields[1].c_str()) - 1)-inStrPos));
                    std::vector<string> genotypeFields = split(fields[i], ':');
                    std::vector<char> genotype;
                    genotype.push_back(genotypeFields[0][0]); genotype.push_back(genotypeFields[0][2]);
                    appendGenotypeBaseToString(scaffoldStrings[i- NUM_NON_GENOTYPE_COLUMNS], fields[3], fields[4], genotype, opt::hetTreatment);
                }
                inStrPos = atoi(fields[1].c_str());
                
#ifdef DEBUG
                if (currentScaffoldReference[inStrPos-1] != fields[3][0]) {
                    std::cerr << "Error!!! Sequence: " << currentScaffoldReference[inStrPos-1] << " vcf-ref: " << fields[3][0] << std::endl;
                }
#endif
            }
            if (processedVariantCounter % 10000 == 0)
                std::cerr << processedVariantCounter << " variants processed..." << std::endl;
        }
    }
    
    print_matrix(statsAllGenes, *statsFile);
    statsAllGenes.clear();
    printPerGeneSummaries(stopsPerGeneSummaryFile, wgAnnotation);
    stopsTranscriptRecord.clear();

    
    return 0;
}

// Return false if the length of the coding sequence is not divisible by three
bool codingSequenceErrorChecks(const string& geneSeq, const string& transcriptName, const std::vector<std::vector<string> >& annotation, const int k, std::ofstream*& badStartStopCodonFile) {
    string::size_type l = geneSeq.length();
    bool divisibleByThree = true;
    if (l % 3 != 0) {
#ifdef DEBUG2
        std::cerr << "Error!!! The length of the gene: " << transcriptName << " is not divisible by three (" << l << ") composed of " << annotation[k].size() << " coding sequences" << std::endl;
        for (std::vector<std::string>::size_type j = 0; j != annotation[k].size(); j++) {
            std::cerr << annotation[k][j] << std::endl;
        }
#endif
        divisibleByThree = false;
    }
    if (!(geneSeq[0] == 'A' && geneSeq[1] == 'T' && geneSeq[2] == 'G')) {
        std::cerr << "Possible error!!! The coding sequence of the gene: " << transcriptName << " does not start with ATG (start codon). Maybe a SNP? Or an annotation error..." << std::endl;
        *badStartStopCodonFile << transcriptName << "\t" << "start" << "\t" << geneSeq.substr(0,3) << std::endl;
    }
    if (geneSeq[l-3] != 'T' && ((geneSeq[l-2] == 'A' && geneSeq[l-1] == 'G') || (geneSeq[l-2] == 'A' && geneSeq[l-1] == 'A') || (geneSeq[l-2] == 'G' && geneSeq[l-1] == 'A'))) {
        std::cerr << "Possible error!!! The coding sequence of the gene: " << transcriptName << " does not end with TAG,TAA, or TGA (stop codons). Maybe a SNP? Or an annotation error..." << std::endl;
        *badStartStopCodonFile << transcriptName << "\t" << "stop" << "\t" << geneSeq.substr(l-3,3) << std::endl;
    }
    return divisibleByThree;
}

// To DO:
// 1) Would also be good to distinguish derived vs. ancestral alleles
// 2) Weigh Ns and synsonymous mutations based on their minor and/or derived allele frequencies to get a single score
//    Or perhaps by how often they are homozygous ('fixed' in a species) vs. heterozygous
void getCodingSequenceStatsPhasedSeq(const std::vector<std::string>& allSeqs, const std::string& refSeq, const string& transcriptName, std::vector<string>& statsThisGene, std::ofstream*& prematureStopCodonFile, const std::vector<string>& sampleNames, Annotation& wgAnnotation) {
    int numSegSites = 0;
    int numNonSynAAchanges = 0;  // NUmber of non-synonymous amino-acid changes
    int numSynAAchanges = 0;  // NUmber of non-synonymous amino-acid changes
    double expectedNumNonSynAAchanges = 0;
    double expectedNumSynAAchanges = 0;
    std::vector<double> derivedAlleleFequencies;
    std::vector<double> synonymousMinorAlleleFequencies;
    std::vector<double> nonsynonymousMinorAlleleFequencies;
    int numCopies = (int)allSeqs.size()*2;
    // int numSynonymous = 0;
    // int numNonSynonymous = 0;
    // bool nonsenseMutation = false;
    if (allSeqs[0].length() != refSeq.length()) {
        std::cerr << "Error!!! allSeqs[0].length() != refSeq.length()" << std::endl;
    }
    std::vector<string> altCodons; altCodons.resize(allSeqs.size());
    std::vector<int> IUPACcounts; IUPACcounts.resize(allSeqs.size());
    for (std::vector<string>::size_type i = 0; i != altCodons.size(); i++) {
        altCodons[i] = "";
    }
    
    
    // std::cerr << "Collecting gene sequence statistics...." << std::endl;
    for (string::size_type i = 0; i != refSeq.length(); i++) {
        char refBase = refSeq[i];
        char altBase;
        int countDerived = 0;
        
        for (std::vector<std::string>::size_type j = 0; j != allSeqs.size(); j++) {
            if (allSeqs[j][i] != refBase) {
                altBase = allSeqs[j][i];
                countDerived++;
                altCodons[j] += altBase;
            } else {
                altCodons[j] += refBase;
            }
        }
        // Find the types of mutation we are dealing with
        string refAA; 
        string altAA;
        if ((i+1)%3 == 0) {
            //if (i == 2)
              //  std::cerr << "Got here; refseq length: " << refSeq.length() << " i-2:" << i-2 << std::endl;
            refAA = getAminoAcid(refSeq.substr(i-2,3));
            expectedNumNonSynAAchanges = expectedNumNonSynAAchanges + getExpectedNumberOfNonsynonymousSites(refSeq.substr(i-2,3));
            expectedNumSynAAchanges = expectedNumSynAAchanges + (3 - getExpectedNumberOfNonsynonymousSites(refSeq.substr(i-2,3)));
            //if (i == 2)
              //  std::cerr << "And here" << std::endl;
            string altAA = "";
            //if (i == 2)
            //    std::cerr << "Now going to loop through codons" << std::endl;
            int nonSyn = 0;
            int nonSynV2 = 0;
            int syn = 0;
            int numStops = 0;
            std::vector<string> haveStop;
            for (std::vector<string>::size_type j = 0; j != altCodons.size(); j++) {
                // only consider this individual if it did not have a premature stop codon
                if(std::find(haveStop.begin(), haveStop.end(), sampleNames[j]) != haveStop.end()) {
                    altAA = getAminoAcid(altCodons[j]);
                    if (altAA != refAA && altAA == "Stop") {
                        // We have a premature stop codon
                        haveStop.push_back(sampleNames[j]);
                        numStops++;
                    } else {
                        int thisIndDistance = getCodonDistance(refSeq.substr(i-2,3),altCodons[j]);
                        if (thisIndDistance == 1) {
                            if (altAA != refAA) {
                                nonSyn++;
                                // std::cerr << "altCodons[i]: " << altCodons[j] << "refSeq.substr(i-2,3) " << refSeq.substr(i-2,3) << std::endl;
                                // std::cerr << "allSeqs[0].substr(i-2,3): " << allSeqs[j].substr(i-2,3) << std::endl;
                            }
                            if (!isSingleChangeSynonymous(refSeq.substr(i-2,3),altCodons[j])) {
                                nonSynV2++;
                            }
                            if (nonSyn != nonSynV2) {
                                std::cerr << "refSeq.substr(i-2,3): " << refSeq.substr(i-2,3) << "; altCodons[j]: " << altCodons[j] << std::endl;
                            }
                            assert(nonSyn == nonSynV2);
                            if (altAA == refAA) {
                                syn++;
                            }
                        } else {
                            double NsThisInd = calculateN(refSeq.substr(i-2,3),altCodons[j], thisIndDistance, false);
                            nonSyn = nonSyn + NsThisInd; nonSynV2 = nonSynV2 + NsThisInd;
                            syn = syn + (thisIndDistance - NsThisInd);
                        }
                    }
                }
                altCodons[j] = "";
            }
            if (numStops > 0) {
                std::string stopTranscriptDetails; stopTranscriptDetails.reserve(500);
                stopTranscriptDetails = transcriptName + "\t" + numToString((i+1)/3) + "\t" + numToString(refSeq.length()/3) + "\t" + numToString((double)numStops/(sampleNames.size()*2)) + "\t";
                for (std::vector<string>::size_type i = 0; i != haveStop.size(); i++) {
                    stopTranscriptDetails.append(haveStop[i]);
                    if (i != (haveStop.size()-1))
                        stopTranscriptDetails.append(",");
                }
                *prematureStopCodonFile << stopTranscriptDetails << std::endl;
                //print_vector(haveStop, *prematureStopCodonFile, ',');
                stopsTranscriptRecord.push_back(stopTranscriptDetails);
            }
            double NRaf;
            if (nonSyn > 0) {
                numNonSynAAchanges++;
               // std::cerr << "nonSyn: " << nonSyn << std::endl;
                NRaf = (double)nonSyn/numCopies;
               // std::cerr << "NRaf: " << NRaf << std::endl;
                if (NRaf > 0.5) {
                    nonsynonymousMinorAlleleFequencies.push_back(1-NRaf);
                } else {
                    nonsynonymousMinorAlleleFequencies.push_back(NRaf);
                }
            }
            if (syn > 0) {
                numSynAAchanges++;
                NRaf = (double)syn/numCopies;
                if (NRaf > 0.5) {
                    synonymousMinorAlleleFequencies.push_back(1-NRaf);
                } else {
                    synonymousMinorAlleleFequencies.push_back(NRaf);
                }
            }
        }
        
        if (countDerived > 0) {
            numSegSites++;
            double daf = (double)countDerived/numCopies;
            derivedAlleleFequencies.push_back(daf);
        }
       // if (i%10==0) {
           // std::cerr << "Processed " << i << "bp" << std::endl;
       // }
    }
    
    statsThisGene.push_back(numToString(refSeq.length()));
    statsThisGene.push_back(numToString(numSegSites));
    statsThisGene.push_back(numToString((double)numSegSites/refSeq.length()));
    statsThisGene.push_back(numToString(refSeq.length()/3));
    statsThisGene.push_back(numToString(numSynAAchanges));
    statsThisGene.push_back(numToString((double)numSynAAchanges/(refSeq.length()/3)));
    statsThisGene.push_back(numToString(numNonSynAAchanges));
    statsThisGene.push_back(numToString((double)numNonSynAAchanges/(refSeq.length()/3)));
    statsThisGene.push_back(numToString(vector_average(synonymousMinorAlleleFequencies)));
    statsThisGene.push_back(numToString(vector_average(nonsynonymousMinorAlleleFequencies)));
    
    // std::map<int,int> dafTable = tabulateVectorTemplate(derivedAlleleFequencies);
    // std::cerr << "Stats Done" << std::endl;
}

void getCodingSequencePhasedPnPs(const std::vector<std::string>& allSeqs, const std::string& refSeq, const string& transcriptName, std::vector<string>& statsThisGene, std::ofstream*& prematureStopCodonFile, const std::vector<string>& sampleNames, Annotation& wgAnnotation) {
    double pN = 0; double pS = 0;
    assert(allSeqs[0].length() == refSeq.length());
    std::vector<string> altCodons; altCodons.resize(allSeqs.size());
    
    std::map<std::string, int> haveStop;
    for (std::vector<string>::size_type i = 0; i != altCodons.size(); i++) {
        altCodons[i] = "";
        haveStop[sampleNames[i]] = 0;
    }
    
    std::vector<std::vector<double> > N_d_jk; initialize_matrix_double(N_d_jk, (int)altCodons.size()); 
    std::vector<std::vector<double> > N_jk; initialize_matrix_double(N_jk, (int)altCodons.size());
    std::vector<std::vector<double> > pN_jk; initialize_matrix_double(pN_jk, (int)altCodons.size());
    std::vector<std::vector<double> > S_d_jk; initialize_matrix_double(S_d_jk, (int)altCodons.size());
    std::vector<std::vector<double> > S_jk; initialize_matrix_double(S_jk, (int)altCodons.size());
    std::vector<std::vector<double> > pS_jk; initialize_matrix_double(pS_jk, (int)altCodons.size());
    // std::cerr << "Collecting gene sequence statistics...." << std::endl;
    for (string::size_type i = 0; i != refSeq.length(); i++) {
        for (std::vector<std::string>::size_type j = 0; j != allSeqs.size(); j++) {
            altCodons[j] += allSeqs[j][i];
        }
        // Find the types of mutation we are dealing with
        if ((i+1)%3 == 0) {
            for (std::vector<std::string>::size_type j = 0; j != altCodons.size(); j++) {
                if (getAminoAcid(altCodons[j]) == "Stop") {
                    haveStop[sampleNames[j]] = 1;
                }
            }
            //std::cerr << "Now going to loop through codons: i = " << i << std::endl;
         //     int n_di = 0; int s_di = 0;
         //     double N_i = 0; double S_i = 0;
            for (std::vector<std::string>::size_type j = 0; j != altCodons.size() - 1; j++) {
                if (haveStop[sampleNames[j]] == 1)
                    continue; // only consider this individual if it did not have a premature stop codon
                for (std::vector<std::string>::size_type k = j+1; k != altCodons.size(); k++) {
                    if (haveStop[sampleNames[k]] == 1)
                        continue;
                    int d = getCodonDistance(altCodons[j],altCodons[k]);
                    //std::cerr << "Got codon distance: d = " << d << std::endl;
                    double n_d_ijk = calculateNd(altCodons[j],altCodons[k], d);
                    //std::cerr << "Calculated Nd; n_d_ijk = " << n_d_ijk << std::endl;
                    double s_d_ijk = d - n_d_ijk;
                    //std::cerr << "altCodons[j] = " << altCodons[j] << "; altCodons[k] = " << altCodons[k] << std::endl;
                    //std::cerr << "N_d_jk[j][k] = " << N_d_jk[j][k] << std::endl;
                    //print_matrix(N_d_jk, std::cout);
                    N_d_jk[j][k] = N_d_jk[j][k] + n_d_ijk;
                    S_d_jk[j][k] = S_d_jk[j][k] + s_d_ijk;
                  //  n_di = n_di + n_d_ijk; s_di = s_di + s_d_ijk;
                    double N_ijk = calculateN(altCodons[j],altCodons[k], d, false);
                    // std::cerr << "Calculated N; N_ijk = " << N_ijk << std::endl;
                    double S_ijk = (3 - calculateN(altCodons[j],altCodons[k], d, false));
                    N_jk[j][k] = N_jk[j][k] + N_ijk; S_jk[j][k] = S_jk[j][k] + S_ijk;
                    //N_i = N_i + N_ijk; S_i = S_i + S_ijk;
                }
            }
            for (std::vector<std::string>::size_type j = 0; j != altCodons.size(); j++) {
                altCodons[j] = "";
            }
        }
    }
    
    double sumPn = 0;
    double sumPs = 0;
    for (std::vector<std::string>::size_type j = 0; j != altCodons.size() - 1; j++) {
        for (std::vector<std::string>::size_type k = j+1; k != altCodons.size(); k++) {
            pN_jk[j][k] = N_d_jk[j][k]/N_jk[j][k];
            pS_jk[j][k] = S_d_jk[j][k]/S_jk[j][k];
            sumPn = sumPn + pN_jk[j][k];
            sumPs = sumPs + pS_jk[j][k];
        }
    }
    pN = (2.0/(altCodons.size()*(altCodons.size()-1)))*sumPn;
    pS = (2.0/(altCodons.size()*(altCodons.size()-1)))*sumPs;
    
    statsThisGene.push_back(numToString(pN));
    statsThisGene.push_back(numToString(pS));
    
}






void getCodingSequenceStatsIUPAC(const std::vector<std::string>& allSeqs, const std::string& refSeq, const string& transcriptName, std::vector<string>& statsThisGene, std::ofstream*& prematureStopCodonFile, const std::vector<string>& sampleNames, Annotation& wgAnnotation) {
    int numSegSites = 0;
    int numNonSynAAchanges = 0;  // NUmber of non-synonymous amino-acid changes
    int numSynAAchanges = 0;  // NUmber of non-synonymous amino-acid changes
    double expectedNumNonSynAAchanges = 0;
    double expectedNumSynAAchanges = 0;
    std::vector<double> derivedAlleleFequencies;
    std::vector<double> synonymousMinorAlleleFequencies;
    std::vector<double> nonsynonymousMinorAlleleFequencies;
    int numCopies = (int)allSeqs.size()*2;
    // int numSynonymous = 0;
    // int numNonSynonymous = 0;
    // bool nonsenseMutation = false;
    if (allSeqs[0].length() != refSeq.length()) {
        std::cerr << "Error!!! allSeqs[0].length() != refSeq.length()" << std::endl;
    }
    std::vector<string> altCodons; altCodons.resize(allSeqs.size());
    std::vector<int> IUPACcounts; IUPACcounts.resize(allSeqs.size());
    for (std::vector<string>::size_type i = 0; i != altCodons.size(); i++) {
        altCodons[i] = "";
        IUPACcounts[i] = 0;
    }
    
    
    // std::cerr << "Collecting gene sequence statistics...." << std::endl;
    for (string::size_type i = 0; i != refSeq.length(); i++) {
        char refBase = refSeq[i];
        char altBase;
        int countDerived = 0;
        
        for (std::vector<std::string>::size_type j = 0; j != allSeqs.size(); j++) {
            if (allSeqs[j][i] != refBase) {
                //if (i == 2)
                //    std::cerr << sampleNames[j] << " " << allSeqs[j][i] << std::endl;
                if (isDNAonly(allSeqs[j][i])) {
                    countDerived = countDerived + 2;
                    altBase = allSeqs[j][i];
                } else {
                    countDerived++;
                    IUPACcounts[j]++;
                    altBase = disambiguateIUPAC(refBase, allSeqs[j][i]);
                }
                altCodons[j] += altBase;
            } else {
                altCodons[j] += refBase;
            }
        }
        // Find the types of mutation we are dealing with
        string refAA;
        string altAA;
        if ((i+1)%3 == 0) {
            //if (i == 2)
            //  std::cerr << "Got here; refseq length: " << refSeq.length() << " i-2:" << i-2 << std::endl;
            refAA = getAminoAcid(refSeq.substr(i-2,3));
            expectedNumNonSynAAchanges = expectedNumNonSynAAchanges + getExpectedNumberOfNonsynonymousSites(refSeq.substr(i-2,3));
            expectedNumSynAAchanges = expectedNumSynAAchanges + (3 - getExpectedNumberOfNonsynonymousSites(refSeq.substr(i-2,3)));
            //if (i == 2)
            //  std::cerr << "And here" << std::endl;
            string altAA = "";
            //if (i == 2)
            //    std::cerr << "Now going to loop through codons" << std::endl;
            int nonSyn = 0;
            int nonSynV2 = 0;
            int syn = 0;
            int numStops = 0;
            std::vector<string> haveStop;
            for (std::vector<string>::size_type j = 0; j != altCodons.size(); j++) {
                if (IUPACcounts[j] > 1) // Multiple hets in a codon under IUPAC sequence encoding
                    continue;
                altAA = getAminoAcid(altCodons[j]);
                if (altAA != refAA && altAA == "Stop") {
                    // We have a premature stop codon
                    if (isDNAonlySeq(allSeqs[j].substr(i-2,3))) {
                        haveStop.push_back(sampleNames[j]+"(hom)");
                        numStops = numStops + 2;
                    } else {
                        haveStop.push_back(sampleNames[j]+"(het)");
                        numStops++;
                    }
                } else {
                    int thisIndDistance = getCodonDistance(refSeq.substr(i-2,3),altCodons[j]);
                    if (thisIndDistance == 1) {
                        if (altAA != refAA) {
                            nonSyn++;
                            // std::cerr << "altCodons[i]: " << altCodons[j] << "refSeq.substr(i-2,3) " << refSeq.substr(i-2,3) << std::endl;
                            // std::cerr << "allSeqs[0].substr(i-2,3): " << allSeqs[j].substr(i-2,3) << std::endl;
                        }
                        if (!isSingleChangeSynonymous(refSeq.substr(i-2,3),altCodons[j])) {
                            nonSynV2++;
                        }
                        assert(nonSyn == nonSynV2);
                        if (altAA == refAA) {
                            syn++;
                        }
                    } else {
                        // findCodonPaths(refSeq.substr(i-2,3),altCodons[j]);
                    }
                }
                altCodons[j] = ""; IUPACcounts[j] = 0;
            }
            if (numStops > 0) {
                std::string stopTranscriptDetails; stopTranscriptDetails.reserve(500);
                stopTranscriptDetails = transcriptName + "\t" + numToString((i+1)/3) + "\t" + numToString(refSeq.length()/3) + "\t" + numToString((double)numStops/(sampleNames.size()*2)) + "\t";
                for (std::vector<string>::size_type i = 0; i != haveStop.size(); i++) {
                    stopTranscriptDetails.append(haveStop[i]);
                    if (i != (haveStop.size()-1))
                        stopTranscriptDetails.append(",");
                }
                *prematureStopCodonFile << stopTranscriptDetails << std::endl;
                //print_vector(haveStop, *prematureStopCodonFile, ',');
                stopsTranscriptRecord.push_back(stopTranscriptDetails);
            }
            double NRaf;
            if (nonSyn > 0) {
                numNonSynAAchanges++;
                // std::cerr << "nonSyn: " << nonSyn << std::endl;
                NRaf = (double)nonSyn/numCopies;
                // std::cerr << "NRaf: " << NRaf << std::endl;
                if (NRaf > 0.5) {
                    nonsynonymousMinorAlleleFequencies.push_back(1-NRaf);
                } else {
                    nonsynonymousMinorAlleleFequencies.push_back(NRaf);
                }
            }
            if (syn > 0) {
                numSynAAchanges++;
                NRaf = (double)syn/numCopies;
                if (NRaf > 0.5) {
                    synonymousMinorAlleleFequencies.push_back(1-NRaf);
                } else {
                    synonymousMinorAlleleFequencies.push_back(NRaf);
                }
            }
        }
        
        if (countDerived > 0) {
            numSegSites++;
            double daf = (double)countDerived/numCopies;
            derivedAlleleFequencies.push_back(daf);
        }
        // if (i%10==0) {
        // std::cerr << "Processed " << i << "bp" << std::endl;
        // }
    }
    
    statsThisGene.push_back(numToString(refSeq.length()));
    statsThisGene.push_back(numToString(numSegSites));
    statsThisGene.push_back(numToString((double)numSegSites/refSeq.length()));
    statsThisGene.push_back(numToString(refSeq.length()/3));
    statsThisGene.push_back(numToString(numSynAAchanges));
    statsThisGene.push_back(numToString((double)numSynAAchanges/(refSeq.length()/3)));
    statsThisGene.push_back(numToString(numNonSynAAchanges));
    statsThisGene.push_back(numToString((double)numNonSynAAchanges/(refSeq.length()/3)));
    statsThisGene.push_back(numToString(vector_average(synonymousMinorAlleleFequencies)));
    statsThisGene.push_back(numToString(vector_average(nonsynonymousMinorAlleleFequencies)));
    
    // std::map<int,int> dafTable = tabulateVectorTemplate(derivedAlleleFequencies);
    // std::cerr << "Stats Done" << std::endl;
}





void parseGetCodingSeqOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'n': opt::bIsCoding = false; break;
            case 'p': opt::bUsePartial = true; break;
            case 'H': arg >> opt::hetTreatment; break;
            case 's': arg >> opt::sampleNameFile; break;
            case OPT_ONLY_STATS: opt::bOnlyStats = true; break;
            case 'h':
                std::cout << CODINGSEQ_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (!opt::bIsCoding && opt::bUsePartial) {
        std::cerr << "The -n and -p options are not compatible\n";
        die = true;
    }
    
    if (opt::hetTreatment != 'r' && opt::hetTreatment != 'p' && opt::hetTreatment != 'i') {
        std::cerr << "The -H (--het-treatment) option can only have the values 'r', 'p', or 'i'\n";
        die = true;
    }
    
    if (argc - optind < 3) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 3)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << CODINGSEQ_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::genomeFile = argv[optind++];
    opt::geneFile = argv[optind++];
}

