//
//  process_vcf_coding_sequences.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 16/09/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include "process_vcf_coding_sequences.h"
#include <ctime>
 
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
"       --no-stats                              do not calculate statistics for coding sequences\n"
"                                               only output the sequences themselves\n"
"       --output-nondiv-3=prefix                output genes whose coding sequence length is not divisible by three\n"
"                                               the specified prefix will be added to output file names\n"
"       -H,   --het-treatment <r|p|b|i>         r: assign het bases randomly (default); p: use the phase information in a VCF outputting haplotype 1 for each individual; b: use both haplotypes as phased; i: use IUPAC codes\n"
"       -s SAMPLES.txt, --samples=SAMPLES.txt   supply a file of sample identifiers to be used for the output\n"
"                                               (default: sample ids from the vcf file are used)\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hpws:cH:";
static std::vector<string> stopsTranscriptRecord;

enum { OPT_ONLY_STATS, OPT_NONDIV_THREE, OPT_USE_PHASE };


static const struct option longopts[] = {
    { "non-coding",   required_argument, NULL, 'n' },
    { "samples",   required_argument, NULL, 's' },
    { "partial",   no_argument, NULL, 'p' },
    { "het-treatment",   required_argument, NULL, 'H' },
    { "output-nondiv-3", required_argument, NULL, OPT_NONDIV_THREE },
    { "no-stats",   no_argument, NULL, OPT_ONLY_STATS },
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
    static bool bNoStats = false;
    static char hetTreatment = 'r';
    static string nonDivPrefix = "";
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
    std::vector<string> scaffoldStringsH2;
    std::vector<std::vector<string> > statsAllGenes;
    string statsFileName = geneFileRoot + "_stats.txt";
    std::ofstream* statsFile = new std::ofstream(statsFileName.c_str());
    string stopsFileName = geneFileRoot + "_prematureStops.txt";
    *statsFile << "transcript" << "\t" << "length_in_nucleotides" << "\t" << "segregating_sites(ss)" << "\t" << "ss_proportion" << "\t" << "length_in_AA" << "\t" << "num_of_AA_with_synonymous_changes(synAAs)" << "\t" << "synAAs_proportion" << "\t" << "non_synonymous_AA_substitutions(nsAAs)" << "\t" << "nsAAs_proportion" << "\t" << "synonymousMAFaverage" << "\t" << "nonsynonymousMAFaverage" << "\t" << "pN" << "\t" << "pS" << std::endl;
    if (opt::hetTreatment == 'b') {
        std::cout << "transcript" << "\t" << "pN" << "\t" << "pS" << std::endl;
    } else {
        std::cout << "transcript" << "\t" << "length_in_nucleotides" << "\t" << "segregating_sites(ss)" << "\t" << "ss_proportion" << "\t" << "length_in_AA" << "\t" << "num_of_AA_with_synonymous_changes(synAAs)" << "\t" << "synAAs_proportion" << "\t" << "non_synonymous_AA_substitutions(nsAAs)" << "\t" << "nsAAs_proportion" << "\t" << "synonymousMAFaverage" << "\t" << "nonsynonymousMAFaverage" << "\t" << "pN" << "\t" << "pS" << std::endl;
    }
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
                    if (opt::hetTreatment == 'b') {
                        for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                            scaffoldStringsH2[i].append(currentScaffoldReference.substr(inStrPos, string::npos));
                        }
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
                        std::vector<string> allSeqs; std::vector<string> allSeqsH2;
                        string geneName = annotLineVec[4];
                        std::vector<string> statsThisGene; statsThisGene.push_back(geneName);
                        string refSeq = getReferenceForThisRegion(annotation[k], annotLineVec[3], currentScaffoldReference);
                        if (opt::bIsCoding)
                            geneLengthDivisibleByThree = codingSequenceErrorChecks(refSeq, geneName, annotation, (int)k, badStartStopCodonFile);
                        if (geneLengthDivisibleByThree)
                            geneOutFiles = new std::ofstream(geneName.c_str());
                        else if (opt::nonDivPrefix != "")
                            geneOutFiles = new std::ofstream((opt::nonDivPrefix+'_'+geneName).c_str());
                        if (geneLengthDivisibleByThree || opt::nonDivPrefix != "") {
                            for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                                string geneSeq = getIndividualSequenceForThisRegion(annotation[k], annotLineVec[3], scaffoldStrings[i]);
                                *geneOutFiles << ">" << sampleNames[i] << std::endl;
                                *geneOutFiles << geneSeq << std::endl;
                                allSeqs.push_back(geneSeq);
                                if (opt::hetTreatment == 'b') {
                                    geneSeq = getIndividualSequenceForThisRegion(annotation[k], annotLineVec[3], scaffoldStringsH2[i]);
                                    *geneOutFiles << ">" << sampleNames[i] << "_H2" << std::endl;
                                    *geneOutFiles << geneSeq << std::endl;
                                    allSeqsH2.push_back(geneSeq);
                                }
                            } geneOutFiles->close();
                            std::cerr << "Got all sequences for: " << annotLineVec[4] << std::endl;
                        }
                        // Get statistics for the sequences
                        if (!opt::bNoStats) {
                            if (opt::bIsCoding && geneLengthDivisibleByThree) {
                                assert(allSeqs[0].length() == refSeq.length());
                                assert(allSeqsH2[0].length() == refSeq.length());
                                if (opt::hetTreatment == 'p' || opt::hetTreatment == 'r') {
                                    getStatsHaploidSeq(allSeqs, statsThisGene);
                                    print_vector_stream(statsThisGene, std::cout);
                                } else if (opt::hetTreatment == 'b') {
                                    std::vector<std::vector<double> > combinedVectorForPCA;
                                    pNsets* sets = nullptr;
                                    getStatsBothPhasedHaps(allSeqs, allSeqsH2, statsThisGene, combinedVectorForPCA, sets);
                                    print_vector_stream(statsThisGene, std::cout);
                                } else {
                                    getStatsIUPAC(allSeqs, refSeq, annotLineVec[4], statsThisGene,stopsFile, sampleNames);
                                }
                                statsAllGenes.push_back(statsThisGene);
                                std::cerr << "Got statistics for: " << annotLineVec[4] << std::endl;
                            }
                        }
                    }
                  //  for (std::vector<std::string>::size_type i = 0; i != numSamples; i++) {
                  //      scaffoldStrings[i] = "";
                  //      if (opt::hetTreatment == 'b') scaffoldStringsH2[i] = "";
                  //  }
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
                
                for (std::vector<string>::size_type i = 0; i != numSamples; i++) {
                    scaffoldStrings[i] = "";
                    scaffoldStrings[i].reserve(currentScaffoldReference.length());
                }
                if (opt::hetTreatment == 'b') {
                    scaffoldStringsH2.resize(numSamples);
                    for (std::vector<string>::size_type i = 0; i != numSamples; i++) {
                        scaffoldStringsH2[i] = "";
                        scaffoldStringsH2[i].reserve(currentScaffoldReference.length());
                    }
                }
                
            }
            //if (wgAnnotation.annotationMap.count(fields[0]) == 1) {
                std::string ref = fields[3];  std::string alt = fields[4];
                if (ref.length() == 1 && alt.length() == 1) { // Only use Biallelic SNPs (or multiallelics split across multiple lines)
                    for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                        //std::cerr << "Going through genotypes1:" << i << std::endl;
                        //std::cerr << scaffoldStrings.size() << " " << inStrPos << " " << fields[1] << " " << currentScaffoldReference.size() << std::endl;
                        int sampleNum = (int)i - NUM_NON_GENOTYPE_COLUMNS;
                        scaffoldStrings[sampleNum].append(currentScaffoldReference.substr(inStrPos, (atoi(fields[1].c_str()) - 1)-inStrPos));
                        if (opt::hetTreatment == 'b')
                            scaffoldStringsH2[sampleNum].append(currentScaffoldReference.substr(inStrPos, (atoi(fields[1].c_str()) - 1)-inStrPos));
                        std::vector<string> genotypeFields = split(fields[i], ':');
                        std::vector<char> genotype;
                        genotype.push_back(genotypeFields[0][0]); genotype.push_back(genotypeFields[0][2]);
                        if (opt::hetTreatment == 'b') {
                            // Somewhat hacky:
                            // 1) add to the first haplotype
                            appendGenotypeBaseToString(scaffoldStrings[sampleNum], fields[3], fields[4], genotype, 'p');
                            // 2) add to the second haplotype
                            appendGenotypeBaseToString(scaffoldStringsH2[sampleNum], fields[3], fields[4], genotype, 'b');
                        } else {
                            appendGenotypeBaseToString(scaffoldStrings[sampleNum], fields[3], fields[4], genotype, opt::hetTreatment);
                        }
                    }
                    inStrPos = atoi(fields[1].c_str());
                    
    #ifdef DEBUG
                    if (currentScaffoldReference[inStrPos-1] != fields[3][0]) {
                        std::cerr << "Error!!! Sequence: " << currentScaffoldReference[inStrPos-1] << " vcf-ref: " << fields[3][0] << std::endl;
                    }
    #endif
                }
            //}
            if (processedVariantCounter % 10000 == 0) {
                std::cerr << processedVariantCounter << " variants processed..." << std::endl;
                //std::cout << "scaffoldStrings[0]:\t" << scaffoldStrings[0] << std::endl;
                //std::cout << "scaffoldStringsH2[0]:\t" << scaffoldStringsH2[0] << std::endl;
                //std::cout << "ScaffoldReference:\t" << currentScaffoldReference.substr(0,inStrPos) << std::endl;
            }
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
 
void getStatsHaploidSeq(const std::vector<std::string>& allSeqs, std::vector<string>& statsThisGene, double tStVratio) {
    std::string::size_type geneLengthNt = allSeqs[0].length();
    int numSamples = (int)allSeqs.size();
    double pN = 0; double pS = 0;
    std::vector<string> altCodons; altCodons.resize(numSamples);
    std::map<std::vector<string>::size_type, int> haveStop;
    
    for (std::vector<string>::size_type i = 0; i != numSamples; i++) {
        assert(allSeqs[i].length() == geneLengthNt);
        altCodons[i] = ""; haveStop[i] = 0; 
    }
    
    CDSComparisonMatrices* pairwiseMatrices = new CDSComparisonMatrices(numSamples);
    for (string::size_type i = 0; i != geneLengthNt; i++) {
        for (std::vector<std::string>::size_type j = 0; j != numSamples; j++) {
            altCodons[j] += allSeqs[j][i];
        }
        // Find the types of mutation we are dealing with
        if ((i+1)%3 == 0) {
            for (std::vector<std::string>::size_type j = 0; j != numSamples; j++) {
                if (getAminoAcid(altCodons[j]) == "Stop") haveStop[j] = 1;
            }
           // std::cerr << "Now going to loop through codons: i = " << i << "altCodons.size():" << altCodons.size() << std::endl;
            addAllPairwiseN_S_Nd_Sd_DifferentIndividuals(altCodons,haveStop, *pairwiseMatrices);
            //std::cerr << "Added pairwise among H1:" << std::endl;
            for (std::vector<std::string>::size_type j = 0; j != numSamples; j++) {
             //   std::cerr << "j = " << j << "; altCodons[j]: " << altCodons[j] << std::endl;
                altCodons[j] = "";
               // std::cerr << "j = " << j << "; haveStop[j]: " << haveStop[j] << std::endl;
            }
        }
    }
     
    double sumPn = 0; double sumPs = 0;
    std::vector<double> pNjackknifeVector;
    std::vector<double> pSjackknifeVector;
    std::vector<double> pN_pSjackknifeVector;
    // Add the within H1 and within H2 comparisons
    for (std::vector<std::string>::size_type j = 0; j != numSamples - 1; j++) {
        //std::cerr << "j = " << j << "; sumPn: " << sumPn << "; sumPs:" << sumPs << std::endl;
        for (std::vector<std::string>::size_type k = j+1; k != numSamples; k++) {
            //double pN_jk = pairwiseMatrices.H1p->N_d_jk[j][k]/pairwiseMatrices.H1p->N_jk[j][k];
            double pN_jk = pairwiseMatrices->N_d_jk[j][k]/((2*tStVratio*pairwiseMatrices->tS_N_jk[j][k])+pairwiseMatrices->tV_N_jk[j][k]);
            if (std::isnan(pN_jk) || pairwiseMatrices->tS_N_jk[j][k] == 0) {
            std::cerr << "j = " << j << "; k = " << k << std::endl;
            std::cerr << "N_d_jk[j][k] = " << pairwiseMatrices->N_d_jk[j][k] << "; tS_N_jk[j][k] = " << pairwiseMatrices->tS_N_jk[j][k] << std::endl;
            std::cerr << "pN_jk = " << pN_jk << "; sumPn = " << sumPn << std::endl;
                print_matrix(pairwiseMatrices->N_d_jk, std::cerr);
                print_matrix(pairwiseMatrices->S_d_jk, std::cerr);
            }
            sumPn = sumPn + pN_jk;
            double pS_jk = pairwiseMatrices->S_d_jk[j][k]/((2*tStVratio*pairwiseMatrices->tS_S_jk[j][k])+pairwiseMatrices->tV_S_jk[j][k]);
            sumPs = sumPs + pS_jk;
    
            // make sure each sample is used only once in the jackknife calculation
            // to keep the observations independent
            if ((j % 2 == 0) && (k == j+1)) {
                pNjackknifeVector.push_back(pN_jk); pSjackknifeVector.push_back(pS_jk);
                pN_pSjackknifeVector.push_back(pN_jk-pS_jk);
            }
        }
    }
    pN = (2.0/(numSamples*(numSamples-1)))*sumPn;
    pS = (2.0/(numSamples*(numSamples-1)))*sumPs;
    
    statsThisGene.push_back(numToString(geneLengthNt));
    statsThisGene.push_back(numToString(pN));
    statsThisGene.push_back(numToString(pS));
    if (pNjackknifeVector.size() > 10) {
        double pNstdErr = jackknive_std_err(pNjackknifeVector);
        double pSstdErr = jackknive_std_err(pSjackknifeVector);
        double pNpSstdErr = jackknive_std_err(pN_pSjackknifeVector);
        statsThisGene.push_back(numToString(pNstdErr));
        statsThisGene.push_back(numToString(pSstdErr));
        statsThisGene.push_back(numToString(pNpSstdErr));
    } else {
        statsThisGene.push_back("NA");
        statsThisGene.push_back("NA");
        statsThisGene.push_back("NA");
    }
}


// To DO:
// 1) Get some sort of inbreeding coefficient (how often are they homozygous ('fixed' in a species) vs. heterozygous)
// 2) Would also be good to distinguish derived vs. ancestral alleles?
void getStatsBothPhasedHaps(const std::vector<std::string>& allSeqs, const std::vector<std::string>& allSeqsH2, std::vector<string>& statsThisGene, std::vector<std::vector<double> >& combinedVectorForPCA, pNsets* sets, double tStVratio, bool nonCodingNull) {
    //std::cerr << "Collecting gene sequence statistics...." << std::endl;
    //clock_t begin = clock();
    double pN = 0; double pS = 0;
    double hetN = 0; double hetS = 0;
    assert(allSeqs.size() == allSeqsH2.size());
    std::string::size_type geneLengthNt = allSeqs[0].length();
    int numSamples = (int)allSeqs.size();
    std::vector<string> altCodons; altCodons.resize(numSamples);
    std::vector<string> altCodonsH2; altCodonsH2.resize(numSamples);
    initialize_matrix_double(combinedVectorForPCA, numSamples);
    
    std::map<std::vector<string>::size_type, int> haveStop;
    std::map<std::vector<string>::size_type, int> haveStopH2;
    for (std::vector<string>::size_type i = 0; i != numSamples; i++) {
        assert(allSeqs[i].length() == geneLengthNt);
        assert(allSeqsH2[i].length() == geneLengthNt);
        altCodons[i] = ""; altCodonsH2[i] = "";
        haveStop[i] = 0; haveStopH2[i] = 0;
    }
    
    CDSH1H2ComparisonMatrices pairwiseMatrices(numSamples);
    for (string::size_type i = 0; i != geneLengthNt; i++) {
        for (std::vector<std::string>::size_type j = 0; j != allSeqs.size(); j++) {
            altCodons[j] += allSeqs[j][i];
            altCodonsH2[j] += allSeqsH2[j][i];
        }
        // Find the types of mutation we are dealing with
        if ((i+1)%3 == 0) {
            if (!nonCodingNull) {
                // If this is an alignment for a coding region:
                // the remainder of the sequence after any stop will be excluded from the calculations for any sequence
                for (std::vector<std::string>::size_type j = 0; j != numSamples; j++) {
                    if (getAminoAcid(altCodons[j]) == "Stop") {
                        //std::cerr << "Stop in : i = " << i << " ; j = " << j << std::endl;
                        haveStop[j] = 1;
                    }
                    if (getAminoAcid(altCodonsH2[j]) == "Stop") haveStopH2[j] = 1;
                    
                }
            }
            //std::cerr << "Now going to loop through codons: i = " << i << std::endl;
            addAllPairwiseN_S_Nd_Sd_DifferentIndividuals(altCodons,haveStop, *pairwiseMatrices.H1p);
            //std::cerr << "Added pairwise among H1:" << std::endl;
            addAllPairwiseN_S_Nd_Sd_DifferentIndividuals(altCodonsH2,haveStopH2, *pairwiseMatrices.H2p);
            //std::cerr << "Added pairwise among H2: " << std::endl;
            addN_S_Nd_Sd_DifferentIndividualsH1againstH2(altCodons, altCodonsH2, haveStop, haveStopH2, *pairwiseMatrices.H1H2p);
            //std::cerr << "Added pairwise beween H1 and H2:" << std::endl;
            addN_S_Nd_Sd_SameIndividualsH1againstH2(altCodons, altCodonsH2, haveStop, haveStopH2, *pairwiseMatrices.H1H2p);
            //std::cerr << "Added heterozygosity beween H1 and H2:" << std::endl;
            
            for (std::vector<std::string>::size_type j = 0; j != numSamples; j++) {
                altCodons[j] = ""; altCodonsH2[j] = "";
            }
        }
    }
     
    double sumPn = 0; double sumPs = 0;
    // double sumS = 0; double sumN = 0;
    std::vector<double> pNjackknifeVector;
    std::vector<double> pSjackknifeVector;
    std::vector<double> pN_pSjackknifeVector;
    std::vector<double> pN_pSjackknifeVector_allComparisons;
    // Add the within H1 and within H2 comparisons
    for (std::vector<std::string>::size_type j = 0; j != numSamples - 1; j++) {
        //std::cerr << "j = " << j << "; sumPn: " << sumPn << "; sumPs:" << sumPs << std::endl;
        for (std::vector<std::string>::size_type k = j+1; k != numSamples; k++) {
            //double pN_jk = pairwiseMatrices.H1p->N_d_jk[j][k]/pairwiseMatrices.H1p->N_jk[j][k];
            double pN_jk = pairwiseMatrices.H1p->N_d_jk[j][k]/((2*tStVratio*pairwiseMatrices.H1p->tS_N_jk[j][k])+pairwiseMatrices.H1p->tV_N_jk[j][k]); sumPn = sumPn + pN_jk;
            //sumN = sumN + (2*tStVratio*pairwiseMatrices.H1p->tS_N_jk[j][k]) + pairwiseMatrices.H1p->tV_N_jk[j][k];
            combinedVectorForPCA[j][k] = combinedVectorForPCA[j][k] + pN_jk;
            //if (isnan(pN_jk) || N_jk[j][k] == 0) {
                //std::cerr << "j = " << j << "; k = " << k << std::endl;
                //std::cerr << "N_d_jk[j][k] = " << N_d_jk[j][k] << "; N_jk[j][k] = " << N_jk[j][k] << std::endl;
                //std::cerr << "pN_jk = " << N_d_jk[j][k] << "; sumPn = " << sumPn << std::endl;
            //}
            double pS_jk = pairwiseMatrices.H1p->S_d_jk[j][k]/((2*tStVratio*pairwiseMatrices.H1p->tS_S_jk[j][k])+pairwiseMatrices.H1p->tV_S_jk[j][k]); sumPs = sumPs + pS_jk;
            //sumS = sumS + pairwiseMatrices.H1p->tS_S_jk[j][k]+pairwiseMatrices.H1p->tV_S_jk[j][k];
            double H2pN_jk = pairwiseMatrices.H2p->N_d_jk[j][k]/((2*tStVratio*pairwiseMatrices.H2p->tS_N_jk[j][k])+pairwiseMatrices.H2p->tV_N_jk[j][k]); sumPn = sumPn + H2pN_jk;
            //sumN = sumN + pairwiseMatrices.H2p->tS_N_jk[j][k]+pairwiseMatrices.H2p->tV_N_jk[j][k];
            combinedVectorForPCA[j][k] = combinedVectorForPCA[j][k] + H2pN_jk;
            //std::cerr << "H2N_d_jk[j][k] = " << H2N_d_jk[j][k] << "; H2N_jk[j][k] = " << H2N_jk[j][k] << std::endl;
            //std::cerr << "H2pN_jk = " << N_d_jk[j][k] << "; sumPn = " << sumPn << std::endl;
            double H2pS_jk = pairwiseMatrices.H2p->S_d_jk[j][k]/((2*tStVratio*pairwiseMatrices.H2p->tS_S_jk[j][k])+pairwiseMatrices.H2p->tV_S_jk[j][k]); sumPs = sumPs + H2pS_jk;
            pN_pSjackknifeVector_allComparisons.push_back(pN_jk-pS_jk); pN_pSjackknifeVector_allComparisons.push_back(H2pN_jk-H2pS_jk);
            if ((j % 2 == 0) && (k == j+1)) {
                pNjackknifeVector.push_back(pN_jk); pNjackknifeVector.push_back(H2pN_jk);
                pSjackknifeVector.push_back(pS_jk); pSjackknifeVector.push_back(H2pS_jk);
                pN_pSjackknifeVector.push_back(pN_jk-pS_jk); pN_pSjackknifeVector.push_back(H2pN_jk-H2pS_jk);
            }
            
            if (sets->initialised == true) {
                if ((sets->set1Loci.count(j) == 1 && sets->set1Loci.count(k) == 1) || (sets->set1Loci.count(k) == 1 && sets->set1Loci.count(j) == 1)) {
                    sets->withinSet1andSet2pN = sets->withinSet1andSet2pN + pN_jk + H2pN_jk;
                    sets->withinSet1pN = sets->withinSet1pN + pN_jk + H2pN_jk;
                } else if ((sets->set2Loci.count(j) == 1 && sets->set2Loci.count(k) == 1) || (sets->set2Loci.count(k) == 1 && sets->set2Loci.count(j) == 1)) {
                    sets->withinSet1andSet2pN = sets->withinSet1andSet2pN + pN_jk + H2pN_jk;
                    sets->withinSet2pN = sets->withinSet2pN + pN_jk + H2pN_jk;
                } else if ((sets->set1Loci.count(j) == 1 && sets->set2Loci.count(k) == 1)) {
                    sets->set1vsSet2pN = sets->set1vsSet2pN + pN_jk + H2pN_jk;
                    sets->withinSet1andSet2pN = sets->withinSet1andSet2pN + pN_jk + H2pN_jk;
                } else if (sets->set1Loci.count(k) == 1 && sets->set2Loci.count(j) == 1) {
                    sets->set1vsSet2pN = sets->set1vsSet2pN + pN_jk + H2pN_jk;
                    sets->withinSet1andSet2pN = sets->withinSet1andSet2pN + pN_jk + H2pN_jk;
                } else if ((sets->set1Loci.count(j) == 1 || sets->set2Loci.count(j) == 1) && sets->set3Loci.count(k) == 1) {
                    sets->sets1and2vsSet3pN = sets->sets1and2vsSet3pN + pN_jk + H2pN_jk;
                } else if ((sets->set1Loci.count(k) == 1 || sets->set2Loci.count(k) == 1) && sets->set3Loci.count(j) == 1) {
                    sets->sets1and2vsSet3pN = sets->sets1and2vsSet3pN + pN_jk + H2pN_jk;
                    //std::cerr << "j = " << j << "; (set 3): " << allSeqs[j] << std::endl;
                    //std::cerr << "k = " << k << "; (set 1 or 2): " << allSeqs[k] << std::endl;
                    //std::cerr << std::endl;
                } else if ((sets->set3Loci.count(j) == 1 && sets->set3Loci.count(k) == 1) || (sets->set3Loci.count(k) == 1 && sets->set3Loci.count(j) == 1)) {
                    sets->withinSet3pN = sets->withinSet3pN + pN_jk + H2pN_jk;
                }
            }
        }
    }
    
    double sumHetPn = 0; double sumHetPs = 0;
    // Add the between H1 and H2 comparisons
    for (std::vector<std::string>::size_type j = 0; j != numSamples; j++) {
        //std::cerr << "j = " << j << "; sumPn: " << sumPn << "; sumPs:" << sumPs << std::endl;
        for (std::vector<std::string>::size_type k = 0; k != numSamples; k++) {
            if (j != k) {
                //std::cerr << ((2*tStVratio*pairwiseMatrices.H1H2p->tS_N_jk[j][k])+pairwiseMatrices.H1H2p->tV_N_jk[j][k]) << std::endl;
                double H1H2pN_jk = pairwiseMatrices.H1H2p->N_d_jk[j][k]/((2*tStVratio*pairwiseMatrices.H1H2p->tS_N_jk[j][k])+pairwiseMatrices.H1H2p->tV_N_jk[j][k]);
                double H1H2pS_jk = pairwiseMatrices.H1H2p->S_d_jk[j][k]/((2*tStVratio*pairwiseMatrices.H1H2p->tS_S_jk[j][k])+pairwiseMatrices.H1H2p->tV_S_jk[j][k]);
                if (j < k)
                    combinedVectorForPCA[j][k] = combinedVectorForPCA[j][k] + H1H2pN_jk;
                else
                    combinedVectorForPCA[k][j] = combinedVectorForPCA[k][j] + H1H2pN_jk;
                sumPn = sumPn + H1H2pN_jk;
                sumPs = sumPs + H1H2pS_jk;
                pN_pSjackknifeVector_allComparisons.push_back(H1H2pN_jk-H1H2pS_jk);
                
                if (sets->initialised == true) {
                    if ((sets->set1Loci.count(j) == 1 && sets->set1Loci.count(k) == 1) || (sets->set1Loci.count(k) == 1 && sets->set1Loci.count(j) == 1)) {
                        sets->withinSet1andSet2pN = sets->withinSet1andSet2pN + H1H2pN_jk;
                        sets->withinSet1pN = sets->withinSet1pN + H1H2pN_jk;
                    } else if ((sets->set2Loci.count(j) == 1 && sets->set2Loci.count(k) == 1) || (sets->set2Loci.count(k) == 1 && sets->set2Loci.count(j) == 1)) {
                        sets->withinSet1andSet2pN = sets->withinSet1andSet2pN + H1H2pN_jk;
                        sets->withinSet2pN = sets->withinSet2pN + H1H2pN_jk;
                    } else if ((sets->set1Loci.count(j) == 1 && sets->set2Loci.count(k) == 1) || (sets->set1Loci.count(k) == 1 && sets->set2Loci.count(j) == 1)) {
                        sets->set1vsSet2pN = sets->set1vsSet2pN + H1H2pN_jk;
                        sets->withinSet1andSet2pN = sets->withinSet1andSet2pN + H1H2pN_jk;
                    } if ((sets->set1Loci.count(j) == 1 || sets->set2Loci.count(j) == 1) && sets->set3Loci.count(k) == 1) {
                        sets->sets1and2vsSet3pN = sets->sets1and2vsSet3pN + H1H2pN_jk;
                    } else if ((sets->set1Loci.count(k) == 1 || sets->set2Loci.count(k) == 1) && sets->set3Loci.count(j) == 1) {
                        sets->sets1and2vsSet3pN = sets->sets1and2vsSet3pN + H1H2pN_jk;
                    } else if ((sets->set3Loci.count(j) == 1 && sets->set3Loci.count(k) == 1) || (sets->set3Loci.count(k) == 1 && sets->set3Loci.count(j) == 1)) {
                        sets->withinSet3pN = sets->withinSet3pN + H1H2pN_jk;
                    }
                }
            } else {
                //std::cerr << ((2*tStVratio*pairwiseMatrices.H1H2p->tS_N_jk[j][j])+pairwiseMatrices.H1H2p->tV_N_jk[j][j]) << std::endl;
                double H1H2pN_jj = pairwiseMatrices.H1H2p->N_d_jk[j][j]/((2*tStVratio*pairwiseMatrices.H1H2p->tS_N_jk[j][j])+pairwiseMatrices.H1H2p->tV_N_jk[j][j]);
                double H1H2pS_jj = pairwiseMatrices.H1H2p->S_d_jk[j][j]/((2*tStVratio*pairwiseMatrices.H1H2p->tS_S_jk[j][j])+pairwiseMatrices.H1H2p->tV_S_jk[j][j]);
                sumHetPn = sumHetPn + H1H2pN_jj;
                sumHetPs = sumHetPs + H1H2pS_jj;
            }
        }
    }

    double pNstdErr = jackknive_std_err(pNjackknifeVector);
    double pSstdErr = jackknive_std_err(pSjackknifeVector);
    double pNpSstdErr = jackknive_std_err(pN_pSjackknifeVector);
    double pNpSstdErrAllComparisons = jackknive_std_err(pN_pSjackknifeVector_allComparisons);
    pN = sumPn/(2*(numSamples*(numSamples-1)));
    pS = sumPs/(2*(numSamples*(numSamples-1)));
    hetN = sumHetPn/(double)numSamples;
    hetS = sumHetPs/(double)numSamples;
    
    statsThisGene.push_back(numToString(geneLengthNt));
    statsThisGene.push_back(numToString(pN));
    statsThisGene.push_back(numToString(pS));
    statsThisGene.push_back(numToString(hetN));
    statsThisGene.push_back(numToString(hetS));
    statsThisGene.push_back(numToString(pNstdErr));
    statsThisGene.push_back(numToString(pSstdErr));
    statsThisGene.push_back(numToString(pNpSstdErr));
    statsThisGene.push_back(numToString(pNpSstdErrAllComparisons));
   // clock_t end = clock();
   // double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //std::cerr << "Time per gene: " << elapsed_secs << std::endl;
}



void getStatsIUPAC(const std::vector<std::string>& allSeqs, const std::string& refSeq, const string& transcriptName, std::vector<string>& statsThisGene, std::ofstream*& prematureStopCodonFile, const std::vector<string>& sampleNames) {
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
                        if (s_mapCodonPairToSynonymous.count(refSeq.substr(i-2,3)+altCodons[j]) == 0) {
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
            case OPT_ONLY_STATS: opt::bNoStats = true; break;
            case OPT_NONDIV_THREE: arg >> opt::nonDivPrefix; break;
            case 'h':
                std::cout << CODINGSEQ_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (!opt::bIsCoding && opt::bUsePartial) {
        std::cerr << "The -n and -p options are not compatible\n";
        die = true;
    }
    
    if (opt::hetTreatment != 'r' && opt::hetTreatment != 'p' && opt::hetTreatment != 'b' && opt::hetTreatment != 'i') {
        std::cerr << "The -H (--het-treatment) option can only have the values 'r', 'p', 'b', or 'i'\n";
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

