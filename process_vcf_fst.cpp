//
//  process_vcf_fst.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 03/11/2013.
//  Copyright (c) 2013 Milan Malinsky. All rights reserved.
//

#include "process_vcf_fst.h"
#include "process_vcf_annotation_tools.h"
#include "process_vcf_stats_utils.h"

#define SUBPROGRAM "fst"

static const char *FST_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS]\n"
"Calculate Fst statistic from a vcf, or ms simuation, or summarise eigensoft Fst output\n"
"The SAMPLE_SETS.txt file should have exactly two lines, with each line defining one of the two sample sets\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -n, --run-name                          run-name will be included in the output file name(s)\n\n"
"       To calculate Fst statistics from a vcf:\n"
"       --vcf=VCF_FILE.vcf                      (required)Fst will be calculated over all variants in this file\n"
"       --sets=SAMPLE_SETS.txt                  (required)Define the sets of samples belonging to the same/different populations;\n"
"                                               the file should have exactly two lines, with each line defining one of the two sample sets\n"
"       -s SAMPLES.txt, --samples=SAMPLES.txt   (optional)supply a file of sample identifiers to be used for the VCF file\n"
"                                               (default: sample ids from the vcf file are used)\n"
"       -w SIZE,STEP --window=SIZE,STEP         (optional) sliding window based computation of Fst, Dxy, and  expected heterozygosity\n"
"                                               windows contain SIZE SNPs and move by STEP\n"
"       --ancSets=ANCESTRAL_SAMPLE_SETS.txt     (optional)two sets of samples the form outgroup populations to the two populations for which Fst is calculated\n"
"                                               for particular Fst levels outputs whether the SNPs are segregating in the outgroups\n"
"       --annot=ANNOTATION.gffExtract           (optional)gene annotation in the same format as for the 'getCodingSeq' subprogram\n"
"                                               outputs the location of SNPs with particular Fst levels with respect to exons, introns, UTRs, non-coding regions\n\n"
"       --regions-above=minFst                  (optional, requires -w) outputs the boundaries of regions whose Fst in windows of size set in -w is at least minFst\n"
"                                               the output file has the suffix '_fst_above_minFst.txt'\n"
"       To calculate Fst statistics from ms simulation output:\n"
"       --ms=MS_OUTPUT.txt                      (required)\n"
"       --set1msSimSize=NUM                     (required) set 1 (population 1) was simulated with size NUM in ms\n"
"       --set2msSimSize=NUM                     (required) set 2 (population 2) was simulated with size NUM in ms\n"
"       --set1FstSample=NUM                     (optional) randomly sample NUM individuals from set 1 (population 1) for the Fst calculation\n"
"                                               (default: NUM = set1msSimSize)\n"
"       --set2FstSample=NUM                     (optional) randomly sample NUM individuals from set 2 (population 2) for the Fst calculation\n"
"                                               (default: NUM = set2msSimSize)\n"
"       --msPvals=CUTOFF                        (optinal) ouptput Fisher Exact Test or chi-sq p-values\n\n"
"       To summarise eigensoft Fst output:\n"
"       --eigen=EIGENSOFT_OUTPUT.fst            (required)Fst results from the eigensoft(v.4.2 or 5.0.1) package (smartpca program)\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum { OPT_HELP = 1, OPT_VCF, OPT_SETS, OPT_ANNOT, OPT_MS, OPT_EIGEN, OPT_MS_SET1_SIZE, OPT_MS_SET1_SAMPLE, OPT_MS_SET2_SIZE, OPT_MS_SET2_SAMPLE, OPT_MS_PVALS, OPT_ANC_SETS, OPT_REG_ABOVE  };

static const char* shortopts = "hn:s:w:";

static const struct option longopts[] = {
    { "vcf",   required_argument, NULL, OPT_VCF },
    { "sets",   required_argument, NULL, OPT_SETS },
    { "ancSets",   required_argument, NULL, OPT_ANC_SETS },
    { "window",   required_argument, NULL, 'w' },
    { "regions-above", required_argument, NULL, OPT_REG_ABOVE },
    { "annot",   required_argument, NULL, OPT_ANNOT },
    { "ms",   required_argument, NULL, OPT_MS },
    { "set1msSimSize",   required_argument, NULL, OPT_MS_SET1_SIZE },
    { "set1FstSample",   required_argument, NULL, OPT_MS_SET1_SAMPLE },
    { "set2msSimSize",   required_argument, NULL, OPT_MS_SET2_SIZE },
    { "set2FstSample",   required_argument, NULL, OPT_MS_SET2_SAMPLE },
    { "msPvals", required_argument, NULL, OPT_MS_PVALS },
    { "eigen",   required_argument, NULL, OPT_EIGEN },
    { "samples",   required_argument, NULL, 's' },
    { "run-name",   required_argument, NULL, 'n' },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string sampleSets;
    static string ancSets;
    static string sampleNameFile;
    static int windowSize = 0;
    static int windowStep = 0;
    static double regAbove = 0;
    static string annotFile;
    static string eigensoftFile;
    static string msFile;
    static int msSet1Size = 0;
    static int msSet1FstSample = 0;
    static int msSet2Size = 0;
    static int msSet2FstSample = 0;
    static double msPvalCutoff = 0;
    static string runName = "";
}


SetCounts getVariantCountsForFst(const std::vector<std::string>& fields, const std::vector<size_t>& set1_loci, const std::vector<size_t>& set2_loci) {
    SetCounts thisVariantCounts;
    thisVariantCounts.individualsWithVariant.assign((fields.size()-NUM_NON_GENOTYPE_COLUMNS),0);
    thisVariantCounts.set1individualsWithVariant.assign(set1_loci.size(),0);
    thisVariantCounts.set2individualsWithVariant.assign(set2_loci.size(),0);
    // std::cerr << fields[0] << "\t" << fields[1] << std::endl;
    int set1i = 0; bool inSet1 = false; int set2i = 0; bool inSet2 = false;
    for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
        if (std::find(set1_loci.begin(), set1_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set1_loci.end()) { inSet1 = true; }
        if (std::find(set2_loci.begin(), set2_loci.end(), i-NUM_NON_GENOTYPE_COLUMNS) != set2_loci.end()) { inSet2 = true; }
        
        if (fields[i][0] == '1') {
            thisVariantCounts.overall++;
            if (inSet1) { thisVariantCounts.set1Count++; thisVariantCounts.set1individualsWithVariant[set1i]++; }
            if (inSet2) { thisVariantCounts.set2Count++; thisVariantCounts.set2individualsWithVariant[set2i]++; }
            thisVariantCounts.individualsWithVariant[i- NUM_NON_GENOTYPE_COLUMNS]++;
        }
        if (fields[i][2] == '1') {
            thisVariantCounts.overall++;
            if (inSet1) { thisVariantCounts.set1Count++; thisVariantCounts.set1individualsWithVariant[set1i]++; }
            if (inSet2) { thisVariantCounts.set2Count++; thisVariantCounts.set2individualsWithVariant[set2i]++; }
            thisVariantCounts.individualsWithVariant[i-NUM_NON_GENOTYPE_COLUMNS]++;
        }
        if (inSet1) { set1i++; inSet1 = false; }
        if (inSet2) { set2i++; inSet2 = false; }
    }
    return thisVariantCounts;
}

double calculateDxy(const SetCounts& thisVarCounts, const int n1, const int n2) {
    double Dxy;
    int sumKij = 0;
    for (std::vector<std::string>::size_type i = 0; i != n1/2; i++) {
        for (std::vector<std::string>::size_type j = 0; j != n2/2; j++) {
            if (thisVarCounts.set1individualsWithVariant[i] == 0 && thisVarCounts.set2individualsWithVariant[j] == 0) {
            } else if (thisVarCounts.set1individualsWithVariant[i] == 1 && thisVarCounts.set2individualsWithVariant[j] == 0) {
                sumKij = sumKij + 2;
            } else if (thisVarCounts.set1individualsWithVariant[i] == 0 && thisVarCounts.set2individualsWithVariant[j] == 1) {
                sumKij = sumKij + 2;
            } else if (thisVarCounts.set1individualsWithVariant[i] == 1 && thisVarCounts.set2individualsWithVariant[j] == 1) {
                sumKij = sumKij + 2;
            } else if (thisVarCounts.set1individualsWithVariant[i] == 2 && thisVarCounts.set2individualsWithVariant[j] == 1) {
                sumKij = sumKij + 2;
            } else if (thisVarCounts.set1individualsWithVariant[i] == 1 && thisVarCounts.set2individualsWithVariant[j] == 2) {
                sumKij = sumKij + 2;
            } else if (thisVarCounts.set1individualsWithVariant[i] == 2 && thisVarCounts.set2individualsWithVariant[j] == 2) {
            } else if (thisVarCounts.set1individualsWithVariant[i] == 2 && thisVarCounts.set2individualsWithVariant[j] == 0) {
                sumKij = sumKij + 4;
            } else if (thisVarCounts.set1individualsWithVariant[i] == 0 && thisVarCounts.set2individualsWithVariant[j] == 2) {
                sumKij = sumKij + 4;
            }
        }
    }
    Dxy = (double)sumKij/(n1*n2);
    return Dxy;
}

double calculateFst(const std::vector<double>& fstNumerators, const std::vector<double>& fstDenominators) {
    double numeratorAverage = vector_average(fstNumerators);
    double denominatorAverage = vector_average(fstDenominators);
    double Fst = numeratorAverage/denominatorAverage;
    if (Fst < 0) Fst = 0;
    return Fst;
}

std::vector<double> getSetHeterozygozities(const SetCounts& thisVarCounts, const int n1, const int n2) {
    std::vector<double> heterozygosities;
    double p1 = (double)thisVarCounts.set1Count/n1;
    double p2 = (double)thisVarCounts.set2Count/n2;
    double het1 = calculateExpectedHeterozygositySimple(p1); heterozygosities.push_back(het1);
    double het2 = calculateExpectedHeterozygositySimple(p2); heterozygosities.push_back(het2);
    double hetNei1 = calculateExpectedHeterozygosityNei78(p1,n1); heterozygosities.push_back(hetNei1);
    double hetNei2 = calculateExpectedHeterozygosityNei78(p2,n2); heterozygosities.push_back(hetNei2);
    return heterozygosities;
}


void getFstFromVCF() {
    std::cerr << "Calculating Fst using variants from: " << opt::vcfFile << std::endl;
    std::cerr << "Between the two 'populations' defined in: " << opt::sampleSets << std::endl;
    if (opt::windowSize > 0) {
        std::cerr << "also using a sliding window of size: " << opt::windowSize << " variants and sliding in steps of: " << opt::windowStep << std::endl;
    }
    string fileRoot = stripExtension(opt::sampleSets);
    // Open connection to read from the vcf file
    std::istream* vcfFile = createReader(opt::vcfFile.c_str());
    std::ifstream* setsFile = new std::ifstream(opt::sampleSets.c_str());
    std::ifstream* annotFile;
    std::ofstream* snpCategoryFstFile;
    std::ofstream* regionsAboveFstFile; bool inRegAbove = false;
    std::ofstream* fstDxyFixedWindowFile;
    std::ifstream* ancSetsFile; std::ofstream* ancSetsOutFile;
    std::vector<string> ancSet1; std::vector<string> ancSet2;
    Annotation wgAnnotation;
    if (!opt::annotFile.empty()) {
        annotFile = new std::ifstream(opt::annotFile.c_str());
        Annotation Annot(annotFile, false); // Does not use transcripts annotated as 5' or 3' partial
        wgAnnotation = Annot;
        string snpCategoryFstFileName = fileRoot + "_" + opt::runName + "SNPcategory_fst.txt";
        snpCategoryFstFile = new std::ofstream(snpCategoryFstFileName.c_str());
        *snpCategoryFstFile << "SNPcategory" << "\t" << "thisSNPFst" << "\t" << "thisSNPDxy" << "\t" << "scaffold" << "\t" << "position" << std::endl;
    }
    if (!opt::ancSets.empty()) {
        ancSetsFile = new std::ifstream(opt::ancSets);
        string ancOutFileName = fileRoot + "_" + opt::runName + "ancestralSNPs_fst.txt";
        ancSetsOutFile = new std::ofstream(ancOutFileName);
        *ancSetsOutFile << "scaffold" << "\t" << "position" << "\t" << "AncAllelePopulation" << "\t" << "Fst" << "\t" << "ancSet1_segregating" << "\t" << "ancSet2_segregating" << std::endl;
        string ancSet1String; string ancSet2String;
        getline(*ancSetsFile, ancSet1String);
        getline(*ancSetsFile, ancSet2String);
        ancSet1 = split(ancSet1String, ','); ancSet2 = split(ancSet2String, ',');
        std::sort(ancSet1.begin(),ancSet1.end()); std::sort(ancSet2.begin(),ancSet2.end());
    }
    
    if (opt::regAbove > 0) {
        string regionsAboveFstFileName = fileRoot + "_w_" + numToString(opt::windowSize) + opt::runName + "_fst_above" + numToString(opt::regAbove) + ".txt";
        regionsAboveFstFile = new std::ofstream(regionsAboveFstFileName.c_str());
    }
    
    string FstResultsFileName = fileRoot + "_w_" + numToString(opt::windowSize) + opt::runName + "_fst.txt";
    std::ofstream* pFst = new std::ofstream(FstResultsFileName.c_str());
    fstDxyFixedWindowFile = new std::ofstream(fileRoot + opt::runName + "dXY_fixedWindow.txt");
    string heterozygositySetsFileName = fileRoot + "_w_" + numToString(opt::windowSize) + opt::runName + "_heterozygosity.txt";
    std::ofstream* pHetSets = new std::ofstream(heterozygositySetsFileName.c_str());
    
    
    string set1String; string set2String;
    getline(*setsFile, set1String);
    getline(*setsFile, set2String);
    std::vector<string> set1 = split(set1String, ',');
    std::vector<string> set2 = split(set2String, ',');
    std::sort(set1.begin(),set1.end());
    std::sort(set2.begin(),set2.end());
    
    int numChromosomes;
    int totalVariantNumber = 0;
    int countedVariantNumber = 0;
    string windowMiddleVariant = "first\tWindow";
    string windowStartEnd = "scaffold_0\t0";
    int windowStart = 0; int windowEnd;
    int fixedWindowStart = 0; std::vector<double> fixedWindowDxyVector;
    std::vector<string> sampleNames;
    std::vector<string> fields;
    std::vector<size_t> set1Loci; std::vector<size_t> set2Loci;
    std::vector<size_t> ancSet1Loci; std::vector<size_t> ancSet2Loci;
    short n1; short n2; short n1anc; short n2anc;
    string line;
    std::map<std::string, double> loc_pval;
    std::vector<double> fstNumerators; fstNumerators.reserve(30000000);
    std::vector<double> fstDenominators; fstDenominators.reserve(30000000);
    std::vector<double> DxyVector; DxyVector.reserve(30000000);
    std::vector<std::vector<double> > heterozygositiesVector; heterozygositiesVector.reserve(30000000);
    std::vector<double> set1heterozygositiesSimple; set1heterozygositiesSimple.reserve(30000000);
    std::vector<double> set2heterozygositiesSimple; set2heterozygositiesSimple.reserve(30000000);
    std::vector<double> set1heterozygositiesNei; set1heterozygositiesNei.reserve(30000000);
    std::vector<double> set2heterozygositiesNei; set2heterozygositiesNei.reserve(30000000);
    while (getline(*vcfFile, line)) {
        if (line[0] == '#' && line[1] == '#') {
            
        } else if (line[0] == '#' && line[1] == 'C') {
            std::vector<std::string> fields = split(line, '\t');
            const std::vector<std::string>::size_type numSamples = fields.size() - NUM_NON_GENOTYPE_COLUMNS;
            numChromosomes = (int)numSamples * 2;
            // std::cerr << "Number of chromosomes: " << numChromosomes << std::endl;
            
            if (opt::sampleNameFile.empty()) {
                for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                    sampleNames.push_back(fields[i]);
                }
            } else {
                sampleNames = readSampleNamesFromTextFile(opt::sampleNameFile);
            }
            set1Loci = locateSet(sampleNames, set1);
            set2Loci = locateSet(sampleNames, set2);
            n1 = set1Loci.size()*2; n2 = set2Loci.size()*2;
            std::cerr << "Set1 loci: " << std::endl;
            print_vector_stream(set1Loci, std::cerr);
            std::cerr << "Set2 loci: " << std::endl;
            print_vector_stream(set2Loci, std::cerr);
            
            if (!opt::ancSets.empty()) {
                ancSet1Loci = locateSet(sampleNames, ancSet1);
                ancSet2Loci = locateSet(sampleNames, ancSet2);
                std::cerr << "Ancestral Set1 loci: " << std::endl;
                print_vector_stream(ancSet1Loci, std::cerr);
                std::cerr << "Ancestral Set2 loci: " << std::endl;
                print_vector_stream(ancSet2Loci, std::cerr);
                n1anc = ancSet1Loci.size() * 2; n2anc = ancSet2Loci.size() * 2;
            }
            
            if (opt::windowSize > 0) {
                if (opt::windowSize == opt::windowStep) {
                    *pHetSets << "scaffold" << "\t" << "Start" << "\t" << "End" << "Set1_heterozygosity" << "\t" << "Set2_heterozygosity" << "\t" << "Set1_heterozygosity_Nei" << "\t" << "Set2_heterozygosity_Nei" << std::endl;
                    *pFst << "var_num" << "\t" << "scaffold" << "\t" << "Start" << "\t" << "End" << "\t" << "Fst" << "\t" << "Dxy" << "\t" << "windowSize" << std::endl;
                    if (opt::regAbove > 0) *regionsAboveFstFile << "scaffold" << "\t" << "Start" << "\t" << "End" << std::endl;
                } else {
                    *pHetSets << "Middle_SNP_position" << "\t" << "Set1_heterozygosity" << "\t" << "Set2_heterozygosity" << "\t" << "Set1_heterozygosity_Nei" << "\t" << "Set2_heterozygosity_Nei" << std::endl;
                }
            }
        } else {
            totalVariantNumber++;
            
            std::vector<std::string> fields = split(line, '\t');
            std::vector<std::string> info = split(fields[7], ';');
            if (info[0] != "INDEL") {  // Without indels
                SetCounts counts = getVariantCountsForFst(fields,set1Loci,set2Loci);
                if ((counts.set1Count > 0 || counts.set2Count > 0) && (counts.set1Count < n1 || counts.set2Count < n2)) {
                    countedVariantNumber++;
                    double FstNumerator = calculateFstNumerator(counts, n1, n2); fstNumerators.push_back(FstNumerator);
                    double FstDenominator = calculateFstDenominator(counts, n1, n2); fstDenominators.push_back(FstDenominator);
                    assert(FstDenominator != 0);
                    double thisSNPDxy = calculateDxy(counts, n1, n2); DxyVector.push_back(thisSNPDxy); fixedWindowDxyVector.push_back(thisSNPDxy);
                    std::vector<double> thisSNPhet = getSetHeterozygozities(counts, n1, n2); heterozygositiesVector.push_back(thisSNPhet);
                    set1heterozygositiesSimple.push_back(thisSNPhet[0]); set2heterozygositiesSimple.push_back(thisSNPhet[1]);
                    set1heterozygositiesNei.push_back(thisSNPhet[2]); set2heterozygositiesNei.push_back(thisSNPhet[3]);
                    if (!opt::annotFile.empty()) {
                        string scaffold = fields[0]; string loc = fields[1]; // Scaffold
                        string SNPcategory = wgAnnotation.getCategoryOfSNP(scaffold, loc);
                        double thisSNPFst = FstNumerator/FstDenominator;
                        *snpCategoryFstFile << SNPcategory << "\t" << thisSNPFst << "\t" << thisSNPDxy << "\t" << scaffold << "\t" << loc << std::endl;
                    }
                    if (!opt::ancSets.empty()) {
                        double thisSNPFst = FstNumerator/FstDenominator;
                        if (thisSNPFst < 0) { thisSNPFst = 0; }
                        string AA = split(info[info.size()-1],'=')[1];
                        //std::cerr << "AA=" << " " << AA << std::endl;
                        FourSetCounts c;
                        if (AA == fields[3]) {
                            c = getFourSetVariantCounts(fields,set1Loci,set2Loci,ancSet1Loci,ancSet2Loci,"ref");
                            *ancSetsOutFile << fields[0] << "\t" << fields[1] << "\t" << c.set1daAF-c.set2daAF << "\t" << thisSNPFst << "\t";
                            if (c.set3daAF > 0 & c.set3daAF < 1) { *ancSetsOutFile << "1" << "\t"; } else { *ancSetsOutFile << "0" << "\t"; }
                            if (c.set4daAF > 0 & c.set4daAF < 1) { *ancSetsOutFile << "1" << std::endl; } else { *ancSetsOutFile << "0" << std::endl; }
                        } else if (AA == fields[4]) {
                            c = getFourSetVariantCounts(fields,set1Loci,set2Loci,ancSet1Loci,ancSet2Loci,"alt");
                            *ancSetsOutFile << fields[0] << "\t" << fields[1] << "\t" << c.set1daAF-c.set2daAF << "\t" << thisSNPFst << "\t";
                            if (c.set3daAF > 0 & c.set3daAF < 1) { *ancSetsOutFile << "1" << "\t"; } else { *ancSetsOutFile << "0" << "\t"; }
                            if (c.set4daAF > 0 & c.set4daAF < 1) { *ancSetsOutFile << "1" << std::endl; } else { *ancSetsOutFile << "0" << std::endl; }
                            // std::cerr << "AA=alt" << " " << c.set1daAF << " " << c.set2daAF << std::endl;
                        } else {
                            c = getFourSetVariantCounts(fields,set1Loci,set2Loci,ancSet1Loci,ancSet2Loci,"N");
                            *ancSetsOutFile << fields[0] << "\t" << fields[1] << "\t" << "-888" << "\t" << thisSNPFst << "\t";
                            if (c.set3AltAF > 0 & c.set3AltAF < 1) { *ancSetsOutFile << "1" << "\t"; } else { *ancSetsOutFile << "0" << "\t"; }
                            if (c.set4AltAF > 0 & c.set4AltAF < 1) { *ancSetsOutFile << "1" << std::endl; } else { *ancSetsOutFile << "0" << std::endl; }
                        }
                        
                        
                    }
                    
                    std::vector<string> s = split(windowStartEnd, '\t');
                    if (s[0] == fields[0]) {
                        if (atoi(fields[1].c_str()) > (fixedWindowStart+10000)) {
                            double thisFixedWindowDxy = vector_average_withRegion(fixedWindowDxyVector, 10000);
                            *fstDxyFixedWindowFile << fields[0] << "\t" << fixedWindowStart << "\t" << fixedWindowStart+10000 << "\t" << thisFixedWindowDxy << std::endl;
                            fixedWindowDxyVector.clear(); fixedWindowStart= fixedWindowStart+10000;
                        }
                    } else {
                        fixedWindowStart = 0;
                    }
                    
                    
                    
                    if (opt::windowSize == 1) {
                        double Fst = FstNumerator/FstDenominator;
                        if (Fst < 0) Fst = 0;
                        *pFst << countedVariantNumber << "\t" << fields[0] + "\t" + fields[1] << "\t" << Fst << "\t" << thisSNPDxy << std::endl;
                        
                    } else if ((opt::windowSize > 0) && (countedVariantNumber % opt::windowStep == 0) && countedVariantNumber >= opt::windowSize) {
                        std::vector<double> windowFstNumerators(fstNumerators.end()-opt::windowSize, fstNumerators.end());
                        std::vector<double> windowFstDenominators(fstDenominators.end()-opt::windowSize, fstDenominators.end());
                        double windowFst = calculateFst(windowFstNumerators, windowFstDenominators); if (windowFst < 0) windowFst = 0;
                        std::vector<double> windowDxyVec(DxyVector.end()-opt::windowSize, DxyVector.end());
                        double windowDxy = vector_average(windowDxyVec);
                        if (opt::windowSize == opt::windowStep) {
                            std::vector<string> s = split(windowStartEnd, '\t');
                            if (s[0] == fields[0]) {
                                windowStartEnd = windowStartEnd + "\t" + fields[1];
                                windowEnd = atoi(fields[1].c_str());
                                double windowDxyIncNonSeg = vector_average_withRegion(windowDxyVec, windowEnd-windowStart);
                                *pFst << countedVariantNumber-opt::windowSize+1 << "\t" << windowStartEnd << "\t" << windowFst << "\t" << windowDxy << "\t" << windowDxyIncNonSeg << "\t" << windowFstDenominators.size() << std::endl;
                                if (opt::regAbove > 0) {
                                    if (windowFst >= opt::regAbove && !inRegAbove) {
                                        inRegAbove = true;
                                        *regionsAboveFstFile << s[0] << "\t" << s[1] << "\t";
                                    } else if (windowFst < opt::regAbove && inRegAbove) {
                                        inRegAbove = false;
                                        *regionsAboveFstFile << s[1] << std::endl;
                                    }
                                }
                            }
                        } else {
                            *pFst << countedVariantNumber-opt::windowSize+1 << "\t" << windowMiddleVariant << "\t" << windowFst << "\t" << windowDxy << "\t" << windowFstDenominators.size() << std::endl;
                        }
                        // Now calculate and output expected heterozygosities for this window
                        std::vector<double> windowHetS1Vec(set1heterozygositiesSimple.end()-opt::windowSize, set1heterozygositiesSimple.end());
                        double windowHetS1 = vector_average(windowHetS1Vec);
                        std::vector<double> windowHetS2Vec(set2heterozygositiesSimple.end()-opt::windowSize, set2heterozygositiesSimple.end());
                        double windowHetS2 = vector_average(windowHetS2Vec);
                        std::vector<double> windowHetNei1Vec(set1heterozygositiesNei.end()-opt::windowSize, set1heterozygositiesNei.end());
                        double windowHetNei1 = vector_average(windowHetNei1Vec);
                        std::vector<double> windowHetNei2Vec(set2heterozygositiesNei.end()-opt::windowSize, set2heterozygositiesNei.end());
                        double windowHetNei2 = vector_average(windowHetNei2Vec);
                        if (opt::windowSize == opt::windowStep) {
                            std::vector<string> s = split(windowStartEnd, '\t');
                            if (s[0] == fields[0]) {
                                *pHetSets << windowStartEnd << "\t" << windowHetS1 << "\t" << windowHetS2 << "\t" << windowHetNei1 << "\t" << windowHetNei2 << std::endl;
                                windowStartEnd = fields[0] + "\t" + fields[1];
                                windowStart = atoi(fields[1].c_str());
                            } else {
                                windowStartEnd = fields[0] + "\t0";
                                windowStart = 0;
                            }
                        } else {
                            *pHetSets << windowMiddleVariant << "\t" << windowHetS1 << "\t" << windowHetS2 << "\t" << windowHetNei1 << "\t" << windowHetNei2 << std::endl;
                            windowMiddleVariant = fields[0] + "\t" + fields[1];     // works only if STEP is half SIZE for the window
                        }
                    }
                }
            }
            if (totalVariantNumber % 100000 == 0) {
                double Fst = calculateFst(fstNumerators, fstDenominators);
                std::cerr << totalVariantNumber << " variants processed... Fst: " << Fst << std::endl;
            }
                
        }
    }
    double Fst = calculateFst(fstNumerators, fstDenominators);
    double overallHetS1 = vector_average(set1heterozygositiesSimple);
    double overallHetS2 = vector_average(set2heterozygositiesSimple);
    double overallHetNei1 = vector_average(set1heterozygositiesNei);
    double overallHetNei2 = vector_average(set2heterozygositiesNei);
    
    std::cerr << "Fst: " << Fst << std::endl;
    std::cerr << "Heterozygosities: " << "\tS1:" << overallHetS1 << "\tS2:" << overallHetS2 << "\tNei1:" << overallHetNei1 << "\tNei2" << overallHetNei2 << std::endl;
    *pHetSets << "#Heterozygosities: " << "\tS1:" << overallHetS1 << "\tS2:" << overallHetS2 << "\tNei1:" << overallHetNei1 << "\tNei2" << overallHetNei2 << std::endl;
}

void getFstFromMs() {
    std::cerr << "Calculating Fst using variants from: " << opt::msFile << std::endl;
    std::cerr << "and outputting chi-sq test p-vals < " << opt::msPvalCutoff << std::endl;
    
    std::ifstream* msFile = new std::ifstream(opt::msFile.c_str());
    string fileRoot = stripExtension(opt::msFile);
    string PvalFileName = fileRoot + "_" + opt::runName + "_pvals.txt";
    std::ofstream* pValFile;
    if (opt::msPvalCutoff > 0) {
        pValFile = new std::ofstream(PvalFileName.c_str());
        *pValFile << "Fisher p-val" << "\t" << "chi-sq pval" << "\t" << "set1Alt" << "\t" << "set1Ref" << "\t" << "set2Alt" << "\t" << "set2Ref" << "\t" << "Fst" << std::endl;
    }
    
    
    std::vector<int> set1_loci;
    std::vector<int> set2_loci;
    srand((int)time(NULL));
    if (opt::msSet1FstSample == 0) {
        opt::msSet1FstSample = opt::msSet1Size;
        for (int i = 0; i != opt::msSet1FstSample; i++) {
            set1_loci.push_back(i);
        }
    } else { // Randomly sample individuals from population 1 for Fst calculation
        for (int i = 0; i != opt::msSet1FstSample; i++) {
            int rand_sample = (rand()%opt::msSet1Size);
            while (std::find(set1_loci.begin(),set1_loci.end(),rand_sample) != set1_loci.end()) {
                rand_sample = (rand()%opt::msSet1Size);
            }
            set1_loci.push_back(rand_sample);
        }
    }
    // Do the same for set2
    if (opt::msSet2FstSample == 0) {
        opt::msSet2FstSample = opt::msSet2Size;
        for (int i = 0; i != opt::msSet2FstSample; i++) {
            set2_loci.push_back(i+opt::msSet1Size);
        }
    } else { // Randomly sample individuals from population 2 for Fst calculation
        for (int i = 0; i != opt::msSet2FstSample; i++) {
            int rand_sample = (rand()%opt::msSet2Size)+opt::msSet1Size;
            while (std::find(set2_loci.begin(),set2_loci.end(),rand_sample) != set2_loci.end()) {
                rand_sample = (rand()%opt::msSet2Size)+opt::msSet1Size;
            }
            set2_loci.push_back(rand_sample);
        }
    }

    std::cerr << "Selected population 1 individuals: "; print_vector_stream(set1_loci, std::cerr);
    std::cerr << "Selected population 2 individuals: "; print_vector_stream(set2_loci, std::cerr);
    
    if (opt::msSet1Size != opt::msSet1FstSample || opt::msSet2Size != opt::msSet2FstSample) {
        std::cerr << "Warning: the Fst column is going to contain '-1' values where the site is not a segregating site in the sampled individuals for Fst calcultation" << std::endl;
    }
    
    std::vector<double> fstNumerators; fstNumerators.reserve(500000000);
    std::vector<double> fstDenominators; fstDenominators.reserve(500000000);
    
    
    string line;
    int numFixedSites = 0;
    int numNearlyFixedSites = 0;
    std::vector<double> nullForChisq;
    std::vector<int> moreSet1;
    std::vector<int> lessSet1;
    std::vector<int> moreSet2;
    std::vector<int> lessSet2;
    while (getline(*msFile, line)) {
        SetCounts counts;
        double thisFst = -1;
        for (std::vector<int>::iterator it = set1_loci.begin(); it != set1_loci.end(); it++) {
            // std::cerr << line[*it] << std::endl;
            if (line[*it] == '1') {
                counts.set1Count++;
            }
        }
        for (std::vector<int>::iterator it = set2_loci.begin(); it != set2_loci.end(); it++) {
            if (line[*it] == '1') {
                counts.set2Count++;
            }
        }
        
        if (counts.set1Count > 0 || counts.set2Count > 0) {
            double FstNum = calculateFstNumerator(counts, opt::msSet1FstSample, opt::msSet2FstSample);
            double FstDenom = calculateFstDenominator(counts, opt::msSet1FstSample, opt::msSet2FstSample);
            thisFst = FstNum/FstDenom; if (thisFst < 0) thisFst = 0;
            fstNumerators.push_back(FstNum);
            fstDenominators.push_back(FstDenom);
            
        }
        
        if ((counts.set1Count == 0 && counts.set2Count == opt::msSet2FstSample) || (counts.set1Count == opt::msSet1FstSample && counts.set2Count == 0)) {
            numFixedSites++;
        }
        
        if ((counts.set1Count == 1 && counts.set2Count == opt::msSet2FstSample) || (counts.set1Count == 0 && counts.set2Count == opt::msSet2FstSample-1) ||
            (counts.set1Count == opt::msSet1FstSample-1 && counts.set2Count == 0) || (counts.set1Count == opt::msSet1FstSample && counts.set2Count == 1)) {
            numNearlyFixedSites++;
        }
        
        
        int set1WithoutVariant = opt::msSet1FstSample-counts.set1Count;
        int set2WithoutVariant = opt::msSet2FstSample-counts.set2Count;
        
        if (counts.set1Count >= set1WithoutVariant) {
            moreSet1.push_back(counts.set1Count);
            lessSet1.push_back(set1WithoutVariant);
            moreSet2.push_back(counts.set2Count);
            lessSet2.push_back(set2WithoutVariant);
        } else {
            moreSet1.push_back(set1WithoutVariant);
            lessSet1.push_back(counts.set1Count);
            moreSet2.push_back(set2WithoutVariant);
            lessSet2.push_back(counts.set2Count);
        }
        
        if ((counts.set1Count != 0 || counts.set2Count != 0) && (set1WithoutVariant != 0 || set2WithoutVariant != 0)) {
            if (opt::msSet1FstSample + opt::msSet2FstSample <= 60) {
                counts.fisher_pval = fisher_exact(counts.set1Count,set1WithoutVariant , counts.set2Count, set2WithoutVariant);
   //             std::cerr << "Fisher: " << counts.fisher_pval << std::endl;
                counts.chi_sq_pval = pearson_chi_sq_indep(counts.set1Count,set1WithoutVariant , counts.set2Count, set2WithoutVariant);
            } else {
                counts.chi_sq_pval = pearson_chi_sq_indep(counts.set1Count,set1WithoutVariant , counts.set2Count, set2WithoutVariant);
            }
        }
        
        if (counts.fisher_pval < opt::msPvalCutoff || counts.chi_sq_pval < opt::msPvalCutoff) {
            *pValFile << counts.fisher_pval << "\t" << counts.chi_sq_pval << "\t" << counts.set1Count << "\t" << set1WithoutVariant << "\t" << counts.set2Count << "\t" << set2WithoutVariant << "\t" << thisFst << std::endl;
        }
    }
    
    double Fst = calculateFst(fstNumerators, fstDenominators);
    std::cerr << "Fst: " << Fst << std::endl;
    std::cerr << "Fixed sites: " << numFixedSites << std::endl;
    std::cerr << "Tier2 sites: " << numNearlyFixedSites << std::endl;
    std::cerr << "Null ChiSq 1:" << vector_average(moreSet1)/opt::msSet1FstSample << "\t" << vector_average(lessSet1)/opt::msSet1FstSample << std::endl;
    std::cerr << "Null ChiSq 2:" << vector_average(moreSet2)/opt::msSet2FstSample << "\t" << vector_average(lessSet2)/opt::msSet2FstSample << std::endl;
}

void summariseEigensoft() {
    std::ifstream* eigenFile = new std::ifstream(opt::eigensoftFile.c_str());
    string fileRoot = stripExtension(opt::eigensoftFile);
    string FstResultsFileName = fileRoot + "_" + opt::runName + "_fst_matrix.forR";
    std::ofstream* pFst = new std::ofstream(FstResultsFileName.c_str());
    std::vector<std::vector<std::string> > fst_matrix;
    
    string line;
    getline(*eigenFile, line); // Get the first description line
    short type;
    if (line == "## Fst statistics between populations:         fst       std error") {
        type = 1;
    } else {
        type = 2;
    }
    if (type == 1) {
        getline(*eigenFile, line);
        std::vector<std::string> fields = split(line, '\t');
        std::vector<std::string> this_indiv_fst;
        std::vector<std::string> all_indiv;
        string this_indiv = fields[0];
        this_indiv_fst.push_back(fields[2]);
        while (getline(*eigenFile, line)) {
            fields = split(line, '\t');
            if (this_indiv == fields[0]) {
                this_indiv_fst.push_back(fields[2]);
            } else {
                fst_matrix.push_back(this_indiv_fst);
                all_indiv.push_back(this_indiv);
                this_indiv = fields[0];
                this_indiv_fst.clear();
                this_indiv_fst.push_back(fields[2]);
            }
        }
        all_indiv.push_back(this_indiv);
        fst_matrix.push_back(this_indiv_fst);
        this_indiv_fst.clear(); this_indiv_fst.push_back("0"); all_indiv.push_back(fields[1]);
        fst_matrix.push_back(this_indiv_fst);
        
        for (std::vector<std::vector<std::string> >::size_type i = 0; i != fst_matrix.size(); i++) {
            std::reverse(fst_matrix[i].begin(), fst_matrix[i].end());
            fst_matrix[i].insert(fst_matrix[i].end(), "0");
            while (fst_matrix[i].size() != fst_matrix[0].size()) {
                fst_matrix[i].insert(fst_matrix[i].end(), "0");
            }
        }
        std::reverse(fst_matrix.begin(), fst_matrix.end());
        std::reverse(all_indiv.begin(), all_indiv.end());
        
        print_vector(all_indiv, *pFst);
        print_matrix(fst_matrix, *pFst);
    } else if (type == 2) {
        std::cerr << "type2" << std::endl;
        std::vector<std::string> fields = split(line, '\t');
        std::vector<std::string> all_indiv(fields.begin()+1,fields.end());
        getline(*eigenFile, line); getline(*eigenFile, line);
        std::vector<std::string> this_indiv_fst;
        while (getline(*eigenFile, line)) {
            fields = split(line, '\t');
            std::copy(fields.begin()+1,fields.end(),std::back_inserter(this_indiv_fst));
            for (std::vector<std::string>::size_type i = 0; i != this_indiv_fst.size(); i++) {
                double fst = convertToDouble(this_indiv_fst[i]) / 1000;
                this_indiv_fst[i] = numToString(fst);
            }
            fst_matrix.push_back(this_indiv_fst);
            this_indiv_fst.clear();
        }
        print_vector(all_indiv, *pFst);
        print_matrix(fst_matrix, *pFst);
    }
}

int fstMain(int argc, char** argv) {
    parseFstOptions(argc, argv);
    
    if (!opt::vcfFile.empty())
        getFstFromVCF();
    if (!opt::eigensoftFile.empty())
        summariseEigensoft();
    if (!opt::msFile.empty())
        getFstFromMs();
    
    return 0;
}

void parseFstOptions(int argc, char** argv) {
    bool die = false;
    std::vector<string> windowSizeStep;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case OPT_VCF: arg >> opt::vcfFile; break;
            case OPT_SETS: arg >> opt::sampleSets; break;
            case OPT_ANC_SETS: arg >> opt::ancSets; break;
            case 'w':
                windowSizeStep = split(arg.str(), ',');
                opt::windowSize = atoi(windowSizeStep[0].c_str());
                opt::windowStep = atoi(windowSizeStep[1].c_str());
                break;
            case OPT_REG_ABOVE: arg >> opt::regAbove; break;
            case OPT_ANNOT: arg >> opt::annotFile; break;
            case OPT_MS: arg >> opt::msFile; break;
            case OPT_EIGEN: arg >> opt::eigensoftFile; break;
            case OPT_MS_SET1_SIZE: arg >> opt::msSet1Size; break;
            case OPT_MS_SET1_SAMPLE: arg >> opt::msSet1FstSample; break;
            case OPT_MS_SET2_SIZE: arg >> opt::msSet2Size; break;
            case OPT_MS_SET2_SAMPLE: arg >> opt::msSet2FstSample; break;
            case OPT_MS_PVALS: arg >> opt::msPvalCutoff; break;
            case 's': arg >> opt::sampleNameFile; break;
            case 'n': arg >> opt::runName; break;
            case '?': die = true; break;
            case 'h': std::cout << FST_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    if (argc - optind > 0)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if ((!opt::vcfFile.empty() && opt::sampleSets.empty()) || (opt::vcfFile.empty() && !opt::sampleSets.empty())) {
        std::cerr << "To calculate Fst statistics from a vcf, you need to use both the --vcf and the --sets options\n";
        die = true;
    }
    
    if ((!opt::msFile.empty() && (opt::msSet1Size == 0 || opt::msSet2Size == 0)) || (opt::msFile.empty() && (opt::msSet1Size != 0 || opt::msSet2Size != 0))) {
        std::cerr << "To calculate Fst statistics from ms output, you need to use the --ms and the --set1msSimSize and --set2msSimSize  options together\n";
        die = true;
    }
    
    if (opt::msPvalCutoff > 1) {
        std::cerr << "The pvalue cutoff cannot be more than one...\n";
        die = true;
    }
    
    if (opt::regAbove > 0 && opt::windowSize == 0) {
        std::cerr << "The window size (-w option) needs to be set for --regions-above option to work (you can set -w 1,1 for per variant Fst)\n";
        die = true;
    }
    
    
    if (die) {
        std::cout << "\n" << FST_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    
    
}