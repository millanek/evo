//
//  process_vcf_filter.cpp
//  vcf_process
//
//  Created by Milan Malinsky on 03/06/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include <iostream>
#include "process_vcf_utils.h"
#include "process_vcf_filter.h"

// Filtering constants
// static const int MIN_OVERALL_VARIANT_PHRED_QUAL=30;   //

#define SUBPROGRAM "filter"

static const char *FILTER_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] VCF_FILE\n"
"Custom file for VCF filtering, output to STD_OUT\n"
"\n"
"       --help                           display this help and exit\n"
"\n"
"       The filtereing parameters that can be changed are:\n"
"       -d, --overall-max-depth (DEFAULT +Inf)  Maximum read depth allowed at the putative variant site - filters out variants due to\n"
"                                               collapsed repeats in the reference (reads from multiple sites are all mapped to one)\n"
"       -m, --overall-min-depth (DEFAULT 0)     Minimum read depth allowed at the putative variant site - filters out strange regions where very few reads align\n"
"       -c, --min-copies=MIN (DEFAULT 1)        The variant needs to be present in at least MIN copies\n"
"                                               (i.e. setting --min_copies==1 gives all segregating sites including singletons)\n"
"                                               ( and setting --min_copies==0 also includes sites where all individuals are different from the reference)\n"
"       -s, --min-depth-per-sample=MIN (DEFAULT 0)        Minimum read depth at the variant position for any sample\n"
"                                               (i.e. all samples need to have at least MIN reads at the variant position)\n"
"       --minOverallQuality=MIN  (DEFAULT 30)   The minimum accepted phred scaled [10*log_10*p(variant)] variant quality (30 corresponds to 0.001 probability of error)\n"
"       --minF=MIN  (DEFAULT -1)                Minimum inbreeding coefficient; Negative values mean an excess of hets compared with Hardy Weinberg expectation (2*p*q)\n"
"       --minMQSB=MIN  (DEFAULT 0)              Minimum p-value for Mann-Whitney U test of Mapping Quality vs Strand Bias (common in samtools calling output)\n"
"       --maxFS=MAX  (DEFAULT +Inf)             Maximum phred-scaled p-value using Fisher's exact test to detect strand bias (common in GATK calling output)\n"
"       --allow-missing                         Allow missing genotypes (./.) where read depth < min-depth-per-sample\n"
"       --keep-triallelic                       Do not filter out sites that are not biallelic\n"
"       --max-het-individuals=N (DEFAULT +Inf)  Filter out sites where more than N individuals are heterozygous\n"
"\n"
"\n"
"       You can also just output statistics on the first 1,000,000 variants:\n"
"       --stats                                 Get stats to help with selecting parameters for filtering\n"
"       --subsample=p                           Don't use every variant, but a random subsample, selecting variants with probability p\n"
"\n\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";


enum { OPT_HELP = 1, OPT_TRIALLELIC, OPT_ALLOW_MISSING_GENOTYPES, OPT_MAX_NUM_HET, OPT_STATS, OPT_MQSB, OPT_FS, OPT_OVERALL_Q, OPT_MIN_F, OPT_SUBS };

static const char* shortopts = "d:c:s:m:";

static const struct option longopts[] = {
    { "overall-max-depth",       required_argument, NULL, 'd' },
    { "overall-min-depth",       required_argument, NULL, 'm' },
    { "min-copies",       required_argument, NULL, 'c' },
    { "min-depth-per-sample",       required_argument, NULL, 's' },
    { "help",   no_argument, NULL, OPT_HELP },
    { "max-het-individuals", required_argument, NULL, OPT_MAX_NUM_HET },
    { "keep-triallelic",   no_argument, NULL, OPT_TRIALLELIC },
    { "allow-missing",   no_argument, NULL, OPT_ALLOW_MISSING_GENOTYPES },
    { "stats",   no_argument, NULL, OPT_STATS },
    { "subsample",   required_argument, NULL, OPT_SUBS },
    { "minMQSB",   required_argument, NULL, OPT_MQSB },
    { "maxFS",   required_argument, NULL, OPT_FS },
    { "minF",   required_argument, NULL, OPT_MIN_F },
    { "minOverallQuality",   required_argument, NULL, OPT_OVERALL_Q },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static int min_copies=1;
    static int max_overall_depth = std::numeric_limits<int>::max();
    static int min_overall_depth = 0;
    static int max_het_indiv = std::numeric_limits<int>::max();
    static int min_depth_in_any_individual = 0;
    static int min_overall_quality = 30;
    static double min_F = -1.0;
    static double min_MQSB = 0;
    static double max_FS = std::numeric_limits<double>::max();
    static bool bBiallelicFilter = true;
    static bool bAllowMissingGenotpyes = false;
    static bool bStats = false;
    static double subsample = 1.0;
    static string vcfFile;
}



int filterMain(int argc, char** argv) {
    parseFilterOptions(argc, argv);
    //std::ios_base::openmode mode = std::ios_base::in;
    string fileName = opt::vcfFile;
    string fileRoot = stripExtension(fileName);
    std::vector<int> depthsHetFailed;
    std::vector<int> depthsHetPassed;
    std::vector<int> numVariantsPerHetCount;
    std::vector<std::vector<int> > num_indiv_het_vs_depth;
    
    std::ofstream* statsFFile;
    std::ofstream* statsVarDepthFile;
    std::ofstream* statsStrandBiasFile;
    std::ofstream* statsVariantQualityFile;
    std::ofstream* statsSegregationMetricFile;
    
    // Collect some numbers about why variants were filtered out:
    int numTotalFilteredOut = 0;
    int numMultiallelic = 0; int numInvariant = 0;
    // SB - strand bias; OD - overall depth; OQ - overall quality; F - inbreeding coefficient;
    int numSB = 0; int numOD = 0; int numOQ = 0; int numF = 0; int numOQ_OD = 0; int numOQ_SB = 0; int numOQ_F = 0; int numSB_F = 0; int numSB_OD = 0; int numF_OD = 0;
    int numOQ_OD_F = 0; int numOQ_OD_SB = 0; int numOQ_SB_F = 0; int numSB_F_OD = 0; int numOQ_OD_F_SB = 0;
    
    
    
    if (opt::bStats) {
        statsFFile = new std::ofstream(fileRoot + ".inbreeding");
        statsVarDepthFile = new std::ofstream(fileRoot + ".varDepth");
        statsStrandBiasFile = new std::ofstream(fileRoot + ".strandBias");
        statsVariantQualityFile = new std::ofstream(fileRoot + ".varQual");
        statsSegregationMetricFile = new std::ofstream(fileRoot + ".SGB");
        std::cerr << "Getting statistics for variants from: " << fileName << " to help with setting filtering criteria" << std::endl;
    } else {
        std::cerr << "Filtering variants from: " << fileName << std::endl;
        std::cerr << "Minimum variant quality score (overall) set to: " << opt::min_overall_quality << std::endl;
        std::cerr << "Maximum read depth (overall) set to: " << opt::max_overall_depth << std::endl;
        std::cerr << "Minimum read depth (overall) set to: " << opt::min_overall_depth << std::endl;
        std::cerr << "Minimum copies for a variant (e.g. 1 for allowing singletons): " << opt::min_copies << std::endl;
        std::cerr << "Minimum read depth at the variant position: " << opt::min_depth_in_any_individual << std::endl;
        if (opt::bAllowMissingGenotpyes)
            std::cerr << "Genotypes for individuals with lower read depth at the variant position will be set to missing (./.)" << std::endl;
        if (opt::max_het_indiv < std::numeric_limits<int>::max())
            std::cerr << "Filter out sites where more than " << opt::max_het_indiv << " individuals are heterozygous"  << std::endl;
    }
    // Open connection to read from the vcf file
    // std::ifstream* inFile = new std::ifstream(fileName.c_str(), mode);
    std::istream* inFile = createReader(fileName.c_str());
    
    
    bool gotChromosomeNumber = false;
    int numChromosomes;
    string line;
    int totalVariantNumber = 0;
    while (getline(*inFile, line)) {
        if (line[0] == '#') {
            if (!opt::bStats) std::cout << line << std::endl;
        } else {
            double r = ((double) rand() / (RAND_MAX));
            if (r > opt::subsample) {
                continue;
            }
            totalVariantNumber++;
            if (totalVariantNumber % 100000 == 0) std::cerr << "Filtered " << totalVariantNumber << " variants" << std::endl;
            FilterResult result;
            std::vector<std::string> fields = split(line, '\t');
            if (!gotChromosomeNumber) {
                const std::vector<std::string>::size_type numSamples = fields.size() - NUM_NON_GENOTYPE_COLUMNS;
                numChromosomes = (int)numSamples * 2;
                std::cerr << "Number of chromosomes: " << numChromosomes << std::endl;
                numVariantsPerHetCount.assign(numSamples + 1, 0);
                gotChromosomeNumber = true;
            }
            result.overallQuality = stringToDouble(fields[5]);
            
            result.counts = getThisVariantCounts(fields);
            
            if (opt::bStats) {
                
                *statsFFile << result.counts.inbreedingCoefficient << std::endl;
                *statsVarDepthFile << result.counts.overallDepth << std::endl;
                if (!result.counts.FSpval.empty()) {
                    *statsStrandBiasFile << result.counts.FSpval << std::endl;
                } else if (!result.counts.MQSBpval.empty()) {
                    *statsStrandBiasFile << result.counts.MQSBpval << std::endl;
                }
                if (result.counts.SGB < std::numeric_limits<double>::max()) {
                    *statsSegregationMetricFile << result.counts.SGB << std::endl;
                }
                *statsVariantQualityFile << result.overallQuality << std::endl;
                if (totalVariantNumber >= 1000000)
                    exit(EXIT_SUCCESS);
                else
                    continue;
            }
            
            
            // Start filtering
            // 1) Throw away stuff that is not bi-allelic
            if (opt::bBiallelicFilter) {
                result.biallelicPassed = testBiallelic(fields[4]);
                if (!result.biallelicPassed) { result.biallelicPassed = false; numTotalFilteredOut++; numMultiallelic++; continue; }
            } else
                result.biallelicPassed = true;
            
            // 2) Filtering on overall quality
            if (result.overallQuality < opt::min_overall_quality)
                result.overallQualityPassed = false;
            
            // 3) Filtering on overall depth
            if (result.counts.overallDepth >= opt::min_overall_depth && result.counts.overallDepth <= opt::max_overall_depth) {
                result.overallDepthPassed = true;
            } else
                result.overallDepthPassed = false;
            
            // 4) Filtering on strand bias
            if (!result.counts.FSpval.empty()) {
                if (stringToDouble(result.counts.FSpval) > opt::max_FS)
                    result.strandBiasPassed = false;
            }
            
            if (!result.counts.MQSBpval.empty()) {
                if (stringToDouble(result.counts.MQSBpval) < opt::min_MQSB)
                    result.strandBiasPassed = false;
            }
            
            // 5) Filter on inbreeding coefficinet
            if (result.counts.inbreedingCoefficient < opt::min_F) {
                result.inbreedingCoeffPassed = false;
            }
            
            // 6) Filtering on the number of copies/allele frequency
            if (result.counts.overall <= (numChromosomes - opt::min_copies) && result.counts.overall >= opt::min_copies) {
                
            } else {
                numInvariant++; continue;
            }
            
            if (result.inbreedingCoeffPassed && result.overallDepthPassed && result.overallQualityPassed && result.strandBiasPassed) {
                std::cout << line << std::endl;
            } else {
                numTotalFilteredOut++;
                if (!result.inbreedingCoeffPassed && result.overallDepthPassed && result.overallQualityPassed && result.strandBiasPassed) {
                    numF++;
                } else if (result.inbreedingCoeffPassed && !result.overallDepthPassed && result.overallQualityPassed && result.strandBiasPassed) {
                    numOD++;
                } else if (result.inbreedingCoeffPassed && result.overallDepthPassed && !result.overallQualityPassed && result.strandBiasPassed) {
                    numOQ++;
                } else if (result.inbreedingCoeffPassed && result.overallDepthPassed && result.overallQualityPassed && !result.strandBiasPassed) {
                    numSB++;
                } else if (!result.inbreedingCoeffPassed && result.overallDepthPassed && !result.overallQualityPassed && result.strandBiasPassed) {
                    numOQ_F++;
                } else if (result.inbreedingCoeffPassed && !result.overallDepthPassed && !result.overallQualityPassed && result.strandBiasPassed) {
                    numOQ_OD++;
                } else if (result.inbreedingCoeffPassed && result.overallDepthPassed && !result.overallQualityPassed && !result.strandBiasPassed) {
                    numOQ_SB++;
                } else if (!result.inbreedingCoeffPassed && result.overallDepthPassed && result.overallQualityPassed && !result.strandBiasPassed) {
                    numSB_F++;
                } else if (result.inbreedingCoeffPassed && !result.overallDepthPassed && result.overallQualityPassed && !result.strandBiasPassed) {
                    numSB_OD++;
                } else if (!result.inbreedingCoeffPassed && !result.overallDepthPassed && result.overallQualityPassed && result.strandBiasPassed) {
                    numF_OD++;
                } else if (!result.inbreedingCoeffPassed && !result.overallDepthPassed && !result.overallQualityPassed && result.strandBiasPassed) {
                    numOQ_OD_F++;
                } else if (result.inbreedingCoeffPassed && !result.overallDepthPassed && !result.overallQualityPassed && !result.strandBiasPassed) {
                    numOQ_OD_SB++;
                } else if (!result.inbreedingCoeffPassed && result.overallDepthPassed && !result.overallQualityPassed && !result.strandBiasPassed) {
                    numOQ_SB_F++;
                } else if (!result.inbreedingCoeffPassed && !result.overallDepthPassed && result.overallQualityPassed && !result.strandBiasPassed) {
                    numSB_F_OD++;
                } else if (!result.inbreedingCoeffPassed && !result.overallDepthPassed && !result.overallQualityPassed && !result.strandBiasPassed) {
                    numOQ_OD_F_SB++;
                }
            }
            
            // 7) Filtering on per-individual criteria
            
            /*if (result.counts.minimumDepthInAnIndividual >= opt::min_depth_in_any_individual) {
                
            } else if (opt::bAllowMissingGenotpyes) {
                // Set genotypes to missing for individuals with too low coverage
                for (int i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                    if (result.counts.depthPerIndividual[i-NUM_NON_GENOTYPE_COLUMNS] < opt::min_depth_in_any_individual) {
                        fields[i][0] = '.'; fields[i][2] = '.';
                    }
                }
                print_vector_stream(fields, std::cout,'\t');
            }
            */
            
        }
    }
    std::cerr << "Total variants filtered out: " << numTotalFilteredOut << std::endl;
    std::cerr << "Reasons for filtering:" << std::endl;
    std::cerr << "Not biallelic: " << numMultiallelic << std::endl;
    std::cerr << "Invariant (not polymorphic in called samples): " << numInvariant << std::endl;
    std::cerr << "Other reasons:" << std::endl;
    std::cerr << "SB - strand bias; OD - overall depth; OQ - overall quality; F - inbreeding coefficient" << std::endl;
    std::cerr << "OQ+OD+F+SB:\t" << numOQ_OD_F_SB << " (" << ((double)numOQ_OD_F_SB/numTotalFilteredOut)*100 << "%)" << std::endl;
    std::cerr << "OQ+OD+F: \t" << numOQ_OD_F << " (" << ((double)numOQ_OD_F/numTotalFilteredOut)*100 << "%)" << std::endl;
    std::cerr << "OQ+OD+SB: \t" << numOQ_OD_SB << " (" << ((double)numOQ_OD_SB/numTotalFilteredOut)*100 << "%)" << std::endl;
    std::cerr << "OQ+SB+F: \t" << numOQ_SB_F << " (" << ((double)numOQ_SB_F/numTotalFilteredOut)*100 << "%)" << std::endl;
    std::cerr << "SB+F+OD: \t" << numSB_F_OD << " (" << ((double)numSB_F_OD/numTotalFilteredOut)*100 << "%)" << std::endl;
    std::cerr << "OQ+OD: \t" << numOQ_OD << " (" << ((double)numOQ_OD/numTotalFilteredOut)*100 << "%)" << std::endl;
    std::cerr << "OQ+F: \t" << numOQ_F << " (" << ((double)numOQ_F/numTotalFilteredOut)*100 << "%)" << std::endl;
    std::cerr << "OQ+SB: \t" << numOQ_SB << " (" << ((double)numOQ_SB/numTotalFilteredOut)*100 << "%)" << std::endl;
    std::cerr << "SB+F: \t" << numSB_F << " (" << ((double)numSB_F/numTotalFilteredOut)*100 << "%)" << std::endl;
    std::cerr << "SB+OD: \t" << numSB_OD << " (" << ((double)numSB_OD/numTotalFilteredOut)*100 << "%)" << std::endl;
    std::cerr << "F+OD: \t" << numF_OD << " (" << ((double)numF_OD/numTotalFilteredOut)*100 << "%)" << std::endl;
    std::cerr << "OQ: \t" << numOQ << " (" << ((double)numOQ/numTotalFilteredOut)*100 << "%)" << std::endl;
    std::cerr << "SB: \t" << numSB << " (" << ((double)numSB/numTotalFilteredOut)*100 << "%)" << std::endl;
    std::cerr << "OD: \t" << numOD << " (" << ((double)numOD/numTotalFilteredOut)*100 << "%)" << std::endl;
    std::cerr << "F: \t" << numF << " (" << ((double)numF/numTotalFilteredOut)*100 << "%)" << std::endl;
    return 0;
}

// filter out sites where more than MAX_NUM_HET individuals are heterozygous
/*
 int num_hets = 0; bool mnhPassed;
 if (opt::max_het_indiv < std::numeric_limits<int>::max()) {
 for (std::vector<std::vector<int> >::size_type i = 0; i < result.counts.individualsWithVariant.size(); i++) {
 if (result.counts.individualsWithVariant[i] == 1)
 num_hets++;
 }
 mnhPassed = (num_hets < opt::max_het_indiv) ? true : false;
 } else {
 mnhPassed = true;
 }
 // bool mnhPassed = testMaxNumHet(result, depthsHetFailed, depthsHetPassed, numVariantsPerHetCount, opt::max_het_indiv,num_indiv_het_vs_depth);
 */


void parseFilterOptions(int argc, char** argv) {
    
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'd': arg >> opt::max_overall_depth; break;
            case 'm': arg >> opt::min_overall_depth; break;
            case 'c': arg >> opt::min_copies; break;
            case 's': arg >> opt::min_depth_in_any_individual; break;
            case '?': die = true; break;
            case OPT_TRIALLELIC: opt::bBiallelicFilter = false; break;
            case OPT_ALLOW_MISSING_GENOTYPES: opt::bAllowMissingGenotpyes = true; break;
            case OPT_MAX_NUM_HET: arg >> opt::max_het_indiv; break;
            case OPT_MQSB: arg >> opt::min_MQSB; break;
            case OPT_FS: arg >> opt::max_FS; break;
            case OPT_STATS: opt::bStats = true; break;
            case OPT_SUBS: arg >> opt::subsample; break;
            case OPT_MIN_F: arg >> opt::min_F; break;
            case OPT_OVERALL_Q: arg >> opt::min_overall_quality; break;
            case OPT_HELP:
                std::cout << FILTER_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    if (argc - optind < 1) {
        std::cerr << "missing arguments\n";
        die = true;
    } 
    else if (argc - optind > 1) 
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (opt::min_F < -1 || opt::min_F > 1) {
        std::cerr << "The values of the inbreeding coefficient F range between -1 and 1 (negative values are likely variant calling errors)\n";
        die = true;
    }
    
    if (opt::min_MQSB < 0 || opt::min_MQSB > 1) {
        std::cerr << "The MQSB p-values range between 0 and 1 (try a value like 0.001)\n";
        die = true;
    }
    
    
    
    
    if (die) {
        std::cout << "\n" << FILTER_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    
}
