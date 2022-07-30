//
//  process_vcf.cpp
//  C++ project
//
//  Created by Milan Malinsky on 18/04/2013.
//  Copyright (c) 2013 University of Cambridge. All rights reserved.
//

#include <iostream>
#include <stdlib.h>          
#include "process_vcf_stats.h"
#include "process_vcf_filter.h"
#include "process_vcf_utils.h"
#include "process_vcf_variant_sharing.h"
#include "process_vcf_testing.h"
#include "process_vcf_massoko.h"
#include "process_vcf_get_sequences.h"
#include "process_vcf_sequenom.h"
#include "process_vcf_use_map.h"
#include "process_vcf_coding_sequences.h"
#include "evo_shared_variation.h"
#include "process_vcf_search_sex.h"
#include "process_vcf_mt_sequences.h"
#include "process_vcf_stats_testing.h"
#include "process_vcf_reorder.h"
#include "process_vcf_vcf_from_sequenom.h"
#include "process_vcf_fst.h"
#include "process_vcf_merge.h"
#include "process_vcf_cbs.h"
#include "evo_abba_baba.h"
#include "process_vcf_get_aa_seq.h"
#include "process_vcf_fill_aa.h"
#include "process_vcf_join_multiFasta.h"
#include "process_vcf_shortRNA.h"
#include "process_vcf_linkGeneNames.h"
#include "evo_codingStats_from_alignment.h"
#include "evo_DNA_to_Protein.h"
#include "evo_protein_SegregatingSites.h"
#include "evo_codingSeqs_fromGenomes.h"
#include "evo_fullAnnotationExtract.h"
#include "evo_diversity_subsampling.h"
#include "remove_lowercase.h"
#include "evo_permute_codons.h"
#include "evo_Dmin.h"
#include "evo_Dmin_combine.h"
#include "evo_PBS.h"
#include "evo_ABS.h"
#include "evo_FstAgainstAll.h"
#include "evo_combineVCFs.h"
#include "evo_AlleleFeq.h"
#include "evo_agpToNewFasta.h"
#include "evo_distanceToOutgroups.h"
#include "evo_diversityFromHaps.h"
#include "evo_getInformativePairs.h"
#include "evo_findDiscordantPairs.h"
#include "evo_getInformativeReadsFromSam.h"
#include "evo_findDiscordantPairsFromSAM.h"


//#define TESTING 1


#define AUTHOR "Milan Malinsky"
#define PACKAGE_VERSION "0.1 r28"


static const char *VERSION_MESSAGE =
"evo software Version " PACKAGE_VERSION "\n"
"Written by Milan Malinsky.\n"
"\n";

static const char *USAGE_MESSAGE =
"Program: " PROGRAM_BIN "\n"
"Version: " PACKAGE_VERSION "\n"
"Contact: " AUTHOR " [" PACKAGE_BUGREPORT "]\n"
"Usage: " PROGRAM_BIN " <command> [options]\n\n"
"Commands:\n"
"           abba-baba           Perform the abba-baba test\n"
"           alleleFreq          Calculate allele frequencies per population\n"
"           cbs                 calculate the lenghts of tracts of 'compatibility by sequence' (cbs) between samples from genotypes\n"
"           codingStats         claculate various interesting statistics about protein coding sequences given multiple alignment\n"
"           Dmin                calculate the minimum value of the D statistic (ABBA-BABA) for all possible trios of samples\n"
"           DminCombine         combine results from different Dmin runs (e.g. per-chromosome)\n"
"           filter              filter a vcf file\n"
"           FstGlobal           Fst statistic for each population compared against the rest of the dataset\n"
"           getWGSeq            obtain whole genome sequences (e.g. for phylogenetic analyses, recombination estimation etc.)\n"
"           getMtSeq            obtain mitochondrial sequences (e.g. for phylogenetic analyses, recombination estimation etc.)\n"
"           getCodingSeq        obtain coding sequences (e.g. for gene/variant classification - missense/nonsense/etc..)\n"
"           massoko             filtering and getting fixed (and nearly fixed) variants from the Lake Massoko VCF file\n"
"           PBS                 Calculating the population branch statistics (looking for positive selection)\n"
"           RegionsDxy          Dxy (or other statistics) for regions in a .bed file; optional subsampling\n"
"           RegionsPiGeneral    Pi for regions in a .bed file; can deal with multiallelic sites and missigness\n"
"           SeqFromGenomes      Get subsequences (usually genes) from whole genome sequences\n"
"           sequenom            prepare a variant list for a sequenom assay in the format used at the Sanger Institute\n"
"           sharing             find out how Massoko variants segregate in Lake Malawi\n"
"           stats               get various statistics from a vcf file\n"
"           test                used for testing/developing new programs\n"
"           VCFfromSequenom     Generate a VCF file from a sequenom output file\n"
"           vcf-comb            Combine VCF files called against different reference genomes\n"
"           agpToNewFasta       Make a new reference fasta file with rearrangements specified in an agp file\n"
"           DistOutgroups       Calculate the genetic distance of populations to given outgroups\n"
" Various utils:    \n"
"           aa-seq              Final steps in generating an ancestral sequence from a multiple alignment (.maf)\n"
"           aa-fill             Adding the AA (ancestral allele) info field to a VCF (only for SNPs), given an ancestral sequence from aa-seq\n"
"           linkGeneNames                                                                   \n"
"           map                 Use a genetic map to link scaffolds in an assembly\n"
"           merge               Merge two VCF files containing the same variants but different samples\n"
"           reorder             Shuffle columns in a vcf file\n"
"           sharedVariation     search for sharing of heterozygous sites (or variation shared between groups)\n"
"           fst                 Calculating Fst values from VCF, ms simulations, or summarising eigensoft output\n"
"           DNAtoProtein        Translate multiple aligments of gene sequences from DNA into protein\n"
"           ProteinSs           Output segregating sites from a (protein) multiple alignment\n"
"           permuteCodons       Generates new alignment files with a random permutation of codons for NULL coding distribution\n"
"           smallRNA            Generate a distribution showing read lengths and the starting nucleotide for a smallRNA library\n"
"           multi-fasta         A utility tool for dealing with a multi-fasta file (e.g. join all sequences)\n"
"           statsTest           Testing statistical routines in development\n"
"           AnnotationPreformat Extract CDS coordinates from ANNOTATION.gtf for use in 'getCodingSeq'\n"
"           InfoPairs           \n"
"           DiscordantPairs           \n"
"           DiscordantPairsFromSAM \n"
"           InfoReadsSam            \n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

int main(int argc, char **argv) {
    
    if(argc <= 1)
    {
        std::cout << USAGE_MESSAGE;
        return 0;
    }
    else
    {
        std::string command(argv[1]);
        if(command == "help" || command == "--help" || command == "-h")
        {
            std::cout << USAGE_MESSAGE;
            return 0;
        }
        else if(command == "version" || command == "--version")
        {
            std::cout << VERSION_MESSAGE;
            return 0;
        }
        
        if(command == "filter")
            filterMain(argc - 1, argv + 1);
        else if(command == "linkGeneNames")
            linkGNMain(argc - 1, argv + 1);
        else if(command == "remove-lowercase")
            removeLcMain(argc - 1, argv + 1);
        else if(command == "abba-baba")
            abbaBabaMain(argc - 1, argv + 1);
        else if(command == "aa-seq")
            getAaSeqMain(argc - 1, argv + 1);
        else if(command == "aa-fill")
            fillAaMain(argc - 1, argv + 1);
        else if(command == "stats")
            statsMain(argc - 1, argv + 1);
        else if(command == "cbs")
            cbsMain(argc - 1, argv + 1);
        else if(command == "sharing")
            variantSharingMain(argc - 1, argv + 1);
        else if(command == "test")
            testMain(argc - 1, argv + 1);
        else if(command == "massoko")
            massokoMain(argc - 1, argv + 1);
        else if(command == "merge")
            mergeMain(argc - 1, argv + 1);
        else if(command == "multi-fasta")
            processMultiFastaMain(argc - 1, argv + 1);
        else if (command == "fst")
            fstMain(argc - 1, argv + 1);
        else if(command == "getWGSeq")
            getSeqMain(argc - 1, argv + 1);
        else if(command == "getCodingSeq")
            getCodingSeqMain(argc - 1, argv + 1);
        else if(command == "sequenom")
            sequenomMain(argc - 1, argv + 1);
        else if(command == "map")
            mapMain(argc - 1, argv + 1);
        else if(command == "sharedVariation")
            sharedVarMain(argc - 1, argv + 1);
        else if(command == "sex-search")
            sexSearchMain(argc - 1, argv + 1);
        else if(command == "smallRNA")
            shortRNAMain(argc - 1, argv + 1);
        else if (command == "getMtSeq")
            getMtSeqMain(argc - 1, argv + 1);
        else if (command == "statsTest")
            statsTestingMain(argc - 1, argv + 1);
        else if (command == "reorder")
            reorderMain(argc - 1, argv + 1);
        else if (command == "VCFfromSequenom")
            VCFfromSequenomMain(argc - 1, argv + 1);
        else if (command == "codingStats")
            getCodingStats(argc - 1, argv + 1);
        else if (command == "SeqFromGenomes")
            SeqFromGenomes(argc - 1, argv + 1);
        else if (command == "DNAtoProtein")
            DNAtoProtein(argc - 1, argv + 1);
        else if (command == "AnnotationPreformat")
            AnnotationPreformat(argc - 1, argv + 1);
        else if (command == "ProteinSs")
            ProteinSs(argc - 1, argv + 1);
        else if (command == "RegionsDxy")
            subsamplingDxy(argc - 1, argv + 1);
        else if (command == "permuteCodons")
            permuteCodons(argc - 1, argv + 1);
        else if (command == "Dmin")
            DminMain(argc - 1, argv + 1);
        else if (command == "DminCombine")
            DminCombineMain(argc - 1, argv + 1);
        else if (command == "PBS")
            PBSmain(argc - 1, argv + 1);
        else if (command == "ABS")
            ABSmain(argc - 1, argv + 1);
        else if (command == "FstGlobal")
            FstGlobalMain(argc - 1, argv + 1);
        else if (command == "vcf-comb")
            VCFcombMain(argc - 1, argv + 1);
        else if (command == "alleleFreq")
            AFmain(argc - 1, argv + 1);
        else if (command == "agpToNewFasta")
            agpFastaMain(argc - 1, argv + 1);
        else if (command == "DistOutgroups")
            DistOutMain(argc - 1, argv + 1);
        else if (command == "RegionsPiGeneral")
            regionPiMain(argc - 1, argv + 1);
        else if (command == "InfoPairs")
            InfoPairsMain(argc - 1, argv + 1);
        else if (command == "DiscordantPairs")
            DiscordPairsMain(argc - 1, argv + 1);
        else if (command == "InfoReadsSam")
            InfoReadsMain(argc - 1, argv + 1);
        else if (command == "DiscordantPairsFromSAM")
            DiscordPairsFromSAMMain(argc - 1, argv + 1);
        else
        {
            std::cerr << "Unrecognized command: " << command << "\n";
            return 1;
        }
        return 0;
    }
}


