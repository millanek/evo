//
//  evo_Dmin.cpp
//  process_vcf
//
//  Created by Milan Malinsky on 18/06/2018.
//  Copyright © 2018 Milan Malinsky. All rights reserved.
//

#include "evo_Dmin.h"

#define SUBPROGRAM "Dmin"

#define DEBUG 1

static const char *DMIN_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] INPUT_FILE.vcf SETS.txt\n"
"Calculate the Dmin-statistic - the ABBA/BABA stat for all trios of species in the dataset (the outgroup being fixed)"
"the calculation is as definded in Durand et al. 2011"
"The SETS.txt should have two columns: SAMPLE_ID    SPECIES_ID\n"
"The outgroup (can be multiple samples) should be specified by using the keywork Outgroup in place of the SPECIES_ID\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       --AAeqO                                 ancestral allele info in the VCF is from the outgroup (e.g. Pnyererei for Malawi)\n"
"                                               the Outgroup setting in the SETS.txt file will be ignored\n"
"       --fixP3=SPECIES                         (optional) fix the P3 individual and only claculate the stats for cominations of P1 and P2\n"
"       -w SIZE, --window=SIZE                  (optional) output D statistics for nonoverlapping windows containing SIZE SNPs with nonzero D (default: 50)\n"
"       -n, --run-name                          run-name will be included in the output file name\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";


enum { OPT_AA_EQ_O };

static const char* shortopts = "hfw:n:";

static const int JK_WINDOW = 5000;

static const struct option longopts[] = {
    { "run-name",   required_argument, NULL, 'n' },
    { "window",   required_argument, NULL, 'w' },
    { "AAeqO",   no_argument, NULL, OPT_AA_EQ_O },
    { "frequency",   no_argument, NULL, 'f' },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string setsFile;
    static string sampleNameFile;
    static string runName = "";
    static bool bAaEqO = false;
    static int minScLength = 0;
    static int windowSize = 50;
    int jkWindowSize = JK_WINDOW;
}

class GeneralSetCounts {
public:
    GeneralSetCounts(const std::map<string, std::vector<size_t>>& setsToPosMap, const int nSamples) : overall(0) {
        for(std::map<string, std::vector<size_t>>::const_iterator it = setsToPosMap.begin(); it != setsToPosMap.end(); ++it) {
            setRefCounts[it->first] = 0; setAltCounts[it->first] = 0; setAlleleCounts[it->first] = 0;
            setAAFs[it->first] = -1.0; setDAFs[it->first] = -1.0;
            setSizes.push_back(it->second.size());
        }
        individualsWithVariant.assign(nSamples, 0);
    };
    
    int overall;
    std::map<string,int> setRefCounts;
    std::map<string,int> setAltCounts;
    std::map<string,int> setAlleleCounts; // The number of non-missing alleles for this set
    std::vector<size_t> setSizes;
    std::map<string,double> setAAFs; // Allele frequencies - alternative allele
    std::unordered_map<string,double> setDAFs; // Allele frequencies - derived allele
    std::vector<int> individualsWithVariant; // 0 homRef, 1 het, 2 homAlt
   // std::vector<int> set1individualsWithVariant; std::vector<int> set2individualsWithVariant;
   // std::vector<int> set3individualsWithVariant; std::vector<int> set4individualsWithVariant;
};


// Works only on biallelic markers
void getSetVariantCounts(GeneralSetCounts* c, const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap) {
    // std::cerr << fields[0] << "\t" << fields[1] << std::endl;
    
    // Go through the genotypes - only biallelic markers are allowed
    for (std::vector<std::string>::size_type i = 0; i != genotypes.size(); i++) {
        // The first allele in this individual
        if (genotypes[i][0] == '1') {
            c->overall++; c->individualsWithVariant[i]++;
            c->setAltCounts[posToSpeciesMap.at(i)]++; c->setAlleleCounts[posToSpeciesMap.at(i)]++;
        } else if (genotypes[i][0] == '0') {
            c->setAlleleCounts[posToSpeciesMap.at(i)]++; c->setRefCounts[posToSpeciesMap.at(i)]++;
        }
        // The second allele in this individual
        if (genotypes[i][2] == '1') {
            c->overall++;
            c->setAltCounts[posToSpeciesMap.at(i)]++; c->setAlleleCounts[posToSpeciesMap.at(i)]++;
            c->individualsWithVariant[i]++;
        } else if (genotypes[i][0] == '0') {
            c->setAlleleCounts[posToSpeciesMap.at(i)]++; c->setRefCounts[posToSpeciesMap.at(i)]++;
        }
    }
    
    // If at least one of the outgroup individuals has non-missing data
    // Find out what is the "ancestral allele" - i.e. the one more common in the outgroup
    int AAint = -1;
    if (c->setAlleleCounts.at("Outgroup") > 0) {
        if (c->setRefCounts.at("Outgroup") > c->setAltCounts.at("Outgroup")) {
            AAint = 0;
        } else {
            AAint = 1;
        }
    }
    
    // Now fill in the allele frequencies
    for(std::map<string,int>::iterator it = c->setAltCounts.begin(); it != c->setAltCounts.end(); ++it) {
        if (c->setAlleleCounts.at(it->first) > 0) {
            c->setAAFs[it->first] = (double)c->setAltCounts.at(it->first)/c->setAlleleCounts.at(it->first);
            if (AAint == 0) { // Ancestral allele seems to be the ref, so derived is alt
                c->setDAFs[it->first] = (double)c->setAltCounts.at(it->first)/c->setAlleleCounts.at(it->first);
            } else { // Ancestral allele seems to be alt, so derived is ref
                c->setDAFs[it->first] = (double)c->setRefCounts.at(it->first)/c->setAlleleCounts.at(it->first);
            }
        }
    }
}

inline unsigned nChoosek( unsigned n, unsigned k )
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;
    
    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

int DminMain(int argc, char** argv) {
    parseDminOptions(argc, argv);
    string line; // for reading the input files
    
    std::istream* vcfFile = createReader(opt::vcfFile.c_str());
    std::ifstream* setsFile = new std::ifstream(opt::setsFile.c_str());
    string setsFileRoot = stripExtension(opt::setsFile);
    std::ofstream* outFileBBAA = new std::ofstream(setsFileRoot+ "_" + opt::runName + "_BBAA.txt");
    std::ofstream* outFileDmin = new std::ofstream(setsFileRoot+ "_" + opt::runName + "_Dmin.txt");
    std::map<string, std::vector<string>> speciesToIDsMap;
    std::map<string, string> IDsToSpeciesMap;
    std::map<string, std::vector<size_t>> speciesToPosMap;
    std::map<size_t, string> posToSpeciesMap;
    
    // Get the sample sets
    while (getline(*setsFile, line)) {
       // std::cerr << line << std::endl;
        std::vector<string> ID_Species = split(line, '\t');
        speciesToIDsMap[ID_Species[1]].push_back(ID_Species[0]);
        IDsToSpeciesMap[ID_Species[0]] = ID_Species[1];
        std::cerr << ID_Species[1] << "\t" << ID_Species[0] << std::endl;
    }
    // Get a vector of set names (usually species)
    std::vector<string> species;
    for(std::map<string,std::vector<string>>::iterator it = speciesToIDsMap.begin(); it != speciesToIDsMap.end(); ++it) {
        if ((it->first) != "Outgroup") {
            species.push_back(it->first);
            // std::cerr << it->first << std::endl;
        }
    } std::cerr << "There are " << species.size() << " sets (excluding the Outgroup)" << std::endl;
    int nCombinations = nChoosek((int)species.size(),3);
    std::cerr << "Going to calculate " << nCombinations << " Dmin values" << std::endl;
    
    
    // first, get all combinations of three sets (species):
    std::vector<std::vector<string>> trios; trios.resize(nCombinations);
    std::vector<std::vector<int>> triosInt; triosInt.resize(nCombinations);
    std::vector<bool> v(species.size()); std::fill(v.begin(), v.begin() + 3, true); // prepare a selection vector
    int pNum = 0;
    do {
        for (int i = 0; i < v.size(); ++i) {
            if (v[i]) { trios[pNum].push_back(species[i]); triosInt[pNum].push_back(i); }
        } pNum++;
    } while (std::prev_permutation(v.begin(), v.end())); // Getting all permutations of the selection vector - so it selects all combinations
    std::cerr << "Done permutations" << std::endl;
    
    // And need to prepare the vectors to hold the D values:
    std::vector<double> init(3,0.0); // Vector of initial values
    std::vector<std::vector<double>> initDs(3); // vector with three empty (double) vectors
    std::vector<std::vector<double>> Dnums; Dnums.assign(nCombinations,init);
    std::vector<std::vector<double>> Ddenoms; Ddenoms.assign(nCombinations,init);
    std::vector<std::vector<double>> localDnums; localDnums.assign(nCombinations,init);
    std::vector<std::vector<double>> localDdenoms; localDdenoms.assign(nCombinations,init);
    std::vector<std::vector<std::vector<double>>> regionDs; regionDs.assign(nCombinations, initDs);
    std::vector<double> ABBAtotals(nCombinations,0); std::vector<double> BABAtotals(nCombinations,0);
    std::vector<double> BBAAtotals(nCombinations,0);
    std::vector<int> usedVars(nCombinations,0); // Will count the number of used variants for each trio
    int totalVariantNumber = 0;
    std::vector<string> sampleNames; std::vector<std::string> fields;
    // Find out how often to report progress, based on the number of trios
    int reportProgressEvery; if (nCombinations < 1000) reportProgressEvery = 100000;
    else if (nCombinations < 100000) reportProgressEvery = 10000;
    else reportProgressEvery = 1000;
    std::clock_t start; std::clock_t startGettingCounts; std::clock_t startCalculation;
    double durationOverall; double durationGettingCounts; double durationCalculation;
    
    while (getline(*vcfFile, line)) {
        if (line[0] == '#' && line[1] == '#')
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            fields = split(line, '\t');
            std::vector<std::string> sampleNames(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            for (std::vector<std::string>::size_type i = 0; i != sampleNames.size(); i++) {
                posToSpeciesMap[i] = IDsToSpeciesMap[sampleNames[i]];
            }
            // Iterate over all the keys in the map to find the samples in the VCF:
            // Give an error if no sample is found for a species:
            for(std::map<string, std::vector<string>>::iterator it = speciesToIDsMap.begin(); it != speciesToIDsMap.end(); ++it) {
                string sp =  it->first;
                //std::cerr << "sp " << sp << std::endl;
                std::vector<string> IDs = it->second;
                std::vector<size_t> spPos = locateSet(sampleNames, IDs);
                if (spPos.empty()) {
                    std::cerr << "Did not find any samples in the VCF for " << sp << std::endl;
                    assert(!spPos.empty());
                }
                speciesToPosMap[sp] = spPos;
            }
            start = std::clock();
           //  std::cerr << " " << std::endl;
           //  std::cerr << "Outgroup at pos: "; print_vector_stream(speciesToPosMap["Outgroup"], std::cerr);
           //  std::cerr << "telvit at pos: "; print_vector_stream(speciesToPosMap["telvit"], std::cerr);
        } else {
            totalVariantNumber++;
            //if (totalVariantNumber % reportProgressEvery == 0) {
                durationOverall = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
                std::cerr << "Processed " << totalVariantNumber << " variants in " << durationOverall << "secs" << std::endl;
                std::cerr << "GettingCounts " << durationGettingCounts << " calculation " << durationCalculation << "secs" << std::endl;
            //}
            fields = split(line, '\t');
            std::vector<std::string> genotypes(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            //std::vector<std::string> info = split(fields[7], ';');
            string refAllele = fields[3]; if (refAllele.length() > 1) continue; // Only consider biallelic SNPs
            string altAllele = fields[4]; if (altAllele.length() > 1) continue; // Only consider biallelic SNPs
            
            startGettingCounts = std::clock();
            GeneralSetCounts* c = new GeneralSetCounts(speciesToPosMap, (int)genotypes.size());
            getSetVariantCounts(c, genotypes, posToSpeciesMap);
            genotypes.clear(); genotypes.shrink_to_fit();
            durationGettingCounts = ( std::clock() - startGettingCounts ) / (double) CLOCKS_PER_SEC;
            
            
            startCalculation = std::clock();
            double p_O = c->setDAFs.at("Outgroup");
            if (p_O == 0) { delete c; continue; } // We need to make sure that the outgroup is defined
            
            
            std::vector<double> allPs;
            for (std::vector<std::string>::size_type i = 0; i != species.size(); i++) {
                allPs.push_back(c->setDAFs.at(species[i]));
            }
            
            
            // Now calculate the D stats:
            double p_S1; double p_S2; double p_S3;
            for (int i = 0; i != trios.size(); i++) {
                usedVars[i]++;
                
                p_S1 = allPs[triosInt[i][0]];  //pS1test = c->setDAFs.at(trios[i][0]); //assert(p_S1 == pS1test);
                p_S2 = allPs[triosInt[i][1]];
                p_S3 = allPs[triosInt[i][2]];

                if (p_S1 == -1 || p_S2 == -1 || p_S3 == -1)
                    continue; // If any member of the trio has entirely missing data, just move on to the next trio
                
                double ABBA = ((1-p_S1)*p_S2*p_S3*(1-p_O)); ABBAtotals[i] += ABBA;
                double BABA = (p_S1*(1-p_S2)*p_S3*(1-p_O)); BABAtotals[i] += BABA;
                double BBAA = ((1-p_S3)*p_S2*p_S1*(1-p_O)); BBAAtotals[i] += BBAA;
             //   if (p_O == 0.0) {
                Dnums[i][0] += ABBA - BABA;
                Dnums[i][1] += ABBA - BBAA;  // Dnums[i][1] += ((1-p_S1)*p_S3*p_S2*(1-p_O)) - (p_S1*(1-p_S3)*p_S2*(1-p_O));
                Dnums[i][2] += BBAA - BABA;  // Dnums[i][2] += ((1-p_S3)*p_S2*p_S1*(1-p_O)) - (p_S3*(1-p_S2)*p_S1*(1-p_O));

                Ddenoms[i][0] += ABBA + BABA;
                Ddenoms[i][1] += ABBA + BBAA;   // ((1-p_S1)*p_S3*p_S2*(1-p_O)) + (p_S1*(1-p_S3)*p_S2*(1-p_O));
                Ddenoms[i][2] += BBAA + BABA; // ((1-p_S3)*p_S2*p_S1*(1-p_O)) + (p_S3*(1-p_S2)*p_S1*(1-p_O));
                
                localDnums[i][0] += ABBA - BABA; localDnums[i][1] += ABBA - BBAA; localDnums[i][2] += BBAA - BABA;
                localDdenoms[i][0] += ABBA + BABA; localDdenoms[i][1] += ABBA + BBAA; localDdenoms[i][2] += BBAA + BABA;
                if (usedVars[i] % opt::jkWindowSize == 0) {
                    double regionD0 = localDnums[i][0]/localDdenoms[i][0]; double regionD1 = localDnums[i][1]/localDdenoms[i][1];
                    double regionD2 = localDnums[i][2]/localDdenoms[i][2];
                    regionDs[i][0].push_back(regionD0); regionDs[i][1].push_back(regionD1); regionDs[i][2].push_back(regionD2);
                    localDnums[i][0] = 0; localDnums[i][1] = 0; localDnums[i][2] = 0;
                    localDdenoms[i][0] = 0; localDdenoms[i][1] = 0; localDdenoms[i][2] = 0;
                }
            // }
            }
            durationCalculation = ( std::clock() - startCalculation ) / (double) CLOCKS_PER_SEC;
            delete c;
        }
    }
    
    
    for (int i = 0; i != trios.size(); i++) { //
        // Get the standard error values:
        double D1stdErr = jackknive_std_err(regionDs[i][0]); double D2stdErr = jackknive_std_err(regionDs[i][1]);
        double D3stdErr = jackknive_std_err(regionDs[i][2]);
        // Get the D values
        double D1 = Dnums[i][0]/Ddenoms[i][0]; double D2 = Dnums[i][1]/Ddenoms[i][1];
        double D3 = Dnums[i][2]/Ddenoms[i][2];
        // Get the Z-scores
        double D1_Z = abs(D1)/D1stdErr; double D2_Z = abs(D2)/D2stdErr;
        double D3_Z = abs(D3)/D3stdErr;
        
        
        // Find which topology is in agreement with the counts of the BBAA, BABA, and ABBA patterns
        if (BBAAtotals[i] >= BABAtotals[i] && BBAAtotals[i] >= ABBAtotals[i]) {
            if (D1 >= 0)
                *outFileBBAA << trios[i][0] << "\t" << trios[i][1] << "\t" << trios[i][2];
            else
                *outFileBBAA << trios[i][1] << "\t" << trios[i][0] << "\t" << trios[i][2];
            *outFileBBAA << "\t" << abs(D1) << "\t" << D1_Z << "\t";
            *outFileBBAA << BBAAtotals[i] << "\t" << BABAtotals[i] << "\t" << ABBAtotals[i] << std::endl;
        } else if (BABAtotals[i] >= BBAAtotals[i] && BABAtotals[i] >= ABBAtotals[i]) {
            if (D2 >= 0)
                *outFileBBAA << trios[i][0] << "\t" << trios[i][2] << "\t" << trios[i][1];
            else
                *outFileBBAA << trios[i][2] << "\t" << trios[i][0] << "\t" << trios[i][1];
            *outFileBBAA << "\t" << abs(D2) << "\t" << D2_Z << "\t";
            *outFileBBAA << BABAtotals[i] << "\t" << BBAAtotals[i] << "\t" << ABBAtotals[i] << std::endl;
        } else if (ABBAtotals[i] >= BBAAtotals[i] && ABBAtotals[i] >= BABAtotals[i]) {
            if (D3 >= 0)
                *outFileBBAA << trios[i][2] << "\t" << trios[i][1] << "\t" << trios[i][0];
            else
                *outFileBBAA << trios[i][1] << "\t" << trios[i][2] << "\t" << trios[i][0];
            *outFileBBAA << "\t" << abs(D3) << "\t" << D3_Z << "\t";
            *outFileBBAA << ABBAtotals[i] << "\t" << BABAtotals[i] << "\t" << BBAAtotals[i] << std::endl;
        }
        
        // Find Dmin:
        if (abs(D1) <= abs(D2) && abs(D1) <= abs(D3)) { // (P3 == S3)
            if (D1 >= 0)
                *outFileDmin << trios[i][0] << "\t" << trios[i][1] << "\t" << trios[i][2] << "\t" << D1 << "\t" << D1_Z << "\t" << std::endl;
            else
                *outFileDmin << trios[i][1] << "\t" << trios[i][0] << "\t" << trios[i][2] << "\t" << abs(D1) << "\t" << D1_Z << "\t"<< std::endl;
            if (BBAAtotals[i] < BABAtotals[i] || BBAAtotals[i] < ABBAtotals[i])
                std::cerr << "\t" << "WARNING: Dmin tree different from DAF tree" << std::endl;
        } else if (abs(D2) <= abs(D1) && abs(D2) <= abs(D3)) { // (P3 == S2)
            if (D2 >= 0)
                *outFileDmin << trios[i][0] << "\t" << trios[i][2] << "\t" << trios[i][1] << "\t" << D2 << "\t" << D2_Z << "\t"<< std::endl;
            else
                *outFileDmin << trios[i][2] << "\t" << trios[i][0] << "\t" << trios[i][1] << "\t" << abs(D2) << "\t" << D2_Z << "\t"<< std::endl;
            if (BABAtotals[i] < BBAAtotals[i] || BABAtotals[i] < ABBAtotals[i])
                std::cerr << "\t" << "WARNING: Dmin tree different from DAF tree" << std::endl;
        } else if (abs(D3) <= abs(D1) && abs(D3) <= abs(D2)) { // (P3 == S1)
            if (D3 >= 0)
                *outFileDmin << trios[i][2] << "\t" << trios[i][1] << "\t" << trios[i][0] << "\t" << D3 << "\t" << D3_Z << "\t"<< std::endl;
            else
                *outFileDmin << trios[i][1] << "\t" << trios[i][2] << "\t" << trios[i][0] << "\t" << abs(D3) << "\t" << D3_Z << "\t" << std::endl;;
            if (ABBAtotals[i] < BBAAtotals[i] || ABBAtotals[i] < BABAtotals[i])
                std::cerr << "\t" << "WARNING: Dmin tree different from DAF tree" << std::endl;
        }
        //std::cerr << trios[i][0] << "\t" << trios[i][1] << "\t" << trios[i][2] << "\t" << D1 << "\t" << D2 << "\t" << D3 << "\t" << BBAAtotals[i] << "\t" << BABAtotals[i] << "\t" << ABBAtotals[i] << std::endl;
       
    }
    return 0;
    
}



void parseDminOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'w': arg >> opt::windowSize; break;
            case 'n': arg >> opt::runName; break;
            case OPT_AA_EQ_O: opt::bAaEqO = true; break;
            case 'h':
                std::cout << DMIN_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind < 2) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 2)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << DMIN_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::setsFile = argv[optind++];
}
