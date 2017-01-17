
CXXFLAGS=-std=c++11
CXX=g++
BIN := Build
LDFLAGS=-lz

all: $(BIN)/evo

$(BIN)/evo: $(BIN)/process_vcf.o $(BIN)/evo_abba_baba.o $(BIN)/evo_codingStats_from_alignment.o $(BIN)/process_vcf_get_aa_seq.o $(BIN)/process_vcf_cbs.o $(BIN)/process_vcf_fill_aa.o $(BIN)/process_vcf_join_multiFasta.o $(BIN)/process_vcf_utils.o $(BIN)/process_vcf_IUPAC.o $(BIN)/process_vcf_annotation_tools.o $(BIN)/process_vcf_print_routines.o $(BIN)/process_vcf_stats.o $(BIN)/process_vcf_stats_functions.o $(BIN)/process_vcf_filter.o $(BIN)/process_vcf_variant_sharing.o $(BIN)/process_vcf_testing.o $(BIN)/process_vcf_massoko.o $(BIN)/process_vcf_get_sequences.o $(BIN)/process_vcf_coding_sequences.o $(BIN)/process_vcf_sequenom.o $(BIN)/process_vcf_use_map.o $(BIN)/process_vcf_fixed_search.o $(BIN)/process_vcf_search_sex.o $(BIN)/process_vcf_mt_sequences.o $(BIN)/process_vcf_stats_testing.o $(BIN)/process_vcf_reorder.o $(BIN)/process_vcf_vcf_from_sequenom.o $(BIN)/process_vcf_fst.o $(BIN)/process_vcf_merge.o $(BIN)/process_vcf_shortRNA.o process_vcf_linkGeneNames.o $(BIN)/gzstream.o $(BIN)/remove_lowercase.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@ 


#$(BIN)/process_vcf_get_aa_seq.o $(BIN)/process_vcf_cbs.o $(BIN)/process_vcf_fill_aa.o $(BIN)/process_vcf_join_multiFasta.o $(BIN)/process_vcf_utils.o $(BIN)/process_vcf_IUPAC.o $(BIN)/process_vcf_annotation_tools.o $(BIN)/process_vcf_print_routines.o $(BIN)/process_vcf_stats.o $(BIN)/process_vcf_stats_functions.o $(BIN)/process_vcf_filter.o $(BIN)/process_vcf_variant_sharing.o $(BIN)/process_vcf_testing.o $(BIN)/process_vcf_massoko.o $(BIN)/process_vcf_get_sequences.o $(BIN)/process_vcf_coding_sequences.o $(BIN)/process_vcf_sequenom.o $(BIN)/process_vcf_use_map.o $(BIN)/process_vcf_fixed_search.o $(BIN)/process_vcf_search_sex.o $(BIN)/process_vcf_mt_sequences.o $(BIN)/process_vcf_stats_testing.o $(BIN)/process_vcf_reorder.o $(BIN)/process_vcf_vcf_from_sequenom.o $(BIN)/process_vcf_fst.o $(BIN)/process_vcf_merge.o $(BIN)/gzstream.o $(BIN)/remove_lowercase.o: $(BIN)

$(BIN)/%.o: %.cpp 
	$(CXX) -c $(CXXFLAGS) $< -o $@ 

#$(BIN):
#	mkdir -p $@

# Dependencies
$(BIN)/evo: $(BIN)/process_vcf.o $(BIN)/evo_abba_baba.o $(BIN)/evo_codingStats_from_alignment.o $(BIN)/process_vcf_get_aa_seq.o $(BIN)/process_vcf_cbs.o $(BIN)/process_vcf_fill_aa.o $(BIN)/process_vcf_join_multiFasta.o $(BIN)/process_vcf_utils.o $(BIN)/process_vcf_IUPAC.o $(BIN)/process_vcf_annotation_tools.o $(BIN)/process_vcf_print_routines.o $(BIN)/process_vcf_stats.o $(BIN)/process_vcf_stats_functions.o $(BIN)/process_vcf_filter.o $(BIN)/process_vcf_variant_sharing.o $(BIN)/process_vcf_testing.o $(BIN)/process_vcf_massoko.o $(BIN)/process_vcf_get_sequences.o $(BIN)/process_vcf_coding_sequences.o $(BIN)/process_vcf_sequenom.o $(BIN)/process_vcf_use_map.o $(BIN)/process_vcf_fixed_search.o $(BIN)/process_vcf_search_sex.o $(BIN)/process_vcf_mt_sequences.o $(BIN)/process_vcf_stats_testing.o $(BIN)/process_vcf_reorder.o $(BIN)/process_vcf_vcf_from_sequenom.o $(BIN)/process_vcf_fst.o $(BIN)/process_vcf_merge.o $(BIN)/gzstream.o $(BIN)/remove_lowercase.o | $(BIN)


