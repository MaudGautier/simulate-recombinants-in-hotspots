#include "Global_variables.h"

//default_random_engine GENERATOR;
//mt19937 GENERATOR(1); // SEED
mt19937 GENERATOR(time(NULL)); // RANDOM

StringVector HEADER_FRAGMENTS = {
	"#READ_ID",
	"#VARIANTS",
	"#B6_VAR",
	"#CAST_VAR",
	"CHR",
	"POSITIONS",
	"GENOTYPES",
	"MUT_TYPES",
	"QUALITIES",
	"ALLELES",
	"REF_B6_ALLELES",
	"REF_CAST_ALLELES",
	"#OVERLAP",
	"VCF_FILTER",
	"VCF_COVERAGE",
	"VCF_FREQ",
	"TARGET",
	"CHR_START_STOP",
	"SAMPLE_FILE"
};