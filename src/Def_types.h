#ifndef NEW_TYPES
#define NEW_TYPES


// ------------------------------------------------------------------------- //
//                                ENUMERATORS                                //
// ------------------------------------------------------------------------- //

// Genotype
enum Genotype { B6, CAST };
const boost::unordered_map < Genotype, const char* > 
  Genot_to_string = map_list_of
	(B6, "B6")
	(CAST, "CAST");
const boost::unordered_map < const char*, Genotype > 
  String_to_genot = map_list_of
	("B6", B6)
	("CAST", CAST);

// Recombinant product
enum Recomb_product { CO, NCO };
const boost::unordered_map < Recomb_product, const char* > 
  Recomb_to_string = map_list_of
	(CO, "CO")
	(NCO, "NCO");


// ------------------------------------------------------------------------- //
//                                  TYPEDEF                                  //
// ------------------------------------------------------------------------- //

// Special vectors
typedef vector < int > IntVector; // Of int
typedef vector < string > StringVector; // Of string
struct Coord_fragments;
typedef vector < Coord_fragments > Vector_fragments; // of coord_fragments

// Special hash maps
typedef map < string, unsigned int > CounterMap;


// ------------------------------------------------------------------------- //
//                                STRUCTURES                                 //
// ------------------------------------------------------------------------- //

// Fragments: name, start, end
struct Coord_fragments { 
	string name;
	int start;
	int end;
};

// Genomic coordinates: chr, start, stop
struct Genomic_coordinates {
	string chrom;
	int start;
	int end;
};

// Bed line: chr, start, stop, ID
struct Bed_line { 
	string chrom;
	int start;
	int end;
	string ID;
};

// Param settings: selection method, 
//				   #CO, CO mean length, CO sd length, CO asymetry,
//				   #NCO, NCO mean length, NCO sd length, NCO asymetry
//				   probability of CAST donor.
struct Param_settings {
	string sel_method;
	unsigned int CO_nb;
	double CO_tl_mean;
	double CO_tl_sd;
	double CO_asym;
	unsigned int NCO_nb;
	double NCO_tl_mean;
	double NCO_tl_sd;
	double NCO_asym;
	double prob_CAST;
	unsigned int nb_frag;
};

// List hash tables: 1) Hash for parameter settings (Param_settings struct)
//					 2) Hash for CO counts (CounterMap typedef)
//					 3) Hash for NCO counts (CounterMap typedef)
struct List_Hash {
	std::unordered_map < string, Param_settings > hash_param_settings;
	std::unordered_map < string, CounterMap > hash_param_COs;
	std::unordered_map < string, CounterMap > hash_param_NCOs;
};

// List variants
struct List_variants {
	IntVector vect_pos;
	vector < string > vect_B6_alleles;
	vector < string > vect_CAST_alleles;
};

#endif
