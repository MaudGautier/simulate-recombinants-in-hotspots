#ifndef DEF_FUNCTIONS
#define DEF_FUNCTIONS

#include "Libraries.h"


// ========================================================================= //
//                         DECLARATION OF FUNCTIONS                          //
// ========================================================================= //
	
// ---------------------------------------------------------------------------
/**
	Reads a BED file and extracts all pieces of info.
	Uses substr_line_from_bed() to extract info from each line.
	
	@param file_name: Path to the BED file.
	
	@return: vector of Bed_line structure.
*/
vector < Bed_line > bed_file_to_vector (string file_name);


// ---------------------------------------------------------------------------
/**
	Extracts chrom, start, end, ID from a BED line.	
	
	@param line: The line from which info must be extracted.
	@param delim: The delimiter between fields.
	
	@return: the four fields under a Bed_line structure.
*/
Bed_line substr_line_from_bed (string line, string delim = "\t");


// ---------------------------------------------------------------------------
/**
	Extracts info from line.

	@param line: Line from which info must be extracted.
	@param nb_elements: nb_elements to extract.
	@param delim: The delimitor between fields.

	@return: vector of strings (containing all fields).
*/
StringVector extract_info_from_line (
		string line, 
		int nb_elements, 
		string delim = "\t"
		);


// ---------------------------------------------------------------------------
/**
	Creates a vector, for which each field corresponds to one line of the 
	input file.

	@param file_name: Path to the input file.
	
	@return: vector of one string per line.
*/
StringVector create_vector_from_file (string file_name); 


// ---------------------------------------------------------------------------
/** 
	Returns a subselection of the input vector.
	
	@param vec_strings: The input vector (of strings).
	@param selection_size: The nb of elements to be extracted.
	
	@return: A random sample vector.
*/ 
StringVector random_sample_from_string_vector (
		StringVector vec_strings, 
		size_t selection_size
		);


// ---------------------------------------------------------------------------
/**
	Returns a subselection of the input vector, with probabilities given by 
	the weights in vec_weights.

	@param vec_strings: The input vector (of strings).
	@param vec_weights: THe weight associated to each element in vec_strings.
	@param selection_size: The nb of elements to be extracted.
	
	@return: A random sample vector (taking into account the given weights)
*/
StringVector random_sample_from_weighted_string_vector (
		StringVector vec_strings, 
		vector < double > vec_weights, 
		size_t selection_size
		);


// ---------------------------------------------------------------------------
/**
	Returns a vector of pairs containing the string (as key) and the counts
	(as values).

	@param vec_elements: The initial vector of elements.
	@return: The counter map.
*/
CounterMap vector_elements_to_counter_map (StringVector vec_elements);


// ---------------------------------------------------------------------------
/**
	Returns a Param_settings structure of the list of parameters.

	@param line: The string containing all parameters to read.
	@param delim: The delimitor between parameters.

	@return: A structure Param_settings containing all parameter settings.
*/
Param_settings extract_settings_from_line (
		string line, 
		string delim = "\t"
		);


// ---------------------------------------------------------------------------
/**
	Gets all parameters from file and puts them in hash tables.
	
	@param file_parameters: Path to the file containing all parameters
	@param set_hotspot_instances: Set of the names of all hotstpot instances.
								  Passed by reference.
	@param path_to_data: Path to the name of the folder "Data".
	@param delim: The delimitor between parameters.

	@return: A strucutre of type List Hash containing three hash tables:
			 - hash_param_settings:   HASH (paramID --> Param_settings)
			 - hash_param_COs:        HASH (paramID --> CounterMap_COs)
			 - hash_param_NCOs:       HASH (paramID --> CounterMap_NCOs)
*/
List_Hash get_parameters (
		string file_parameters,
		set < string > & set_hotspot_instances,
		string path_to_data,
		string delim = "\t",
		bool verbose = false
		);


// ---------------------------------------------------------------------------
/** 
	Modifies both a vector of list of hotspots (generally empty) and a vector
	of list of weights (generally empty) given in parameters.

	@param list_hotspots: A vector of list of hotspots (generally empty).
						  Passed by reference (to be created).
	@param list_weights: A vector of list of weights (generally empty).
						 Passed by reference (to be created).
	@param file_name: The name of the file (that must contain one column 
					  for the name of the hotspot, and one column for the
					  intensity).
	@param delim: The delimitor between fields in the file.
*/
void modify_weight_vector_from_file (
		StringVector & list_hotspots,
		vector < double > & list_weights,
		string file_name,
		string delim = "\t"
		);


// ---------------------------------------------------------------------------
/**
	Returns a vector of 2 vectors containing the list of hotspots for COs and 
	NCOs.
	
	@param CO_nb: The number of COs required.
	@param NCO_nb: The number of NCOs required.
	@param selection_method: The method for the selection of hotspots.
	@param path_to_data: Path to the name of the folder "Data".
	@param verbose: Whether or not things should be printed.

	@return: A vector containing 2 vectors:
			 - vec_hotspot_ID_COs: The vector of all hotspot IDs for COs.
			 - vec_hotspot_ID_NCOs: The vector of all hotspot IDs for NCOs.
*/
pair < StringVector, StringVector > select_hotspots (
		unsigned int CO_nb,
		unsigned int NCO_nb,
		string selection_method,
		string path_to_data
		);


// ---------------------------------------------------------------------------
/**
	Empties output folders or create them.

	@param output_folder: The name of the folder in which the results will 
						  be written.
	@param hash_param_settings: The HASH map containing settings IDs.
*/
void empty_output_folder (
		string output_folder,
		std::unordered_map < string, Param_settings > hash_param_settings
		);


// ---------------------------------------------------------------------------
/**
	Initialise output files (in all sets of parameters).

	@param output_folder: The name of the folder in which the results will 
						  be written.
	@param hash_param_settings: The HASH map containing settings IDs.
	@param name_output_file: The name of the file that must be opened.
	@param header_vector: list of column names (in a vector format).
*/
void initialise_output_files (
		string output_folder,
		std::unordered_map < string, Param_settings > hash_param_settings,
		StringVector list_output_files,
		StringVector header_vector
		);


// ---------------------------------------------------------------------------
/**
	Writes a line into a file.
	
	@param file_name: The path to the file.
	@param vec_values: The list of values that should be written.
	@param delim: The delimitor between fields.
*/
void write_line_to_file (
		ofstream & file_name,
		StringVector vec_values,
		string delim = "\t"
		);


// ---------------------------------------------------------------------------
/**
	Converts a vector of string into a string with field delimitators.

	@param vec_values: The list of fields that should be outputted.
	@param delim: The delimitor between fields.

	@return: A string containing all fields from the input vector, all
			 separated by the specified delimitor.
*/
// TODO: FAIRE UN TEMPLATE
string vector_to_string (
		StringVector vec_values,
		string delim = ";"
		);
string vector_to_string (
		IntVector vec_values,
		string delim = ";"
		);
string vector_to_string (
		vector < Genotype > vec_values,
		string delim = ";"
		);


// ---------------------------------------------------------------------------
/**
	Counts the occurrences of a particular value in a vector.

	@param vec_values: The vector of all values.
	@param searched_val: The searched value.

	@return: The number of occurrences of the searched value in the list.
*/
int count_in_vector (
		vector < Genotype > vec_values,
		Genotype searched_val
		);


// ---------------------------------------------------------------------------
/**
	Returns the indices of vector all_positions that are found in 
	chosen_positions.

	@param all_positions: The vector of all positions.
	@param chosen_positions: The vector of chosen positions.

	@return: A vector of indices.
*/
IntVector extract_indexes_from_vector (
		IntVector all_positions,
		IntVector chosen_positions
		);


// ---------------------------------------------------------------------------
/**
	Returns a subset of the vec_original containing only positions of the 
	chosen_indices.

	@param vec_original: The original vector that must be subsetted.
	@param chosen_indices: The vector of indices.

	@return: The subsetted vector.
*/
StringVector subset_vector_with_indices (
		StringVector vec_original,
		IntVector chosen_indices
		);


// ---------------------------------------------------------------------------
/**
	Extracts the right alleles.

	@param vec_genotypes: The vector of genotypes (either B6 or CAST).
	@param B6_alleles: The alleles corresponding to the B6 parent.
	@param CAST_alleles: The alleles corresponding to the CAST parent.

	@return: The vector of alleles.
*/
StringVector extract_alleles (
		vector < Genotype > vec_genotypes,
		StringVector B6_alleles,
		StringVector CAST_alleles
		);


// ---------------------------------------------------------------------------
/**
	Define what the mutation types are for all alleles of the fragment.

	@param alleles: The vector of alleles corresponding to the fragment.
	@param B6_alleles: The alleles corresponding to the B6 parent.
	@param CAST_alleles: The alleles corresponding to the CAST parent.

	@return: The vector of mutation types.
*/
StringVector get_mut_types (
		StringVector alleles,
		StringVector B6_alleles,
		StringVector CAST_alleles
		);

// ---------------------------------------------------------------------------
/**
	Creates vector of values to print to output file.

	@param hotspot: Hotspot object of the hotspot corresponding to the
					sequenced fragment.
	@param recomb_event: Recombination_event object of the recombination 
						 event corresponding to the sequenced fragment.
	@param fragment: Fragment object of the fragment corresponding to the
					sequenced fragment.
	@param parameters_name: Name of the parameters written.

	@return: Vector of strings, corresponding to all fields to be printed
			 out.
*/
StringVector create_vec_values_vector (
		Hotspot& hotspot,
		Recombination_event& recomb_event, // TODO: Remove this feature (unused)
		Fragment& fragment,
		string parameters_name
		);


// ---------------------------------------------------------------------------
/**
	Extracts the type of the recombination event (CO or NCO), based on a 
	vector of genotypes.

	@param vec_genotypes: Vector of genotypes containing the genotypes of all
						  sequenced positions.
	
	@return: A type (CO or NCO).
*/
string extract_type_from_genot (vector <Genotype> vec_genotypes);

#endif



