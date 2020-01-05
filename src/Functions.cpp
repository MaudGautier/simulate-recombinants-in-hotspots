#include "Functions.h"


// ========================================================================= //
//                                                                           //
//                          DEFINITION OF FUNCTIONS                          //
//                                                                           //
// ========================================================================= //

// ---------------------------------------------------------------------------
vector < Bed_line > bed_file_to_vector (string file_name)
{
	// Initialise output vector
	vector<Bed_line> vect_out;
	
	// Open file
	ifstream file (file_name, ios::in);
		
	if (file) 
	{       
		// Read file line by line and extract information
		string line;
		while (getline (file, line) ) {
			vect_out.push_back (substr_line_from_bed (line) );
		}
		
		// Closing file
		file.close();
	
	} else { // sinon
		cerr << "ERROR: The file cannot be opened ! "
			<< "(function 'bed_file_to_vector')" << endl;
	}

	return vect_out;
	
}


// ---------------------------------------------------------------------------
Bed_line substr_line_from_bed (string line, string delim) 
{
	// Get chromosome
	size_t pos { line.find_first_of (delim) };
	string chrom { line.substr (0, pos) };
	string rest { line.substr (pos + 1) };
	
	// Get start
	pos = rest.find_first_of (delim);
	string start { rest.substr (0, pos) };
	rest = rest.substr (pos + 1);
	int int_start { stoi (start) };

	// Get end
	pos = rest.find_first_of (delim);
	string end { rest.substr (0, pos) };
	rest = rest.substr (pos + 1);
	int int_end { stoi (end) };

	// Get name
	string name { rest.substr (0, rest.find_first_of (delim) ) };

	return {chrom, int_start, int_end, name};
}


// ---------------------------------------------------------------------------
StringVector extract_info_from_line (
		string line, 
		int nb_elements, 
		string delim) 
{
	// Declaration of variables (useful in the loop)
	size_t pos;
	string value;
	StringVector vect_out;

	for (int i = 0; i < nb_elements; i++) {
		// Extract element
		pos = line.find_first_of (delim);
		value = line.substr (0, pos);
		line = line.substr (pos + 1);
		// Add it to vector
		vect_out.push_back (value);
	}
	return vect_out;
}


// ---------------------------------------------------------------------------
StringVector create_vector_from_file (string file_name)
{
	// Open file
	ifstream file(file_name, ios::in);

	if (file) {
		// Initialise 
		StringVector vec_hotspot_dist; // output vector
		string line;
		string delim { "\t" };
		
		// Add to vector
		while (getline(file, line)) {
			vec_hotspot_dist.push_back (line);
		}

		// return
		return vec_hotspot_dist;

	} else {
		cerr << "ERROR: File missing (Function 'create_vector_from_file').\n";
		throw exception();
	}
}


// ---------------------------------------------------------------------------
StringVector random_sample_from_string_vector (
		StringVector vec_strings, 
		size_t selection_size
		)
{
	// Define type for output
	StringVector chosen_strings;

	// Size of the input vector
	size_t const max_size { vec_strings.size() };

	// Random choice of selection_size
	int rand_index;
	for (int i { 0 }; i < selection_size; ++i) { 
		// Selects index (randomly)
		rand_index = rand() % max_size;
		// Add selected string to vector
	  	chosen_strings.push_back (vec_strings [rand_index] ); 
	}
	return chosen_strings;
}


// ---------------------------------------------------------------------------
StringVector random_sample_from_weighted_string_vector (
		StringVector vec_strings, 
		vector < double > vec_weights, 
		size_t selection_size
		)
{
	// Creation of intervals vector
	IntVector vec_intervals = { 0 };
	for (size_t i { 1 }; i <= vec_strings.size(); i++) {
		vec_intervals.push_back(i);
	}
	// TODO: constructeur avec range plus rapide ?

	// Initialisation of output vector of strings
	StringVector vec_output;

	// Building distribution
	piecewise_constant_distribution<> distribution (
			vec_intervals.begin(), // Address of the first element of vec_intervals
			vec_intervals.end(), // Address of the last element of vec_intervals
			vec_weights.begin() // Address of the first element of vec_weights
			);
	
	// Extraction of selection_size elements
	int rand_index;
	for (int i { 0 }; i < selection_size; i++) {
		rand_index = distribution(GENERATOR);
		vec_output.push_back (vec_strings [rand_index] );
	}
	return vec_output;
}


// ---------------------------------------------------------------------------
CounterMap vector_elements_to_counter_map (
		StringVector vec_elements
		)
{
	// Iterate over all elements to create CounterMap
	CounterMap counts;
	for (int i { 0 }; i < vec_elements.size(); ++i)
	{
		CounterMap::iterator it (counts.find (vec_elements [i] ) );
		// Increment counter if already seen, else: initialise to one.
		if ( it != counts.end() ) {
			it->second++;
		} else {
			counts [vec_elements [i] ] = 1;
		}
	}
	return counts;
}


// ---------------------------------------------------------------------------
Param_settings extract_settings_from_line (string line, string delim) 
{
	
	// 1. Get hotspot selection
	size_t pos { line.find_first_of (delim) };
	string hotspot_sel { line.substr (0, pos) };
	string rest { line.substr (pos + 1)  };
	
	// 2. Get infos for CO
	//	  nb_COs
	pos = rest.find_first_of (delim);
	string val { rest.substr (0, pos) };
	rest = rest.substr (pos + 1);
	unsigned int nb_COs { static_cast < unsigned int > ( stoi (val) ) };
	
	//	  CO_tl_mean
	pos = rest.find_first_of (delim);
	val = rest.substr (0, pos);
	rest = rest.substr (pos + 1);
	double CO_tl_mean { stod (val) };

	//	  CO_tl_sd
	pos = rest.find_first_of (delim);
	val = rest.substr (0, pos);
	rest = rest.substr (pos + 1);
	double CO_tl_sd { stod (val) };
	
	//	  CO_asym
	pos = rest.find_first_of (delim);
	val = rest.substr (0, pos);
	rest = rest.substr (pos + 1);
	double CO_asym { stod (val) };

	// 2. Get infos for NCO
	//	  nb_NCOs
	pos = rest.find_first_of (delim);
	val = rest.substr (0, pos);
	rest = rest.substr (pos + 1);
	unsigned int nb_NCOs { static_cast < unsigned int > ( stoi (val) ) };

	//	  NCO_tl_mean
	pos = rest.find_first_of (delim);
	val = rest.substr (0, pos);
	rest = rest.substr (pos + 1);
	double NCO_tl_mean { stod (val) };

	//	  NCO_tl_sd
	pos = rest.find_first_of (delim);
	val = rest.substr (0, pos);
	rest = rest.substr (pos + 1);
	double NCO_tl_sd { stod (val) };
	
	//	  NCO_asym
	pos = rest.find_first_of (delim);
	val = rest.substr (0, pos);
	rest = rest.substr (pos + 1);
	double NCO_asym { stod (val) };
	
	// 3. Get CAST probability
	pos = rest.find_first_of (delim);
	val = rest.substr (0, pos);
	rest = rest.substr (pos + 1);
	double prob_CAST { stod (val) };
	
	// 4. Number fragments
	pos = rest.find_first_of (delim);
	val = rest.substr (0, pos);
	rest = rest.substr (pos + 1);
	unsigned int nb_frag { static_cast < unsigned int > ( stoi (val) ) };

	// Create output
	Param_settings param_list { 
		hotspot_sel,		
		nb_COs, CO_tl_mean, CO_tl_sd, CO_asym, 
		nb_NCOs, NCO_tl_mean, NCO_tl_sd, NCO_asym,
		prob_CAST, nb_frag
   	};
	return param_list;
}


// ---------------------------------------------------------------------------
List_Hash get_parameters(
		string file_parameters,
		set < string > & set_hotspot_instances,
		string path_to_data,
		string delim, 
		bool verbose
		)
{
	// Initialise hash tables
	std::unordered_map < string, Param_settings > hash_param_settings;
	std::unordered_map < string, CounterMap > hash_param_COs;
	std::unordered_map < string, CounterMap > hash_param_NCOs;

	// Open file
	ifstream file (file_parameters, ios::in);
		
	if (file) 
	{
		// Read file line by line and extract parameters
		string line;
		while (getline (file, line) ) {
			
			// Check that different from header and from endline
			if (line.substr(0,1) != "#" && line != "") {

				// Delineate name of parameters and real parameters
				size_t pos { line.find_first_of (delim) };
				string paramID { line.substr (0, pos) };
			
				// Extract parameters
				string rest_line { line.substr (pos + 1) };
				Param_settings param_list { 
					extract_settings_from_line (rest_line, delim) 
				};
				
				// Add parameters to hash table
				hash_param_settings [ paramID ] = param_list;

				// Define which hotspots should be chosen for each CO and NCO
				pair < StringVector, StringVector > vec_hotspot_IDs_COs_and_NCOs;
				
				// Write the required parameters
				if (verbose) {
					cout << paramID << endl;
					cout << " " << param_list.sel_method << endl;
					cout << " " << param_list.CO_nb 
						<< " " << param_list.CO_tl_mean 
						<< " " << param_list.CO_tl_sd 
						<< " " << param_list.CO_asym << endl;
					cout << " " << param_list.NCO_nb 
						<< " " << param_list.NCO_tl_mean 
						<< " " << param_list.NCO_tl_sd 
						<< " " << param_list.NCO_asym << endl;
					cout << " " << param_list.prob_CAST << endl << endl;
				}
				
				// Select the hotspots for all the numbers of COs and NCOs
				vec_hotspot_IDs_COs_and_NCOs = select_hotspots (
						param_list.CO_nb,
						param_list.NCO_nb,
						param_list.sel_method,
						path_to_data
						);
				
				// Extract the list of hotspots for COs and NCOs separately
				StringVector hotspot_IDs_COs { 
					get <0> (vec_hotspot_IDs_COs_and_NCOs) 
				};
				StringVector hotspot_IDs_NCOs { 
					get <1> (vec_hotspot_IDs_COs_and_NCOs) 
				};
				
				// Add list of hotspot instances to set of hotspots
				copy (	hotspot_IDs_COs.begin(), 
						hotspot_IDs_COs.end(), 
						inserter( set_hotspot_instances, 
								  set_hotspot_instances.end() )  
						);
				copy (	hotspot_IDs_NCOs.begin(), 
						hotspot_IDs_NCOs.end(), 
						inserter( set_hotspot_instances, 
								  set_hotspot_instances.end() )  
						);

				// Create CounterMaps
				CounterMap hotspots_COs_counter { 
					vector_elements_to_counter_map (hotspot_IDs_COs) 
				};
				CounterMap hotspots_NCOs_counter { 
					vector_elements_to_counter_map (hotspot_IDs_NCOs) 
				};
				
				// Add CounterMaps to HASH tables
				hash_param_COs [ paramID ] = hotspots_COs_counter;
				hash_param_NCOs [ paramID ] = hotspots_NCOs_counter;
			}
		}
		
		// Closing file
		file.close();
	
	} else { // sinon
		cerr << "ERROR: The file cannot be opened ! "
			<< "(function 'get_parameters')" << endl;
		throw exception();
	}
	
	// Return
	List_Hash output_list_hash { 
		hash_param_settings, 
		hash_param_COs,
		hash_param_NCOs};

	return output_list_hash;
}


// ---------------------------------------------------------------------------
void modify_weight_vector_from_file (
		StringVector & list_hotspots,
		vector < double > & list_weights,
		string file_name,
		string delim
		) 
{
	// Open file
	ifstream file (file_name, ios::in);
	
	if (file) 
	{
		// Read file line by line and extract parameters
		string line;
		while (getline (file, line) ) {
			
			// Check that different from header and from endline
			if (line.substr(0,1) != "#" && line != "") {

				// Extract hotspot name
				size_t pos { line.find_first_of (delim) };
				string hotspot_name { line.substr (0, pos) };
				list_hotspots.push_back (hotspot_name);

				// Extract intensity
				double intensity { stod (line.substr (pos + 1) ) };
				list_weights.push_back (intensity);
			}
		}
		
		// Closing file
		file.close();
	
	} else { // sinon
		cerr << "ERROR: The file cannot be opened ! "
			<< "(function 'create_weight_vector_from_file')" << endl;
		throw exception();
	}
}


// ---------------------------------------------------------------------------
pair < StringVector, StringVector > select_hotspots (
		unsigned int CO_nb,
		unsigned int NCO_nb,
		string hotspot_selection_method,
		string path_to_data
		)
{
	
	// DECLARATION OF OUTPUTS
	StringVector hotspot_IDs_COs;
	StringVector hotspot_IDs_NCOs;
	

	// -------------------------------------------------------------------- //
	//                                CASE 0:                               //
	//                 ALL HOTSPOTS ON ONE SPECIFIED HOTSPOT                //
	// -------------------------------------------------------------------- //
	
	if (hotspot_selection_method == "TEST") {

		for (int i { 0 }; i < CO_nb; i++) {
			hotspot_IDs_COs.push_back("P9peak.chr15_93973100_93974147");
		}
		for (int i { 0 }; i < NCO_nb; i++) {
			hotspot_IDs_NCOs.push_back("P9peak.chr15_93973100_93974147");
		}

	}
	// TODO: check if other method (more efficient ?)
	

	// -------------------------------------------------------------------- //
	//                 CASE 1: RANDOM SELECTION OF HOTSPOTS                 //
	//     BASED ON REAL DISTRIBUTION OF COs/NCOs (DETECTABILITY-BIASED)    //
	// -------------------------------------------------------------------- //
	
	if (hotspot_selection_method == "observed") {
		
		// Files containing the distributions of COs and NCOs per hotspot
		string CO_distrib { path_to_data + "distribution_CO.txt" };
		string NCO_distrib { path_to_data + "distribution_NCO.txt" };

		// Create Vector of hotspots (based on CO and NCO distributions)
		StringVector hotspot_dist_CO { create_vector_from_file (CO_distrib) };
		StringVector hotspot_dist_NCO { 
			create_vector_from_file (NCO_distrib) 
		};
		
		// Select hotspots for COs and NCOs
		hotspot_IDs_COs = 
			random_sample_from_string_vector ( hotspot_dist_CO, CO_nb );
		hotspot_IDs_NCOs = 
			random_sample_from_string_vector ( hotspot_dist_NCO, NCO_nb );
	
	}
	
	
	// -------------------------------------------------------------------- //
	//                 CASE 2: RANDOM SELECTION OF HOTSPOTS                 //
	//             BASED ON INTENSITY (DETECTABILITY-UNBIASED)              //
	// -------------------------------------------------------------------- //
	
	if (hotspot_selection_method == "expected") {
	
		// TODO: 
		// Créer le fichier des "expected" (basé sur l'intensité des hotspots)
		// Tester la fonction random_sample_from_weighted_string_vector.
		
		// Files containing the distributions of COs and NCOs per hotspot
		string expected_distrib_file { path_to_data + "intensities.txt" };
		// Expected dist: #hotspot_name	  #intensity

		// Create list of hotspots 
		// and weights for hotspots (based on intensities of hotspots)
		StringVector list_hotspots = {};
		vector < double > list_weights = {};
		modify_weight_vector_from_file (
				list_hotspots,
				list_weights, 
				expected_distrib_file,
				"\t"
				);
		
		// Select hotspots for COs and NCOs
		hotspot_IDs_COs = random_sample_from_weighted_string_vector ( 
				list_hotspots, 
				list_weights, 
				CO_nb 
				);
		hotspot_IDs_NCOs = random_sample_from_weighted_string_vector ( 
				list_hotspots, 
				list_weights, 
				NCO_nb 
				);
	
	}
	
	
	// -------------------------------------------------------------------- //
	//                 CASE 3: RANDOM SELECTION OF HOTSPOTS                 //
	//          NO WEIGHT GIVEN (SAME PROBABILITY FOR EACH HOTSPOT)         //
	// -------------------------------------------------------------------- //
	
	if (hotspot_selection_method == "random") {
		
		// Create vector of hotspots (no weight)
		string file_list_hotspots { path_to_data + "list_hotspots.txt" };
		StringVector hotspot_list { 
			create_vector_from_file (file_list_hotspots) 
		};
		
		// Sample from the vector of hotspots
		hotspot_IDs_COs = 
			random_sample_from_string_vector (hotspot_list, CO_nb);
		hotspot_IDs_NCOs = 
			random_sample_from_string_vector (hotspot_list, NCO_nb);
		
	}
	
	return make_pair ( hotspot_IDs_COs, hotspot_IDs_NCOs );
}


// ---------------------------------------------------------------------------
void empty_output_folder(
		string output_folder,
		std::unordered_map < string, Param_settings > hash_param_settings
		)
{
	// 1. Check if output folder exists. If not, create it
	struct stat sb;
	if (!(stat(output_folder.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))) { // Returns 0 if path exists. Then check if path leads to folder.
		mkdir(output_folder.c_str(), 0755); // Creates new folder
	}

	// 2. Add all folders for new settings
	//	  (Remove directory if already exists)
	std::unordered_map < string, Param_settings > :: iterator param_it;
	for (	param_it = hash_param_settings.begin();
			param_it != hash_param_settings.end();
			param_it++
			) {

		// Extract name
		string param_name { param_it->first };
		string folder_param_name { output_folder + param_name };

		// Create folder for param_name
		if (	stat (folder_param_name.c_str(), &sb) == 0 && 
				S_ISDIR(sb.st_mode) && folder_param_name != "/"
		   ) {
			cout << folder_param_name << " already exists. "
				<< "Please remove the folder before carrying on." << endl;
			throw exception();
			//boost::filesystem::remove_all(folder_param_name.c_str());
			//string command_line { "exec rm -r " + folder_param_name
			//system("exec rm -r /pandata/gautier/2_dBGC/5_Simulations_CO_NCO/Results_simultations/*");
		}
		mkdir (folder_param_name.c_str(), 0755);
	}
}


// ---------------------------------------------------------------------------
void initialise_output_files (
		string output_folder,
		std::unordered_map < string, Param_settings > hash_param_settings,
		StringVector list_output_files,
		StringVector header_vector
		)
{
	// Iterate over all parameters
	std::unordered_map < string, Param_settings > :: iterator param_it;
	for (	param_it = hash_param_settings.begin();
			param_it != hash_param_settings.end();
			param_it++	   
		) {

		// Extract name
		string param_name { param_it->first };
		string folder_param_name { output_folder + param_name };

		// Initialise output files
		for (string name_output_file : list_output_files) {
			string path_output_file { folder_param_name + name_output_file };
			ofstream output_file (path_output_file); // Opens in 'write' mode.
			write_line_to_file (
					output_file, 
					header_vector
					);
		}
	}
}


// ---------------------------------------------------------------------------
void write_line_to_file (
		ofstream & output_file,
		StringVector vec_values,
		string delim
		)
{
	
	for (size_t i { 0 }; i < vec_values.size() - 1; i++) {
		output_file << vec_values [i] + "\t";
	}
	output_file << vec_values [vec_values.size() - 1] + "\n";
}


// ---------------------------------------------------------------------------
// TODO: BESOIN DE FAIRE UN TEMPLATE
string vector_to_string (
		StringVector vec_values,
		string delim
		)
{
	// Declare output string
	string output_string { "" };

	// Add vec_values one by one
	for (string el : vec_values) {
		output_string += el + delim;
	}

	if (output_string == "") {
		//cout << "ERROR in function 'vector_to_string': output_string is empty." << endl;
		//throw exception();
		output_string = "NULL";
	}

	return output_string.substr (0, output_string.size() - 1);

}
string vector_to_string (
		IntVector vec_values,
		string delim
		)
{
	// Declare output string
	string output_string { "" };

	// Add vec_values one by one
	for (int el : vec_values) {
		output_string += to_string(el) + delim;
	}

	if (output_string == "") {
		//cout << "ERROR in function 'vector_to_string': output_string is empty." << endl;
		//throw exception();
		output_string = "NULL";
	}
	return output_string.substr (0, output_string.size() - 1 );
}
string vector_to_string (
		vector < Genotype > vec_values,
		string delim
		)
{
	// Declare output string
	string output_string { "" };

	// Add vec_values one by one
	for (Genotype el : vec_values) {
		output_string += Genot_to_string.at(el) + delim;

	}

	if (output_string == "") {
		//cout << "ERROR in function 'vector_to_string': output_string is empty." << endl;
		//throw exception();
		output_string = "NULL";
	}
	return output_string.substr (0, output_string.size() - 1);
}


// ---------------------------------------------------------------------------
int count_in_vector (
		vector < Genotype > vec_values,
		Genotype searched_val
		)
{
	// Add 1 every time you find the value
	int nb { 0 };
	for (Genotype val : vec_values) {
		if (val == searched_val) {
			nb += 1;
		}
	}
	return nb;
}


// ---------------------------------------------------------------------------
IntVector extract_indexes_from_vector (
		IntVector all_positions,
		IntVector chosen_positions
		) 
{
	// Read the "all_positions" vector. 
	// Compare to "chosen_positions". If included => keep. 
	IntVector vec_indices = {};
	for (int index = 0; index < all_positions.size(); index++) {
		if (find (chosen_positions.begin(), 
					chosen_positions.end(), 
					all_positions[index]) != chosen_positions.end() 
				) {
			vec_indices.push_back (index) ;
		}
	}
	return vec_indices;
}


// ---------------------------------------------------------------------------
StringVector subset_vector_with_indices (
		StringVector vec_original,
		IntVector chosen_indices
		)
{
	// Create vector with same size as output.	
	StringVector vec_result (chosen_indices.size(), "");
	
	// Modify vector to get values of vec_originals with indices corresponding
	// to "chosen_indices"
	transform (
			chosen_indices.begin(), 
			chosen_indices.end(), 
			vec_result.begin(), 
			[ vec_original ] (size_t pos) { return vec_original [pos] ; }
			);
	return vec_result;
}


// ---------------------------------------------------------------------------
StringVector extract_alleles (
		vector < Genotype > vec_genotypes, 
		StringVector B6_alleles, 
		StringVector CAST_alleles
		)
{
	// Read vec genotypes and extract the allele (A, T, C, G) corresponding
	// to the right genotype.
	StringVector vect_out;
	for ( size_t index = 0 ; index < vec_genotypes.size() ; index++ ) {
		if (vec_genotypes [index] == B6) {
			vect_out.push_back (B6_alleles [index] ) ;
		}
		else if ( vec_genotypes [index] == CAST ) {
			vect_out.push_back (CAST_alleles [index] ) ;
		}
		else {
			cout << "WARNING in 'extract_alleles' function" << endl;
			throw exception();
		}
	}
	return vect_out;
}


// ---------------------------------------------------------------------------
StringVector get_mut_types (
		StringVector alleles,
		StringVector B6_alleles,
		StringVector CAST_alleles
		)
{
	StringVector vect_out;
	for ( size_t index = 0 ; index < alleles.size() ; index++  ) {
		if (B6_alleles[index].length() == 1 and CAST_alleles[index].length() == 1) {
			vect_out.push_back ("SNP") ;
		} else {
			if (alleles[index].length() > 1) {
				vect_out.push_back ("INS") ;
			} else {
				vect_out.push_back ("DEL") ;
			}
		}
	}
	return vect_out;
}


// ---------------------------------------------------------------------------
// TODO: PROBABLEMENT REVOIR CETTE FONCTION (VOIR BIEN CE QU'ON SORT ETC...)
StringVector create_vec_values_vector (
		Hotspot & hotspot,
		Recombination_event & recomb_event, //TODO: Remove this feature (unused)
		Fragment & fragment,
		string parameters_name
		)
{

	// 1. Hotspot
	Genomic_coordinates hot_coord { hotspot.get_coordinates() };
	string chrom { hot_coord.chrom };
	int hot_start { hot_coord.start };
	int hot_end { hot_coord.end };

	// 2. Variants positions
	IntVector vec_variant_positions { fragment.get_variant_positions() };
	StringVector vec_variant_chromosomes (vec_variant_positions.size(), chrom);
	vector < Genotype > vec_genotypes { fragment.get_genotyped_variants() };
		
	// 3. Alleles
	// Extract from hotspots: reference B6, CAST and positions
	StringVector hotspot_B6_alleles { hotspot.get_B6_alleles() };
	StringVector hotspot_CAST_alleles { hotspot.get_CAST_alleles() };
	IntVector hotspot_position_alleles { hotspot.get_variants() };
	
	// Extract from fragments: only sequenced positions
	IntVector chosen_indices { extract_indexes_from_vector ( 
			hotspot_position_alleles,
			vec_variant_positions 
			) };
	StringVector B6_alleles { subset_vector_with_indices (
			hotspot_B6_alleles, 
			chosen_indices
			) };
	StringVector CAST_alleles { subset_vector_with_indices (
			hotspot_CAST_alleles, 
			chosen_indices
			) };
	StringVector alleles { extract_alleles (
			vec_genotypes, 
			B6_alleles, 
			CAST_alleles
			) };
	StringVector mut_types { get_mut_types (
			alleles,
			B6_alleles,
			CAST_alleles
			) };

	// Write the output string vector with all data
	StringVector vec_values = {
		to_string (fragment.get_counter()), // READ_ID
		to_string (vec_variant_positions.size()), // #VARIANTS
		to_string (count_in_vector (vec_genotypes, B6)), // #B6_VAR
		to_string (count_in_vector (vec_genotypes, CAST)), // #CAST_VAR
		vector_to_string (vec_variant_chromosomes), // CHR
		vector_to_string (vec_variant_positions), // POSITIONS
		vector_to_string (vec_genotypes), // GENOTYPES
		vector_to_string (mut_types), // MUT_TYPES AU PIF ???
		"NA", // QUALITIES AU PIF ???
		vector_to_string (alleles), // ALLELES
		vector_to_string (B6_alleles), // REF_B6_ALLELES
		vector_to_string (CAST_alleles), // REF_CAST_ALLELES
		"NA", // #OVERLAP
		"NA", // VCF_FILTER AU PIF ???
		"NA", // VCF_COVERAGE AU PIF ???
		"NA", // VCF_FREQ AU PIF ???
		hotspot.get_name_to_write(), // TARGET
		chrom + ":" + to_string (hot_start) + "-" + to_string (hot_end), // CHR:START-STOP
		"Simulations_" + parameters_name // SAMPLE_FILE
	};
	return vec_values;
}


// ---------------------------------------------------------------------------
string extract_type_from_genot (
		vector < Genotype > vec_genotypes
		) 
{
	// Count number of changes. 
	// If > 1 => NCO (double-switch), else => CO (single-switch).
	int change = 0;
	auto previous_el = vec_genotypes[0];
	for (auto el : vec_genotypes) {
		if (el != previous_el) {
			change++;
		}
		if (change > 1) {
			return "NCO";
		}
		previous_el = el;
   }
	return "CO";
}

// RAJOUTER
// - sortie avec les "réels"
// - script pour comparer les "réels" et les "obtenus" => avec une métrique pour savoir comment comparer les deux
// - lancer avec des tests de paramètres
//
// A FINIR:
// fonction de choix des hotspots
//
//
// Sortie réels:
// Chaque fragment => est-ce que CO/NCO, longueur exacte
// + sortie réels globale: nombre de CO/NCO + nombre de fragments séquencés + nombre de recombinants + 
// Sortie observés:
// idem.
//
//
