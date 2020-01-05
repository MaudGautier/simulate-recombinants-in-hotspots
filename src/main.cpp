#include "Libraries.h"

// ==== MAIN ====
int main(int argc, char* argv[])
{
	// ===================================================================== //
	//                                STEP 0:                                // 
	//                                -------                                //
	//																		 //
	//                              PARAMETERS                               //
	//																		 //
	// ===================================================================== //
	
	// Get parameters
	if (argc != 4) {
		cout << "ERROR: The number of arguments must be 3." << endl;
		throw exception();
	}

	// ARG 1 = File Parameters
	// string file_parameters { "./examples/Parameters.txt" };
	string file_parameters { argv[1] };

	// ARG 2 = Path to the Data folder
	// string data_folder { "./examples/data/" };
	string path_to_data { argv[2] };

	// ARG 3 = Path to the output folder
	// string output_folder { "./examples/output/"  };
	string output_folder { argv[3] };
	

	// List of chosen parameters
	cout << "----------------------------------" << endl;
	cout << "CHOSEN PARAMETERS:" << endl;
	cout << "----------------------------------" << endl;
	cout << "File for parameters: " << file_parameters << endl;
	cout << "Path to data       : " << path_to_data << endl;
	cout << "Output folder      : " << output_folder << endl;
	cout << "----------------------------------" << endl;
	cout << endl << endl;
	
	// Verbose
	bool const verbose { true };
	

	// ===================================================================== //
	//                                STEP 0:                                // 
	//                                -------                                //
	//																		 //
	//                 INITIALISATION OF FILE PATH VARIABLES                 //
	//																		 //
	// ===================================================================== //

	// 1. Baits file
	string baits_file { path_to_data + "baits_mm10.onP9.bed" };
	
	// 2. Fragments folder
	string fragments_folder { path_to_data + "Fragments/" };
	
	// 3. Variants folder
	string variants_folder { path_to_data + "Variants_positions/" };
	
	


	// ===================================================================== //
	//                                STEP 1:                                // 
	//                                -------                                //
	//																		 //
	//         READ LIST OF PARAMETERS SETTINGS AND SELECT HOTSPOTS          //
	//																		 //
	// ===================================================================== //
	
	// 1.   Get list of parameter settings
	// 1.a. Read parameters and select hotspots
	set < string > set_hotspot_instances;
	List_Hash list_hashes; // List of hash parameters
	list_hashes = get_parameters (
			file_parameters,
			set_hotspot_instances,
			path_to_data,
			"\t",
			verbose
			);
	
	// 1.b. Extract hash tables containing parameter settings 
	//      and numbers of CO/NCO per hotspot
	std::unordered_map < string, Param_settings > hash_param_settings { 
		list_hashes.hash_param_settings 
	};
	std::unordered_map < string, CounterMap > hash_param_COs { 
		list_hashes.hash_param_COs 
	};
	std::unordered_map < string, CounterMap > hash_param_NCOs { 
		list_hashes.hash_param_NCOs 
	};
	
	// 1.c. Empty output folder
	empty_output_folder (
			output_folder,
			hash_param_settings
			);

	// 1.d. Write header to output files
	//		Initialise lists of fragments
	StringVector list_files = { "/list_recombinants.txt", "/list_fragments.txt" };
	initialise_output_files (
			output_folder,
			hash_param_settings,
			list_files,
			HEADER_FRAGMENTS
			);
	
	//		Initialise evaluation files
	list_files = { "/test_CO_NCO.txt" };
	initialise_output_files (
			output_folder,
			hash_param_settings,
			list_files,
			{ "#Hotspot_ID", "expectation", "obervation" }
			);




	// ===================================================================== //
	//                                STEP 2:                                // 
	//                                -------                                //
	//																		 //
	//              PROCESS EACH HOTSPOT SUCCESSIVELY TO CREATE              //
	//               CO AND NCO MOLECULES AS WELL AS FRAGMENTS               //
	//																		 //
	// ===================================================================== //
	
	// 1. Extract list of hotspot information (chr, start, stop, ID)
	vector < Bed_line > list_hotspot_bed_lines { 
		bed_file_to_vector (baits_file) 
	};
	
	// 2. Print list of hotspot instances (SET) - Verification.
	if (verbose) {
		cout << "----------------------------------" << endl;
		cout << "NB OF HOTSPOTS TO INITIALISE: ";
		cout << set_hotspot_instances.size() << endl;
		cout << "----------------------------------" << endl << endl;
	}
	
	// 3. Read all lines from BED file, but work solely on selected hotspots
	int count { 0 };
	for (auto bed_l : list_hotspot_bed_lines) {

		// Remove non-selected hotspots
		if (set_hotspot_instances.find(bed_l.ID) 
				== set_hotspot_instances.end() ) {
			continue;
		}
		
		// Verbose
		if (verbose) {
			++count;
			cout << "Processing hotspot number " << count
				 << " (" << bed_l.ID << ")..." << endl;
		}
		

		// ----------------------------------------------------------------- //
		//                  SUB-STEP A: INITIALISE HOTSPOT                   //
		// ----------------------------------------------------------------- //
		
		// Get required files for hotspot creation
		string fragments_file { fragments_folder + bed_l.ID }; // Fragments
		string variants_file { variants_folder + bed_l.ID }; // Variants
		
		// Create hotspot
		Hotspot hotspot (
				bed_l.ID, 
				bed_l.chrom, 
				bed_l.start, 
				bed_l.end,
				variants_file,
				fragments_file
				);
		

		// ----------------------------------------------------------------- //
		//                SUB-STEP B: ITERATION ON PARAM_IDs                 //
		// ----------------------------------------------------------------- //
		
		// Iteration on all parameters
		std::unordered_map < string, Param_settings > :: iterator param_it;
		for (	param_it = hash_param_settings.begin();
				param_it != hash_param_settings.end();
				param_it++
			) {
			
			
			// ------------------------------------------------------------- //
			// ---------------------- GET PARAMETERS ----------------------- //
			// ------------------------------------------------------------- //
			
			// Retrieve parameters from HASH table
			Param_settings chosen_parameters { param_it->second };

			// Get number of COs and NCOs for the given hotspot
			unsigned int nb_COs { static_cast < unsigned int > ( 
						hash_param_COs [param_it->first] [bed_l.ID] 
						) };
			unsigned int nb_NCOs { static_cast < unsigned int > ( 
						hash_param_NCOs [param_it->first] [bed_l.ID] 
						) };
			unsigned int nb_frag { chosen_parameters.nb_frag };


			// ------------------------------------------------------------- //
			// --------------------- OPEN OUTPUT FILES --------------------- //
			// ------------------------------------------------------------- //

			// Extract name of the folder
			string param_name { param_it->first };
			string folder_param_name { output_folder + param_name };
			
			// Open file list_recombinants
			string path_output_file_recomb { folder_param_name + "/list_recombinants.txt" };
			ofstream output_file_recomb;
			output_file_recomb.open(path_output_file_recomb, ios::out | ios::app);
			
			// Open file list_fragments
			string path_output_file_frag { folder_param_name + "/list_fragments.txt" };
			ofstream output_file_frag;
			output_file_frag.open(path_output_file_frag, ios::out | ios::app);

			// Open file CO_NCO
			string path_output_file_CO { folder_param_name + "/test_CO_NCO.txt" };
			ofstream output_file_CO;
			output_file_CO.open(path_output_file_CO, ios::out | ios::app);


			// ------------------------------------------------------------- //
			// ------------------ CREATE AND SEQUENCE COs ------------------ //
			// ------------------------------------------------------------- //

			for (int CO_nb = 0; CO_nb < nb_COs; CO_nb++) {

				// Create COs
				Recombination_event crossover (
						hotspot,
						chosen_parameters.CO_asym,
						chosen_parameters.CO_tl_mean,
						chosen_parameters.CO_tl_sd,
						chosen_parameters.prob_CAST,
						CO
						);

				// Sequence COs
				for (int frag_num = 0 ; frag_num < nb_frag ; frag_num++) {
					Fragment fragment (
							hotspot,
							crossover,
							250
							);

					// Extract results (for vec_values)
					StringVector vec_values { create_vec_values_vector (
							hotspot,
							crossover,
							fragment,
							param_name
							) };

					// Write fragments
					write_line_to_file (
							output_file_frag,
							vec_values
							);

					// Write recombinants
					string type = "NA";
					if ( 	(stoi (vec_values[2]) >= 2) 
						and (stoi (vec_values[3]) >= 2)
					   ) {
						write_line_to_file (
								output_file_recomb,
								vec_values
								);

						vector < Genotype > vec_genotypes { 
							fragment.get_genotyped_variants() 
						};
						type = extract_type_from_genot (vec_genotypes);
					}

					// Write CO_file
					StringVector vec_CO_file = { 
						vec_values[0], hotspot.get_ID(), "CO", type 
					};
					write_line_to_file (output_file_CO, vec_CO_file);

				}
			}
			

			// ------------------------------------------------------------- //
			// ----------------- CREATE AND SEQUENCE NCOs ------------------ //
			// ------------------------------------------------------------- //
	
			for (int NCO_nb = 0; NCO_nb < nb_NCOs; NCO_nb++) {

				// Create NCOs
				Recombination_event noncrossover (
						hotspot,
						chosen_parameters.NCO_asym,
						chosen_parameters.NCO_tl_mean,
						chosen_parameters.NCO_tl_sd,
						chosen_parameters.prob_CAST,
						NCO
						);

				// Sequence NCOs
				for (int frag_num = 0 ; frag_num < nb_frag ; frag_num++) {
					Fragment fragment (hotspot,
							noncrossover,
							250
							);

					// Extract results (for vec_values)
					StringVector vec_values { create_vec_values_vector(
							hotspot,
							noncrossover,
							fragment,
							param_name) 
					};

					// Write fragments
					write_line_to_file (
							output_file_frag,
							vec_values
							);

					// Write recombinants
					string type = "NA";
					if (	(stoi(vec_values[2]) >= 2) 
						and (stoi(vec_values[3]) >= 2) 
						) {
						write_line_to_file (
								output_file_recomb,
								vec_values
								);

						vector < Genotype > vec_genotypes { fragment.get_genotyped_variants() };
						type = extract_type_from_genot(vec_genotypes);
					}

					// Write CO_file
					StringVector vec_CO_file = { vec_values[0], hotspot.get_ID(), "NCO", type };
					write_line_to_file(output_file_CO, vec_CO_file);

				}
			}
			

			// ------------------------------------------------------------- //
			// ------------------------ CLOSE FILES ------------------------ //
			// ------------------------------------------------------------- //
	
			output_file_recomb.close();
			output_file_frag.close();
			output_file_CO.close();
		}
	}
	return 0;
}


