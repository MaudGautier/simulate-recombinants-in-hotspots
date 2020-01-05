#include "Hotspot.h"


// ========================================================================= //
//                                                                           //
//                        DEFINITION OF CONSTRUCTORS                         //
//                                                                           //
// ========================================================================= //

// ---------------------------------------------------------------------------
/* // Default constructor */
// Hotspot::Hotspot() {
//     cout << "ATTENTION: Hotspot pas initialisé !" << endl;
// }
//
//
// // ---------------------------------------------------------------------------
// // Constructor for TESTING PURPOSES ONLY
// // A ameliorer en utilisant une méthode pour lire les infos
// Hotspot::Hotspot(string ID) :
//     ID { ID },
//     chr { "chr1" },
//     start { 1 },
//     end { 1000 },
//     variant_positions { {100,200,300,400,500,600,700,800,900} },
//     DSB { 500 },
//     vec_frag_pos {
//         {	{ "f1", 150, 500  },
//             { "f2", 250, 600  },
//             { "f3", 350, 700  },
//             { "f4", 450, 800  },
//             { "f5", 550, 900  }
//         }
//     }
// {
// }
/*  */

// ---------------------------------------------------------------------------
// Constructor that should be used
Hotspot::Hotspot (string ID,
		string chr,
		int start, int end,
		string variants_file,
		string fragments_file,
		int subset_fragments_factor
		) :
	ID { ID },
	chr { chr },
	start { start },
	end { end },
	DSB { start + (end - start)/2 + 1 }, // + 1 ???
	vec_frag_pos {
		extract_fragments_coordinates(
				fragments_file,
				subset_fragments_factor
				)
	}
{
	List_variants list_var { extract_variants(variants_file) };
	variant_positions = list_var.vect_pos; 
	B6_alleles = list_var.vect_B6_alleles;
	CAST_alleles = list_var.vect_CAST_alleles;
}




// ========================================================================= //
//                                                                           //
//                          DEFINITION OF ACCESSORS                          //
//                                                                           //
// ========================================================================= //

// ---------------------------------------------------------------------------
int Hotspot::get_DSB() const {
	return DSB;
}


// ---------------------------------------------------------------------------
Genomic_coordinates Hotspot::get_coordinates() const {
	return {chr, start, end};
}


// ---------------------------------------------------------------------------
const vector<Bed_line>& Hotspot::get_coord_fragments() const {
	return vec_frag_pos;
}


// ---------------------------------------------------------------------------
const IntVector& Hotspot::get_variants() const {
	return variant_positions;
}


// ---------------------------------------------------------------------------
string Hotspot::get_ID() const {
	return ID;
}

// ---------------------------------------------------------------------------
string Hotspot::get_name_to_write() const {
	return ID + "/" + to_string(start) + "/" + to_string(end);
}


// ---------------------------------------------------------------------------
const StringVector& Hotspot::get_B6_alleles() const {
	return B6_alleles;
}


// ---------------------------------------------------------------------------
const StringVector& Hotspot::get_CAST_alleles() const {
	return CAST_alleles;
}




// ========================================================================= //
//                                                                           //
//                           DEFINITION OF METHODS                           //
//                                                                           //
// ========================================================================= //

// ---------------------------------------------------------------------------
vector < Bed_line > Hotspot::extract_fragments_coordinates (
		string file_name,
		int subset_factor
		)
{
	// Use previously written function
	vector < Bed_line > vect_tmp { bed_file_to_vector (file_name) };

	// Subset vector
	vector < Bed_line > vect_out;
	for (int i { 0 }; i < vect_tmp.size(); i += subset_factor) {
		vect_out.push_back (vect_tmp [i]) ;
	}
	
	return vect_out;
}


// ---------------------------------------------------------------------------
void Hotspot::print_fragments(IntVector vec_indexes) 
{
	cout << "List of fragments:\n";
	for (int index : vec_indexes) {
		cout << "ID: " << vec_frag_pos[index].ID << endl;
	}
}


// ---------------------------------------------------------------------------
List_variants Hotspot::extract_variants (
		string file_name
		)
{
	
	// Declare output vectors
	IntVector vect_pos;
	StringVector vect_B6;
	StringVector vect_CAST;
	
	// Read input file
	ifstream file (file_name, ios::in);
	if (file)
	{
		string line;
		while (getline (file, line) ) {
			// Extract pieces of information
			StringVector vect_line { 
				extract_info_from_line (line, 3, "\t") 
			};
			
			// Increment all output vectors
			vect_pos.push_back (static_cast < int > (stoi (vect_line [0])) );
			vect_B6.push_back ( (string) vect_line [1] );
			vect_CAST.push_back ( (string) vect_line [2] );
		}
		return { vect_pos , vect_B6, vect_CAST };
	}
	else {
		cout << "\nERROR: File not found." << endl;
		throw exception();
	}
}
