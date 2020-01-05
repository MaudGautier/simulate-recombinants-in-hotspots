#ifndef DEF_HOTSPOT
#define DEF_HOTSPOT

#include "Libraries.h"

class Hotspot
{
	public:
		// ================================================================= //
		//                   DECLARATION OF CONSTRUCTORS                     //
		// ================================================================= //
		
		//Hotspot();
		//Hotspot(string ID);
		Hotspot(string ID,
				string chr,
				int start,
				int end,
				string variants_file,
				string fragments_file, int subset_fragments_factor = 1);
		// ~Hotspot() { cout << "Hotspot " << ID << " détruit." << endl; };


	private:
		// ================================================================= //
		//                    DECLARATION OF ATTRIBUTES                      //
		// ================================================================= //

		// ID and Genomic coordinates
		string ID;
		string chr;
		int start;
		int end;

		// SNP positions
		IntVector variant_positions;
		StringVector B6_alleles;
		StringVector CAST_alleles;
		
		// DSB position
		int DSB;

		// Fragments positions
		//Vector_fragments vec_frag_pos;
		vector < Bed_line > vec_frag_pos;

		// Proportion de CO et de NCO à ajouter dedans ?
		

	public:
		// ================================================================= //
		//                     DECLARATION OF ACCESSORS                      //
		// ================================================================= //
		
		int get_DSB() const;
		const vector<Bed_line>& get_coord_fragments() const;
		Genomic_coordinates get_coordinates() const;
		const IntVector& get_variants() const;
		string get_ID() const;
		string get_name_to_write() const;
		const StringVector& get_B6_alleles() const;
		const StringVector& get_CAST_alleles() const;

		
		// ================================================================= //
		//                      DECLARATION OF METHODS                       //
		// ================================================================= //
		
		// --------------------------------------------------------------------	
		/** 
			Reads a bed file containing fragments and puts all genomic 
			coordinates in a vector.

			@param file_name: Path to the file containing all fragments
			@param subset_factor: factor by which the original dataset must 
								  be divided.
								  e.g.: 3 -> 1 line out of 3 will be included
			@return: a vector of Bed_line (corresponding to the ID and 
					 positions of fragments intersecting that particular 
					 hotspot)

		*/
		vector < Bed_line > extract_fragments_coordinates (
			string file_name, 
			int subset_factor = 1
			);
		/*vector<Bed_line> extract_fragments_coordinates(
			string folder, 
			string hotspot_ID, 
			int subset_factor = 1);*/

		// --------------------------------------------------------------------	
		void print_fragments (IntVector vec_indexes);
		
		// --------------------------------------------------------------------	
		/**
			Extracts variants info from file.

			@param file_name: name of the file from which info must be 
				   extracted.
			@return: List_variants of three vectors (position, B6_alleles, 
					 CAST_alleles)
		*/
		List_variants extract_variants (string file_name) ;

};

#endif
