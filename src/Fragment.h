#ifndef DEF_SEQUENCED_FRAGMENT
#define DEF_SEQUENCED_FRAGMENT

#include "Libraries.h"

class Fragment
{
	private:
		static int counter;
	
	
	public:
		// ================================================================= //
		//                   DECLARATION OF CONSTRUCTORS                     //
		// ================================================================= //
		// (forcement explicite - sinon message)
		
		Fragment(
				Hotspot const& hotspot,
				Recombination_event const& recomb_event,
				unsigned int read_length
				);
		// ~Fragment() { cout << "Fragment détruit." << endl; };

	
	protected:
		// ================================================================= //
		//                    DECLARATION OF ATTRIBUTES                      //
		// ================================================================= //

		// FRAGMENT
		Bed_line fragment_bed;
		int fragment_start; 
		int fragment_end;
		
		// READS
		int r1_start;
		int r1_end;
		int r2_start;
		int r2_end;
		
		unsigned int read_length;
		
		// gamete
		int gamete;

		// Variants
		IntVector variant_positions;
		vector < Genotype > genotyped_variants;
		
		// HOTSPOT
		Hotspot const & ref_hotspot;
		
		// RECOMBINATION EVENT
		Recombination_event const & ref_recomb_event;


	public:
		// ================================================================= //
		//                     DECLARATION OF ACCESSORS                      //
		// ================================================================= //
		
		const string & get_readID() const;
		const IntVector & get_variant_positions() const;
		const vector < Genotype > & get_genotyped_variants() const;
		int get_counter() const;


		// ================================================================= //
		//                      DECLARATION OF METHODS                       //
		// ================================================================= //
		
		// --------------------------------------------------------------------	
		void print_hotspot_info();
		
		// --------------------------------------------------------------------	
		/**
			Chooses the genotype, using a Bernouilli distribution with 
			probability p of choosing CAST.

			@param probability: probability of choosing CAST.
			@return genotype: for the giver.
		*/
		//Genotype choose_giver(double probability);
		
		// --------------------------------------------------------------------	
		/**
			Selects randomly the fragment, among all "real" fragments.
			@param hotspot: Hotspot from which the fragment is selected.
			@return: one particular fragment (with Bed_line coordinates).
		*/
		Bed_line choose_fragment(Hotspot const & hotspot);
		
		// --------------------------------------------------------------------	
		/**
			Selects randomly the gamete (may be 1 or 2).
			@return: The number of strand (1 or 2).
		*/
		int choose_gamete();
		
		// --------------------------------------------------------------------	
		/**
			Extracts variant positions from hotspot.
			@param hotspot : Hotspot to which the fragment belongs.
			@param r1_start: start of read 1.
			@param r1_end  : end of read 1.
			@param r2_start: start of read 2.
			@param r2_end  : end of read 2.
		*/
		IntVector extract_variants_from_hotspot (
				Hotspot const & hotspot,
				int r1_start, int r1_end, 
				int r2_start, int r2_end
				);
		
		// --------------------------------------------------------------------	
		/**
			Extract the exact genotypes, based on the position of variants and 
			conversion tract, the type of recombination product (CO/NCO), the 
			chosen gamete and which gamete was the donor in the conversion 
			event.
			@param variant_positions: List of polymorphic positions within the 
									  sequenced fragment.
			@param gamete			: Number of the gamete (1 or 2).
			@param tract_start		: Start position of the conversion tract.
			@param tract_end		: End position of the conversion tract.
			@param recomb_product	: Type of recombination product (CO or 
									  NCO).
			@param donor			: Donor in the conversion event.
			@return: a vector of genotypes (corresponding to the sequenced 
					 genotypes of the selected fragment).
		*/
		vector <Genotype> extract_genotypes(
				IntVector variant_positions,
				int gamete, 
				int tract_start, 
				int tract_end,
				Recomb_product recomb_product,
				Genotype donor
				);
};
 
#endif

