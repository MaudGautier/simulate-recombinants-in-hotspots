#include "Fragment.h"

// Global variable for the count of fragments
int Fragment::counter(0);

// ========================================================================= //
//                                                                           //
//                        DEFINITION OF CONSTRUCTORS                         //
//                                                                           //
// ========================================================================= //

// ---------------------------------------------------------------------------
Fragment::Fragment(
		Hotspot const & hotspot,
		Recombination_event const & recomb_event,
		unsigned int read_length
		) 
	: fragment_bed { this->choose_fragment (hotspot) },
	fragment_start { this->fragment_bed.start },
	fragment_end { this->fragment_bed.end },
	read_length { read_length },
	r1_start { this->fragment_start },
	r1_end { static_cast < int > (this->fragment_start + read_length - 1) },
	r2_start { static_cast < int > (this->fragment_end - read_length + 1) },
	r2_end { this->fragment_end },
	gamete { this->choose_gamete() },
	variant_positions { this->extract_variants_from_hotspot (
			hotspot, 
			this->r1_start, this->r1_end, 
			this->r2_start, this->r2_end
			) },
	genotyped_variants { this->extract_genotypes (
			this->variant_positions,
			this->gamete,
			recomb_event.get_tract_start(),
			recomb_event.get_tract_start() + recomb_event.get_tract_length(),
			recomb_event.get_recomb_product(),
			recomb_event.get_giver()
			) },
	ref_hotspot { hotspot },
	ref_recomb_event { recomb_event }
{
	counter++;
}




// ========================================================================= //
//                                                                           //
//                          DEFINITION OF ACCESSORS                          //
//                                                                           //
// ========================================================================= //

// ---------------------------------------------------------------------------
const string & Fragment::get_readID() const
{
	return fragment_bed.ID;
}


// ---------------------------------------------------------------------------
const IntVector & Fragment::get_variant_positions() const
{
	return variant_positions;
}


// ---------------------------------------------------------------------------
const vector < Genotype > & Fragment::get_genotyped_variants() const
{
	return genotyped_variants;
}


// ---------------------------------------------------------------------------
int Fragment::get_counter() const
{
	return counter;
}




// ========================================================================= //
//                                                                           //
//                           DEFINITION OF METHODS                           //
//                                                                           //
// ========================================================================= //

// ---------------------------------------------------------------------------
Bed_line Fragment::choose_fragment (
		Hotspot const & hotspot
		) 
{
	// Define list of fragments
	vector < Bed_line > hotspot_fragments { hotspot.get_coord_fragments() };
	
	// Choose randomly within the list
	uniform_int_distribution < int > distribution(0, hotspot_fragments.size() - 1);
	int rand_index = distribution (GENERATOR);
	
	return hotspot_fragments [rand_index];

}


// ---------------------------------------------------------------------------
int Fragment::choose_gamete()
{
	// Uniform distribution to choose between gamete 1 and 2
	uniform_int_distribution < int > distribution (0, 1) ;
	int rand_index = distribution (GENERATOR);
	
	return rand_index;
}


// ---------------------------------------------------------------------------
IntVector Fragment::extract_variants_from_hotspot (
		Hotspot const & hotspot,
		int r1_start, int r1_end,
		int r2_start, int r2_end
		)
{
	// Get all variants in the hotspot
	IntVector variants_of_hotspot { hotspot.get_variants() };
	
	// Extract only variants overlapping either read 1 or read 2
	IntVector variants_out;
	for (
			size_t variant_nb = 0; 
			variant_nb < variants_of_hotspot.size(); 
			variant_nb++
		) {
		if ( ( 
				(variants_of_hotspot [variant_nb] >= r1_start) 
				and (variants_of_hotspot [variant_nb] <= r1_end) 
			 ) or (
				 (variants_of_hotspot [variant_nb] >= r2_start) 
				 and (variants_of_hotspot [variant_nb] <= r2_end)
			 )
		   ) {
			variants_out.push_back (variants_of_hotspot [variant_nb] );
		}
	}
	return variants_out;
}


// ---------------------------------------------------------------------------
vector < Genotype > Fragment::extract_genotypes (
		IntVector variant_positions,
		int gamete,
		int tract_start,
		int tract_end,
		Recomb_product recomb_product,
		Genotype donor
		)
{
	// NOTE: 
	//
	// NON-CROSSOVER - C is DONOR
	// BBBBBB|CCCC|BBBB - gamete 0
	// CCCCCCCCCCCCCCCC - gamete 1
	//
	// NON-CROSSOVER - B is DONOR
	// BBBBBBBBBBBBBBBB - gamete 0
	// CCCCCC|BBBB|CCCC - gamete 1
	//
	//
	// CROSSOVER - C is DONOR
	// BBBBB|CCCCCCCCCC - gamete 0
	// CCCCCCCCCCC|BBBB - gamete 1
	//
	// CROSSOVER - B is DONOR
	// BBBBBBBBBBB|CCCC - gamete 0
	// CCCCC|BBBBBBBBBB - gamete 1
	
	
	// Declare
	Genotype genot_1, genot_2;
	Genotype bef_tract, in_tract, aft_tract;

	// According to the gamete
	if (gamete == 0) {
		genot_1 = B6;
		genot_2 = CAST;
	}
	else if (gamete == 1) {
		genot_1 = CAST;
		genot_2 = B6;
	}
	else {
		cout << "ERROR in 'Fragment::extract_genotypes': gamete number ";
		cout << "is different from both 0 and 1." << endl;
		throw exception();
	}
	
	// According to recomb_product
	if (recomb_product == CO) {
		bef_tract = genot_1 ;
		aft_tract = genot_2 ;
		if (donor == B6) {
			in_tract = B6 ;
		}
		else if (donor == CAST) {
			in_tract = CAST ;
		}
	}
	else if (recomb_product == NCO) {
		bef_tract = genot_1 ;
		aft_tract = genot_1 ;
		if ( ( (donor == CAST) and (gamete == 1) )
			 or 
			 ( (donor == B6) and (gamete == 0) )
		   ) {
			in_tract = genot_1 ;
		}
		else {
			in_tract = genot_2 ;
		}
	}

	// Read all variant positions and extract the righteous genotyoe
	vector < Genotype > vect_genotypes;
	for (const int & variant : variant_positions) {
		if (variant < tract_start) {
			vect_genotypes.push_back(bef_tract);	
		}
		else if (variant > tract_end) {
			vect_genotypes.push_back(aft_tract);
		}
		else {
			vect_genotypes.push_back(in_tract);
		}
	}
	
	return vect_genotypes;
}


