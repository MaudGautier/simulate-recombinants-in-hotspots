#include "Recombination_event.h"

int Recombination_event::counter(0);

// ========================================================================= //
//                                                                           //
//                        DEFINITION OF CONSTRUCTORS                         //
//                                                                           //
// ========================================================================= //

// ---------------------------------------------------------------------------
/* Recombination_event::Recombination_event() :  */
//     ref_hotspot { Hotspot() }
// {
//     cout << "ATTENTION: Recombinant pas initialisé !" << endl;
/* } */


// ---------------------------------------------------------------------------
/* Recombination_event::Recombination_event(Hotspot const& hotspot, */
	//     int fragment_start, int fragment_end,
	//     int read_length,
	//     double asymetrie, int longueur_tract,
	//     Genotype donneur
	//     )
	// : fragment_start(fragment_start),
	// fragment_end(fragment_end),
	// r1_start(fragment_start),
	// r1_end(r1_start + read_length - 1),
	// r2_start(fragment_end - read_length + 1),
	// r2_end(fragment_end),
	// longueur_tract(longueur_tract),
	// tract_start( hotspot.get_DSB() - longueur_tract * asymetrie ),
	// tract_end(tract_start + longueur_tract - 1),
	// donneur(donneur),
	// pt_hotspot(&hotspot)
    //
/* {} */


// ---------------------------------------------------------------------------
Recombination_event::Recombination_event(Hotspot const& hotspot,
		double asymmetry, double param1, double param2,
		double probability_CAST, Recomb_product recomb_product
		) 
	: tract_length { this->choose_length(param1, param2, recomb_product) },
	tract_start { static_cast<int>(hotspot.get_DSB() - this->tract_length * asymmetry) },
	tract_end { tract_start + tract_length - 1 },
	giver { this->choose_giver(probability_CAST) },
	//pt_hotspot(&hotspot),
	ref_hotspot { hotspot },
	recomb_product { recomb_product }
{
	counter++;
	// Initialise tract_length based on the boolean (CO or not CO ?)
	
	//tract_length = 
	
	// Initialise tract_start and tract_end
	//tract_start = hotspot.get_DSB() - tract_length * asymmetry;
	//tract_end = tract_start + tract_length - 1;
	
	// Choose the giver (based on the hotspot)
	// Method to choose giver
	

	
}





// ========================================================================= //
//                                                                           //
//                          DEFINITION OF ACCESSORS                          //
//                                                                           //
// ========================================================================= //

// ---------------------------------------------------------------------------
int Recombination_event::get_tract_length() const
{
	return tract_length;
}


// ---------------------------------------------------------------------------
int Recombination_event::get_tract_start() const
{
	return tract_start;
}


// ---------------------------------------------------------------------------
Genotype Recombination_event::get_giver() const 
{
	return giver;
}


// ---------------------------------------------------------------------------
Recomb_product Recombination_event::get_recomb_product() const
{
	return recomb_product;
}


// ---------------------------------------------------------------------------
int Recombination_event::get_counter() const
{
	return counter;
}




// ========================================================================= //
//                                                                           //
//                           DEFINITION OF METHODS                           //
//                                                                           //
// ========================================================================= //

// ---------------------------------------------------------------------------
void Recombination_event::print_hotspot_info()
{
	cout << "DSB: " << ref_hotspot.get_DSB() << endl;
}


// ---------------------------------------------------------------------------
int Recombination_event::choose_length (
		double tract_length_mean, 
		double tract_length_sd, 
		Recomb_product recomb_product
		)
{ 
	// If CO: choose length of CO from a normal distribution.
	if (recomb_product == CO) {
		normal_distribution < double > distribution (
				tract_length_mean, 
				tract_length_sd
				);
		double chosen_length { distribution (GENERATOR) };
		return (int) chosen_length;
	}

	// If NCO: choose length of NCO from a gamma distribution.
	else if (recomb_product == NCO) {
		// Define parameters for the gamma distribution
		//The estimates obtained this way are method of moments estimates. In particular, we know that E(X)=αθE(X)=αθ and Var[X]=αθ2Var[X]=αθ2 for a gamma distribution with shape parameter αα and scale parameter θθ (see wikipedia). Solving these equations for αα and θθ yields α=E[X]2/Var[X]α=E[X]2/Var[X] and θ=Var[X]/E[X]θ=Var[X]/E[X]. Now substitute the sample estimates to obtain the method of moments estimates α̂ =x¯2/s2α^=x¯2/s2 and θ̂ =s2/x¯θ^=s2/x¯.
		// Estimate alpha and beta parameters
		// alpha = E(X)**2/Var(X)
		// beta = Var(X)/E(X)
		double alpha { pow (tract_length_mean, 2) / pow (tract_length_sd, 2) };
		double beta { pow (tract_length_sd, 2) / tract_length_mean };

		gamma_distribution < double > distribution (alpha, beta) ;
		double chosen_length { distribution (GENERATOR) };
		return (int) chosen_length;
	}

	// If none: Error
	else {
		cout << "\nERROR: wrong type in Recombination_event::choose_length\n";
		throw exception();
	}
}


// ---------------------------------------------------------------------------
Genotype Recombination_event::choose_giver (
		double probability
		)
{
	// Random generate from the binomial distribution
	bernoulli_distribution distribution (probability);
	if (distribution (GENERATOR) ) {
		return CAST;
	}
	return B6;
	// ATTENTION AJOUTER: SI -1 => LIRE DANS FICHIER "EXPECTED"
}


