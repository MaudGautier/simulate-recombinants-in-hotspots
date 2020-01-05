#ifndef DEF_RECOMBINANT_MOLECULE
#define DEF_RECOMBINANT_MOLECULE

#include "Libraries.h"

class Recombination_event
{
	private:
		static int counter;
	
	
	public:
		// ================================================================= //
		//                   DECLARATION OF CONSTRUCTORS                     //
		// ================================================================= //
		// (forcement explicite - sinon message)
		
		// Recombination_event();
	    /*  Recombination_event( */
				// Hotspot const& hotspot,
				// int fragment_start, int fragment_end,
				// int read_length,
				// double asymetrie, int longueur_tract,
				// Genotype donneur
				/* ); */
		/* Recombination_event( */
				// Hotspot const& hotspot,
				// double asymmetry, double param1, double param2,
				// Genotype giver, string type
				/* ); */
		Recombination_event(
				Hotspot const& hotspot,
				double asymmetry, double param1, double param2,
				double probability_CAST, Recomb_product recomb_product
				);


	protected:
		// ================================================================= //
		//                    DECLARATION OF ATTRIBUTES                      //
		// ================================================================= //

		// CONVERSION TRACT
		int tract_length;
		int tract_start;
		int tract_end;
		Genotype giver;
		double asymmetry;
		Recomb_product recomb_product;

		// HOTSPOT
		Hotspot const& ref_hotspot;


	public:
		// ================================================================= //
		//                     DECLARATION OF ACCESSORS                      //
		// ================================================================= //
		
		int get_tract_length() const;
		int get_tract_start() const;
		Genotype get_giver() const;
		Recomb_product get_recomb_product() const;
		int get_counter() const;


		// ================================================================= //
		//                      DECLARATION OF METHODS                       //
		// ================================================================= //
		
		void print_hotspot_info();
		
		// -------------------------------------------------------------------
		/**
			Chooses length of the conversion tract, depending on type of event
			(CO or NCO), and mean and sd tract lengths for the conversion 
			tract product.
			
			@param tract_length_mean: Mean tract length.
			@param tract_length_sd  : Standard error of the tract length.
			@param recomb_product   : CO or NCO.
			
			@return: length of the conversion tract.
		*/
		int choose_length (
				double tract_length_mean, 
				double tract_length_sd, 
				Recomb_product recomb_product
				);
		
		// -------------------------------------------------------------------
		/**
			Chooses the genotype, using a Bernouilli distribution with 
			probability p of choosing CAST.

			@param probability: probability of choosing CAST.
		    @return genotype: for the giver.
		*/
		Genotype choose_giver (double probability) ;

};
 
#endif

