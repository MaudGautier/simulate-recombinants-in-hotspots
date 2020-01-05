Description
===========

The ``simulate-recombinants-in-hotspots`` repository is a C++ project allowing to simulate a sequencing experiment of recombination products in a defined set of hotspots.

In practice, the simulation is divided into two consecutive steps: 
1. Creating hotspots based on:
		* The position and nature of the variants it hosts
		* The alleles carried by the two parents at each variant
		* The position of the DSB
		* The list of fragments effectively sequenced on that hotspot (from the real experiment) 
2. Creating recombination events with several parameters
		* Hotspot hosting the recombination event
		* Conversion tract position and length
		* Asymetry in the position of the DSB
		* Choice of the parental haplotype that is the donor in the gene conversion event
3. Creating sequenced fragments
		* Choice of the gamete sequenced
		* Exact position (start-stop) of the fragment
		* Read length



Installation
============

Dependencies
------------

Hereunder is the list of dependencies that are necessary for this simulator:

* C++14
* C++ library lboost_filesystem
* C++ library lboost_system


Clone the repository
--------------------

To clone the repository, use this command line:

```
git clone git@github.com:MaudGautier/simulate-recombinants-in-hotspots.git
```



Layout of the repository
========================

Source scripts to compile the simulator are located in the ``src/`` folder.

The executable is in the ``bin/`` folder.

Examples to run the simulator are in the ``examples/`` folder, which contains the ``examples/data/`` folder containing input datasets necessary to the simulator and the ``examples/output/`` folder containing outputs of a simulation run.


Usage
=====

To compile the simulator, use the following command line:

```bash
g++ src/Fragment.cpp src/Functions.cpp src/Hotspot.cpp src/Recombination_event.cpp src/main.cpp src/Global_variables.cpp -o bin/simulator -std=c++14 -lboost_system -lboost_filesystem
```

To run the simulator, use the following command line:

```bash
./bin/simulator ./examples/Parameters.txt ./examples/data/ ./examples/output/
```


Documentation
=============

The documentation for this project is available here.

