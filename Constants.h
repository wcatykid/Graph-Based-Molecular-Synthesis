/*
 *  This file is part of synth.
 *
 *  synth is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  synth is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with synth.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _CONSTANTS_GUARD
#define _CONSTANTS_GUARD 1


#include <string>


const unsigned int null = 0;

const std::string COMPLIANT_EXE = "compliantwriter";

const int NOT_FOUND = -1;


// Debugging constants
const bool DEBUG = true;
const bool HYPERGRAPH_CONSTR_DEBUG = true;
const bool g_debug_output = false;


const int THREAD_POOL_SIZE = 10;


// skip the entire synthesis, just output lipinski descriptors for
//  the input fragments to "initial_fragments_logfile.txt" and exit
const bool g_calculate_lipinski_descriptors_for_input_fragments_only = false;

// upper bound for level-based threading synthesis
extern unsigned int HIERARCHICAL_LEVEL_BOUND;

// Limiting factor on molecule generation.
extern double MOLWT_UPPERBOUND;
extern double HBD_UPPERBOUND;
extern double HBA1_UPPERBOUND;
extern double LOGP_UPPERBOUND;

extern std::string DEFAULT_OUTPUT_DIR;

extern bool VALIDATE;

#endif
