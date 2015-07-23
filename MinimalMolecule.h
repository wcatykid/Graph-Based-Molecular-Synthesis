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

#ifndef _MINIMAL_MOLECULE_GUARD
#define _MINIMAL_MOLECULE_GUARD 1


#include <string>


#include "Options.h"
#include "SimpleFragmentGraph.h"


//
// An aggregator that contains the minimal amount of information we require per molecule. 
//
class MinimalMolecule
{
  public:

    // More compact, string representation of the counters of the fragment
    // used in the molecule.
    char* key;
    SimpleFragmentGraph* fingerprint;
    unsigned int uniqueIndexID;
	
    MinimalMolecule(SimpleFragmentGraph* const g,
                    unsigned short int* const fcounter, unsigned sz) : fingerprint(g)
    {
        // Create the key value for this molecule.
        // Do so by converting each counter to a two-digit string
        key = new char[2 * sz + 1];
        for (int i = 0, kcounter = 0; i < sz; i++)
        {
            key[kcounter++] = intToChar(fcounter[i] / 10); // units digit
            key[kcounter++] = intToChar(fcounter[i]); // tens digit
        }

        // Indicate end of string
        key[2 * sz] = '\0';
    }
      
    ~MinimalMolecule()
    {
        delete[] key;
        delete fingerprint;
    }

    virtual bool equals(MinimalMolecule* const that) const
    {
        //
        // The fragment counter maintains the number of instances of each specific fragment;
        // if any of those counts differ, we have non-isomorphism.
        //
        if (strcmp(this->key, that->key) != 0) return false;

        //
        // If we reach this point in the code, we can expect the two molecules to have the
        // same number of (1) linkers, (2) rigids, (3) unique rigids, (4) unique linkers,
        // (5) bonds, and (6) # atoms
        //
        //
        // Fingerprint verification.
        //
        // Fingerprint checking is last since it is slow; check other characteristics first.
        //
        return this->fingerprint->IsIsomorphicTo(that->fingerprint);
    }

    virtual std::string toString() const
    {
        std::ostringstream oss;

        oss << "Key |" << *this->key << "|" << std::endl;

        oss << "Graph: " << this->fingerprint->toString() << std::endl;

        return oss.str();
    }

    friend std::ostream& operator<< (std::ostream& os, const MinimalMolecule& mol)
    {
        os << mol.toString();

        return os;
    }
    
  protected:
    unsigned char intToChar(unsigned int i)
    {
        switch(i % 10)
        {
            case 0: return '0';
            case 1: return '1';
            case 2: return '2';
            case 3: return '3';
            case 4: return '4';
            case 5: return '5';
            case 6: return '6';
            case 7: return '7';
            case 8: return '8';
            case 9: return '9';
        }
    }
};
#endif
