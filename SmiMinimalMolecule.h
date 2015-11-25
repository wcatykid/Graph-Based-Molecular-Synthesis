/*
 *  This file is part of esynth.
 *
 *  esynth is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  esynth is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with esynth.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _SMI_MINIMAL_MOLECULE_GUARD
#define _SMI_MINIMAL_MOLECULE_GUARD 1


#include <string>


#include "Options.h"
#include "SimpleFragmentGraph.h"


//
// An aggregator that contains the minimal amount of information we require per molecule. 
//
class SmiMinimalMolecule : public MinimalMolecule
{
  public:

    std::string smi;
	
    SmiMinimalMolecule(const std::string& smiStr, SimpleFragmentGraph* const g,
                       unsigned short int* const fcounter, unsigned sz) 
        : MinimalMolecule(g, fcounter, sz), smi(smiStr)
    {
    }

    ~SmiMinimalMolecule() { }

    bool equals(MinimalMolecule* const that) const
    {
        //
        // Compare molecules using the SMILES format: string comparison.
        //
        return this->smi == static_cast<SmiMinimalMolecule* const>(that)->smi;
    }

    std::string toString() const
    {
        std::ostringstream oss;

        oss << "Key |" << *this->key << "|" << std::endl;
        oss << "SMI |" << this->smi << "|" << std::endl;

        oss << "Graph: " << this->fingerprint->toString() << std::endl;

        return oss.str();
    }

    friend std::ostream& operator<< (std::ostream& os, const SmiMinimalMolecule& mol)
    {
        os << mol.toString();

        return os;
    }
};

#endif
