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

#ifndef _VALIDATOR_GUARD
#define _VALIDATOR_GUARD 1


#include <vector>


#include <openbabel/mol.h>


//
// A class to perform validation on the synthesis:
//    Given (1) a hypergraph (formed from linkers and rigids)
//          (2) a set of Open Babel molecules,
//    Verify that the set of molecules are vertices in the hypergraph.
//
class Validator
{
  public:
    Validator(std::vector<OpenBabel::OBMol*>& mols) : molecules(mols) {}
    bool Validate(OpenBabel::OBMol&);
    void Validate(const std::string& fileName);
    void Validate(std::vector<OpenBabel::OBMol>&);

  private:

    std::vector<OpenBabel::OBMol*>& molecules; 
};

#endif
