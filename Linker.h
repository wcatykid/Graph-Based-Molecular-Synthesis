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

#ifndef _LINKER_GUARD
#define _LINKER_GUARD 1

#include <vector>

#include <openbabel/mol.h>

#include "EdgeAggregator.h"
#include "Molecule.h"

class Linker : public Molecule
{
  public:
    Linker(OpenBabel::OBMol*, const std::string& name);
    Linker() {}

    ~Linker() {}

    virtual bool IsLinker() const { return true; }
    virtual bool IsComplex() const { return false; }
    virtual bool IsRigid() const { return false; }

    bool operator==(const Linker& that) const
    {
        return this->getUniqueIndexID() == that.getUniqueIndexID();
    }

  protected:
    virtual void parseAppendix(std::string& suffix, int numAtoms = -1);
};

#endif
