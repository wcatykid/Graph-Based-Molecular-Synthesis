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

/**********************************************************************

Adopted this source code from obgen to use NOT as a stand-alone,
but as source code directly linked to the synthetic chemsitry project.

                                                    C. Alvin 6-20-2014

***********************************************************************/

/**********************************************************************
obgen.cpp - test program for SMILES 3D coordinate generation
          - using systematic rotor search

Copyright (C) 2006 Tim Vandermeersch
Some portions Copyright (C) 2006 Geoffrey R. Hutchison
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/


#include <sstream>
#include <iostream>

#include <openbabel/babelconfig.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/builder.h>


#include "obgen.h"


//
// A lightweight set of classes that ensures all information thrown to the streami
// is never buffered, nor retained.
//

class NullBuffer : public std::streambuf
{
  public:
    int overflow(int c) { return c; }
};

class NullStream : public std::ostream
{
  public:
    NullStream() : std::ostream(&m_sb) {}

  private:
    NullBuffer m_sb;
};


//
// Original OBGEN: Generate rough 3D coordinates for SMILES (or other 0D files).
//
/*
void OBGen::obgen(std::string& smiMol)
{
    OpenBabel::OBConversion conv;
    conv.SetInFormat("SMI");

    OpenBabel::OBMol smiOBMol;
    conv.ReadString(&smiOBMol, smiMol);

    return obgen(smiOBMol);
}
*/

//
// Modified OBGEN: Minimizes SteepestDescent and WeightedRotorSearch steps.
//
bool OBGen::fast_obgen(OpenBabel::OBMol* mol)

{
    std::string ff = "MMFF94";

    OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField(ff);
    if (!pFF)
    {
        std::cerr << "obgen: could not find forcefield '" << ff << "'." << std::endl;
        return false;
    }

    NullStream null_stream;
    pFF->SetLogFile(&null_stream);
    pFF->SetLogLevel(OBFF_LOGLVL_LOW);

    OpenBabel::OBBuilder builder;
    builder.Build(*mol);

    // hydrogens must be added before Setup(mol) is called
    mol->AddHydrogens(false, true);
    if (!pFF->Setup(*mol))
    {
        std::cerr << "obgen: could not setup force field." << std::endl;
        return false;
    }

    pFF->SteepestDescent(1, 1.0e-4);
    pFF->WeightedRotorSearch(250, 1);
    pFF->SteepestDescent(1, 1.0e-6);

    pFF->UpdateCoordinates(*mol);

    return true;
}


//
// Original OBGEN: Generate rough 3D coordinates for SMILES (or other 0D files).
//
bool OBGen::obgen(OpenBabel::OBMol* mol)

{
    std::string ff = "MMFF94";

    OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField(ff);
    if (!pFF)
    {
        std::cerr << "obgen: could not find forcefield '" << ff << "'." << std::endl;
        return false;
    }
      
    NullStream null_stream;
    pFF->SetLogFile(&null_stream);
    pFF->SetLogLevel(OBFF_LOGLVL_LOW);
   
    OpenBabel::OBBuilder builder;
    builder.Build(*mol);

    // hydrogens must be added before Setup(mol) is called 
    mol->AddHydrogens(false, true);
    if (!pFF->Setup(*mol)) 
    {
        std::cerr << "obgen: could not setup force field." << std::endl;
        return false;
    }
 
    pFF->SteepestDescent(500, 1.0e-4); 
    pFF->WeightedRotorSearch(250, 50);
    pFF->SteepestDescent(500, 1.0e-6);

    pFF->UpdateCoordinates(*mol);

return true;

    //
    // Write the molecule to a string in SDF format.
    //
/*
    OpenBabel::OBConversion conv;
    conv.SetOutFormat("SDF");

    std::ostringstream oss;
    conv.Write(&mol, &oss);

    return oss.str();
*/
}
