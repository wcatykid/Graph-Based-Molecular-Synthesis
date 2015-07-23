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

#ifndef _OBGEN_GUARD
#define _OBGEN_GUARD 1

/**********************************************************************

Adopted this source code from obgen to use NOT as a stand-alone,
but as source code directly linked to the synthetic chemsitry project.


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


#include <string>

#include <openbabel/mol.h>


class OBGen
{
  public:

//    static std::string obgen(std::string& mol);
    static bool obgen(OpenBabel::OBMol* mol);
    static bool fast_obgen(OpenBabel::OBMol* mol);

  private:
    OBGen() {}
};

#endif
