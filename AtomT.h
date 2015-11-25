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

#ifndef _ATOM_TYPE_GUARD
#define _ATOM_TYPE_GUARD 1

#include <string>
#include <cctype>
#include <cstdlib>

const unsigned ATOM_T_CARBON     = 0;
const unsigned ATOM_T_CHLORINE   = 1;
const unsigned ATOM_T_HYDROGEN   = 2;
const unsigned ATOM_T_NITROGEN   = 3;
const unsigned ATOM_T_OXYGEN     = 4;
const unsigned ATOM_T_PHOSPHORUS = 5;
const unsigned ATOM_T_SULFUR     = 6;
const unsigned ATOM_T_FLUORINE   = 7;
const unsigned ATOM_T_BROMINE    = 8;
const unsigned ATOM_T_BORON      = 9;
const unsigned ATOM_T_IODINE     = 10;
const unsigned ATOM_T_UNKNOWN    = 11;


const unsigned SPECIAL_T_AM       = 0;
const unsigned SPECIAL_T_AROMATIC = 1;
const unsigned SPECIAL_T_NONE     = 2;
const unsigned SPECIAL_T_PL       = 3;
const unsigned SPECIAL_T_CO       = 4;
const unsigned SPECIAL_T_O        = 5;
const unsigned SPECIAL_T_CAT      = 6;

typedef unsigned AtomEnumT;
typedef unsigned SpecialEnumT;

//
// Light-weight aggregator for Atoms and their types
//
class AtomT
{
  public:
    AtomT(const AtomT* const);
    AtomT(AtomEnumT type = ATOM_T_UNKNOWN, int val = -1,
          SpecialEnumT spec = SPECIAL_T_NONE);
    AtomT(const std::string&);
 
    std::string toString() const;
    std::string getAtomType() const;
    friend std::ostream& operator<< (std::ostream& os, const AtomT& atomType);
    bool operator==(const AtomT& that) const;
    bool operator!=(const AtomT& that) const { return !(*this == that); }

  private:
    static AtomEnumT convertToAtomEnum(const std::string&);
    static SpecialEnumT convertToSpecialEnum(const std::string&);

    typedef struct CompressedAtomT
    {
        AtomEnumT atomType : 8;
        int specificNum   : 16;
        SpecialEnumT specialT : 8;
    } SmallAtomT;

    SmallAtomT theAtomT;
};

#endif
