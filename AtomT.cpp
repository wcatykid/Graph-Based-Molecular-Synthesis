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

#include <string>
#include <sstream>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <algorithm>


#include "AtomT.h"


AtomT::AtomT(const AtomT* const that)
{
    this->theAtomT.atomType = that->theAtomT.atomType;
    this->theAtomT.specificNum = that->theAtomT.specificNum;
    this->theAtomT.specialT = that->theAtomT.specialT;
}


AtomT::AtomT(AtomEnumT type, int val, SpecialEnumT spec)
{
    this->theAtomT.atomType = type;
    this->theAtomT.specificNum = val;
    this->theAtomT.specialT = spec;
}

AtomEnumT AtomT::convertToAtomEnum(const std::string& s)
{
    if (s == "C")       return ATOM_T_CARBON;
    else if (s == "Cl") return ATOM_T_CHLORINE;
    else if (s == "H")  return ATOM_T_HYDROGEN;
    else if (s == "N")  return ATOM_T_NITROGEN;
    else if (s == "O")  return ATOM_T_OXYGEN;
    else if (s == "P")  return ATOM_T_PHOSPHORUS;
    else if (s == "S")  return ATOM_T_SULFUR;
    else if (s == "F")  return ATOM_T_FLUORINE;
    else if (s == "Br") return ATOM_T_BROMINE;
    else if (s == "B")  return ATOM_T_BORON;
    else if (s == "I")  return ATOM_T_IODINE;


    // In case of emergency
    std::cerr << "Element |" << s << "| not recognized." << std::endl;
    throw "Atom type not recognized...attention required.";

    return ATOM_T_UNKNOWN;
}

SpecialEnumT AtomT::convertToSpecialEnum(const std::string& s)
{
    //
    // Force the contents to be lower case
    //
    std::string loweredS = s;
    std::transform(loweredS.begin(), loweredS.end(), loweredS.begin(), ::tolower);


    if (loweredS == "am")        return SPECIAL_T_AM;
    else if (loweredS == "ar")   return SPECIAL_T_AROMATIC;
    else if (loweredS == "pl")   return SPECIAL_T_PL;
    else if (loweredS == "co")   return SPECIAL_T_CO;
    else if (loweredS == "o")    return SPECIAL_T_O;
    else if (loweredS == "cat")  return SPECIAL_T_CAT;

    // In case of emergency
    if (loweredS.length() > 0)
    {
        std::cerr << "Special Enumeration |" << s << "| not recognized." << std::endl;
        throw "Atom type not recognized...attention required.";
    }

    return SPECIAL_T_NONE;
}

//
// Split the string into the constituent elements: <ELEMENT>.(<Specification> | <VAL>)
//
// returns <VAL>
//
int parseString(const std::string& s, std::string& prefix, std::string& suffix)
{
    suffix = "";
    
    // Look for the period
    int pos = s.find(".");

    // If no period, then we have just an element
    if (pos == std::string::npos)
    {
        prefix = s;
        return -1;
    }

    prefix = s.substr(0, pos);

    suffix = s.substr(pos+1);

    // Check if the last position is a numeric value
    std::string lastChar = s.substr(s.size() - 1);

    int val = -1;
    if (isdigit(lastChar[0]))
    {
        val = atoi(lastChar.c_str());
        // Cut down the suffix string by one character.
        suffix = suffix.substr(0, suffix.size() - 1);
    }

    return val;
}

//
// An AtomT is of the form: <Atom-Type>.<int-type> OR
//                          <Atom-Type>.<string-type>
AtomT::AtomT(const std::string& in)
{
    std::string prefix = "";
    std::string suffix = "";

    this->theAtomT.specificNum = parseString(in, prefix, suffix);

    this->theAtomT.atomType = convertToAtomEnum(prefix);
    this->theAtomT.specialT = convertToSpecialEnum(suffix);
}

// **************************************************************************************

bool AtomT::operator==(const AtomT& that) const
{
    if (this->theAtomT.atomType != that.theAtomT.atomType) return false;

    if (this->theAtomT.specificNum != that.theAtomT.specificNum) return false;

    if (this->theAtomT.specialT != that.theAtomT.specialT) return false;

    return true;
}


// **************************************************************************************

std::string AtomT::toString() const
{
    std::ostringstream oss;

    switch(theAtomT.atomType)
    {
      case ATOM_T_CARBON:      oss << "C"; break;
      case ATOM_T_CHLORINE:    return "Cl";
      case ATOM_T_HYDROGEN:    oss << "H"; break;
      case ATOM_T_NITROGEN:    oss << "N"; break;
      case ATOM_T_OXYGEN:      oss << "O"; break;
      case ATOM_T_PHOSPHORUS:  oss << "P"; break;
      case ATOM_T_SULFUR:      oss << "S"; break;
      case ATOM_T_FLUORINE:    oss << "F"; break;
      case ATOM_T_BROMINE:     oss << "Br"; break;
      case ATOM_T_BORON:       oss << "B"; break;
      case ATOM_T_IODINE:      oss << "I"; break;
    }

    oss << ".";

    switch(theAtomT.specialT)
    {
      case SPECIAL_T_AM:       oss << "am"; break;
      case SPECIAL_T_AROMATIC: oss << "ar"; break;
      case SPECIAL_T_PL:       oss << "pl"; break;
      case SPECIAL_T_CO:       oss << "co"; break;
      case SPECIAL_T_O:        oss << "O"; break;
      case SPECIAL_T_CAT:      oss << "cat"; break;
    }

    if (theAtomT.specificNum != -1) oss << theAtomT.specificNum;

    return oss.str();
}

// **************************************************************************************

std::string AtomT::getAtomType() const
{
    switch(theAtomT.atomType)
    {
      case ATOM_T_CARBON:     return "C";
      case ATOM_T_CHLORINE:   return "Cl";
      case ATOM_T_HYDROGEN:   return "H";
      case ATOM_T_NITROGEN:   return "N";
      case ATOM_T_OXYGEN:     return "O";
      case ATOM_T_PHOSPHORUS: return "P";
      case ATOM_T_SULFUR:     return "S";
      case ATOM_T_FLUORINE:   return "F";
      case ATOM_T_BROMINE:    return "Br";
      case ATOM_T_BORON:      return "B";
      case ATOM_T_IODINE:     return "I";
    }

    std::cerr << theAtomT.atomType << std::endl;

    return "Error";
}

// **************************************************************************************

std::ostream& operator<< (std::ostream& os, const AtomT& atomType)
{
    os << atomType.toString();

    return os;
}
