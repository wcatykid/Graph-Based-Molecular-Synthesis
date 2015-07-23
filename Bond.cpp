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

#include <iostream>
//#include <cstring>
#include <sstream>
#include <string>
//#include <vector>
#include <cstdio>


#include "Bond.h"


/******************************************************************************************/

Bond::Bond(const Bond& that, unsigned offset)
{
    theBond.originAtomID = offset + that.theBond.originAtomID;
    theBond.targetAtomID = offset + that.theBond.targetAtomID;
    theBond.order = that.theBond.order;
}



/******************************************************************************************/

Bond::Bond(unsigned origin, unsigned target, unsigned ord)
{
    theBond.originAtomID = origin;
    theBond.targetAtomID = target;
    theBond.order = ord;
}

/*****************************************************************************************/

std::string Bond::toString() const
{
    std::ostringstream oss;

    oss << " from atom " << theBond.originAtomID
        << " to atom " << theBond.targetAtomID;

    return oss.str();
}


/*****************************************************************************************/

std::ostream& operator<< (std::ostream& os, Bond& bond)
{
    os << bond.toString();

    return os;
}

/*****************************************************************************************/
