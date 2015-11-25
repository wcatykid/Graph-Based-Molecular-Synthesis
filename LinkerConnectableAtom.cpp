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

#include <iostream>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>


#include "LinkerConnectableAtom.h"
#include "RigidConnectableAtom.h"


LinkerConnectableAtom::LinkerConnectableAtom(const LinkerConnectableAtom* const that)
{
    this->atomType = that->atomType;
    this->ownerFragment = that->ownerFragment;

    this->theAtom.connectionID = that->theAtom.connectionID;
    this->theAtom.maxConnect = that->theAtom.maxConnect;
    this->theAtom.numExternalConnections = that->theAtom.numExternalConnections;
    this->theAtom.numAllowConns = 1;
}

/**********************************************************************************/

LinkerConnectableAtom::LinkerConnectableAtom(int maxConn, const std::string& aType,
                                             const Linker* const owner)
{
    this->atomType = AtomT(aType);
    this->ownerFragment = (Molecule*)owner;

    this->theAtom.connectionID = 0;
    this->theAtom.maxConnect = maxConn;
    this->theAtom.numExternalConnections = 0;
    this->theAtom.numAllowConns = 1;
}

/**********************************************************************************/

LinkerConnectableAtom::~LinkerConnectableAtom()
{
}

/**********************************************************************************/

bool LinkerConnectableAtom::CanConnectTo(const Atom& that) const
{
    if (that.IsSimple()) return false;

    //
    // Disallow Linker-Linker connections.
    //
    if (that.IsLinkerAtom()) return false;

    //
    // Are there any allowable spots in the atoms to connect?
    //
    if(!this->SpaceToConnect()) return false;

    if(!that.SpaceToConnect()) return false;

    //
    // Does that atom allow the connection to this?
    //
    const RigidConnectableAtom& rAtom = static_cast<const RigidConnectableAtom&>(that);
    AtomT** thatAllowableTypes = rAtom.getAllowableTypes();
    unsigned numAllowableTypes = rAtom.getNumAllowableTypes();

    for (int t = 0; t < numAllowableTypes; t++)
    {
        if ((*thatAllowableTypes[t]) == this->atomType)
        {
            return true;
        }
    }

    return false;
}

/****************************************************************************************/

std::string LinkerConnectableAtom::toString() const
{
    std::ostringstream oss;

    oss << " Connections{ Max: " << theAtom.maxConnect;

    oss << " Allow: All Conn";

    oss << "\tNum ExtBonds: (" << theAtom.numExternalConnections << ")";

    oss << "  Conn Id: (" << theAtom.connectionID << ")";

    oss << " }";

    return oss.str();
}

/****************************************************************************************/

std::ostream& operator<< (std::ostream& os, ConnectableAtom& atom)
{
    os << atom.toString();

    return os;
}

/****************************************************************************************/
