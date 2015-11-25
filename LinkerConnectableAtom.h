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

#ifndef _LINKER_CONNECTABLE_ATOM_GUARD
#define _LINKER_CONNECTABLE_ATOM_GUARD 1


#include <string>
#include <iostream>


#include "AtomT.h"
#include "ConnectableAtom.h"


class Linker;


class LinkerConnectableAtom : public ConnectableAtom
{
  protected:

  public:
    bool CanConnectToAny() const { return true; }
    bool CanConnectTo(const Atom& that) const;

    bool IsLinkerAtom() const { return true; }
    bool IsRigidAtom() const { return false; }

    LinkerConnectableAtom(const LinkerConnectableAtom* const);
    LinkerConnectableAtom(int maxConn, const std::string&, const Linker* const owner);
    ~LinkerConnectableAtom();

    std::string toString() const;
    friend std::ostream& operator<< (std::ostream& os, LinkerConnectableAtom& atom);
};

#endif
