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

#ifndef _ATOM_GUARD
#define _ATOM_GUARD 1


#include <string>


#include "AtomT.h"


typedef unsigned MoleculeT;
typedef unsigned BooleanT;


const unsigned RIGID   = 0;
const unsigned LINKER  = 1;
const unsigned COMPLEX = 2;


class Molecule;


class Atom
{
  protected:
    // This atom's type
    AtomT atomType;

  public:
    const AtomT& getAtomType() const { return this->atomType; }

    virtual bool CanConnectToAny() const { return false; }
    virtual bool SpaceToConnect() const { return false; }
    virtual bool CanConnectTo(const Atom& that) const { return false; }
    virtual int getConnectionID() const { return -1; }
    virtual void setConnectionID(unsigned) { throw "Should not be called."; }
    virtual int getMaxConnect() const { return 0; }
    virtual void addExternalConnection() { throw "Should not be called."; }

    static Atom* ConstructAtom(const Atom& atom);

    virtual bool IsSimple() const { return true; }
    virtual bool IsConnectable() const { return false; }

    virtual bool IsLinkerAtom() const { return false; }
    virtual bool IsRigidAtom() const { return false; }

    virtual const Molecule* getOwnerFragment() const { throw "Should not be called."; }

    Atom(const Atom&);
    Atom() {}
    Atom(const AtomT& aType);
    Atom(const std::string& aType);
    virtual ~Atom();

    virtual std::string toString() const;
    friend std::ostream& operator<< (std::ostream& os, Atom& atom);
    bool operator==(const Atom& that) const;
};

#endif
