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

#ifndef _BOND_GUARD
#define _BOND_GUARD 1


class Bond
{
  private:
    typedef struct CompressedBondT
    {
        unsigned originAtomID : 15;
        unsigned targetAtomID : 15;
        unsigned order        : 2;
    } BondT;

    BondT theBond;

  public:
    //
    // Get Functions
    //
    int getOriginAtomID() const { return this->theBond.originAtomID; }
    int getTargetAtomID() const { return this->theBond.targetAtomID; }
    unsigned int getOrder() const { return this->theBond.order; }

    Bond(const Bond&, unsigned offset = 0);
    Bond(unsigned origin, unsigned target, unsigned order);
    ~Bond() {}

    std::string toString() const;
    friend std::ostream& operator<< (std::ostream& os, Bond& bond);

    bool operator==(const Bond& that) const
    {
        if (this->theBond.originAtomID != that.theBond.originAtomID) return false;

        if (this->theBond.targetAtomID != that.theBond.targetAtomID) return false;

        return true;
    }
};

#endif
