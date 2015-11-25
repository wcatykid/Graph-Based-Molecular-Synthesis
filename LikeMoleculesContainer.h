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

#ifndef _LIKE_MOLECULES_CONTAINER_GUARD
#define _LIKE_MOLECULES_CONTAINER_GUARD 1

#include <cstring>

#include "MinimalMolecule.h"
 
class LikeMoleculesContainer
{
  private:
    std::vector<MinimalMolecule* > table;

  public:
    LikeMoleculesContainer() : table() { }
	
    //
    // Delete all entries
    //
    ~LikeMoleculesContainer()
    {
        for (int i = 0; i < table.size(); i++)
        {
            // std::cerr << "Deleting MinimalMol " << i << std::endl;
            delete table[i];
        }
        table.clear();
    }
	
    MinimalMolecule* contains(MinimalMolecule* const that) const
    {
        for (std::vector<MinimalMolecule*>::const_iterator it = table.begin();
	     it != table.end();
	     it++)
        {
	    if ((*it)->equals(that)) return *it;
        }

        return 0;
    }
	
    // We assume a call to contains has already been made (as to limit the
    // number of graph isomorphisms)
    void add(MinimalMolecule* const that)
    {
        table.push_back(that);
    }
 
    // Are all the molecules in this container defined by the given key value?
    bool definesKey(char* const thatKey)
    {
        return strcmp((*(table.begin()))->key, thatKey) == 0;
    }
};

#endif
