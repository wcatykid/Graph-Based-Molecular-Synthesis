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

#ifndef _LEVEL_HASH_MAP_GUARD
#define _LEVEL_HASH_MAP_GUARD 1

#include <utility>
#include <map>


#include "LikeMoleculesContainer.h"


class LevelHashMap
{
  private:
    LikeMoleculesContainer** like_molecules;
    static const unsigned int SIZE = 10000;
    unsigned int sz;

    // Trusted hash function for strings (found online)
    unsigned int hash(char* str) const
    {
        unsigned long hash = 5381;
        int c;

        while (c = *str++)
        {
            hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
        }

        return hash % SIZE;
    }

    //
    // Handle collisions with linear probing
    //
    std::pair<unsigned int, bool> privateAdd(MinimalMolecule* const that)
    {
        unsigned int hashVal = hash(that->key);

// std::cerr << "Hash: " << hashVal << std::endl;

        // Add to the table directly or seek the next list.
        for (int index = 0; index < SIZE; index++)
        {
            //
            // Does this list exist?
            // If not, create the list, add the molecule, and indicate we added successfully.
            //
            unsigned int properIndex = (index + hashVal) % SIZE;
            if (like_molecules[properIndex] == 0)
            {
                like_molecules[properIndex] = new LikeMoleculesContainer();
                like_molecules[properIndex]->add(that);


                return std::make_pair(-1, true);
            }
            // The list exists, inquire if the key matches our key.
            if (like_molecules[properIndex]->definesKey(that->key))
            {
                //
                // Determine if this molecule is in the list; if it is, don't add.
                //
                MinimalMolecule* mol = like_molecules[properIndex]->contains(that);
                if (mol) return std::make_pair(mol->uniqueIndexID, false);

                // Add to the list with success.
                like_molecules[properIndex]->add(that);
                return std::make_pair(-1, true);
            }
        }

        return std::make_pair(-1, false);
    }

	
  public:
    LevelHashMap() : sz(0)
    {
        like_molecules = new LikeMoleculesContainer*[SIZE];
		
        for (int i = 0; i < SIZE; i++)
        {
            like_molecules[i] = 0;
        }
    }
	
    //
    // Delete all entries
    //
    ~LevelHashMap()
    {
        for (int i = 0; i < SIZE; i++)
        {
            if (like_molecules[i]) delete like_molecules[i];
        }

        delete[] like_molecules;
    }

    unsigned int size() const { return sz; }

    std::pair<unsigned int, bool> add(MinimalMolecule* const that)
    {
        std::pair<unsigned int, bool> added = privateAdd(that);

        if (added.second) sz++;

        return added;
    }


    //
    // This function is provided ONLY for debugging purposes.
    // This should NOT be called at all to prevent duplicate isomorphism checks.
    //
    bool contains(MinimalMolecule* const that) const
    {
        unsigned int hashVal = hash(that->key);
		
        // Probe linearly.
        for (int index = 0; index < SIZE; index++)
        {
            //
            // Does this list exist?
            // If not, we don't have containment.
            //
            unsigned int properIndex = (index + hashVal) % SIZE;

            if (like_molecules[properIndex] == 0) return false;
			
            // The list exists, inquire if the key matches our key.
            if (like_molecules[properIndex]->definesKey(that->key))
            {
                return like_molecules[properIndex]->contains(that);
            }
        }

        return false;
    }
};

#endif
