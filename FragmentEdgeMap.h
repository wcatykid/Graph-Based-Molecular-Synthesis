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

#ifndef _FRAGMENT_EDGE_MAP_GUARD
#define _FRAGMENT_EDGE_MAP_GUARD 1
 

#include <map>
#include <algorithm>
#include <sstream>


#include "FixedSortedList.h"
#include "Constants.h"

 
class FragmentEdgeMap
{
  private:
    std::map<short, FixedSortedList*> edgeMap;
    static unsigned LIST_SIZE;
    
    //
    // Create the new container for the edge-list; add to the map
    //
    void addNewEdgeList(short from, short to, short id)
    {
        FixedSortedList* newTo = new FixedSortedList(LIST_SIZE);
        newTo->add(to, id);

        edgeMap[from] = newTo;    
    }

    bool addEdge(short from, short to, unsigned id)
    {
        // Find the correct list in the map to search.
        std::map<short, FixedSortedList*>::const_iterator vec_it = edgeMap.find(from);
        
        // Vector (edge-list) not found; create and add.
        if (vec_it == edgeMap.end())
        {
            addNewEdgeList(from, to, id);
        }
        else
        {
            // from-to; if duplicate, just indicate we did not add
            return vec_it->second->add(to, id);
        }
        
        return true;
    }

  public:
    FragmentEdgeMap() : edgeMap() {}
    
    //
    // Delete all entries
    //
    ~FragmentEdgeMap()
    {
        //  Traverse the map to remove all lists
        for (std::map<short, FixedSortedList*>::const_iterator pair_it = edgeMap.begin();
             pair_it != edgeMap.end();
             pair_it++)
        {
            delete pair_it->second; 
        }
    }

    short add(short from, short to, int id)
    {
        short ret_id = contains(from, to);
        if (ret_id != NOT_FOUND) return ret_id;
        
        addEdge(from, to, id);
        
        return id;
    }

    // Returns the id of the edge
    short contains(short from, short to) const
    {
        if (edgeMap.empty()) return NOT_FOUND;

        // Find the correct list in the map to search.
        std::map<short, FixedSortedList*>::const_iterator it =  edgeMap.find(from);

        // Vector (edge-list) not found
        if (it == edgeMap.end()) return NOT_FOUND;

        // Otherwise, search the vector to see if the 'to' value is in the list.
        return it->second->contains(to);
    }

    std::string toString() const
    {
        std::ostringstream oss;

        //  Traverse the map to print all maps
        for (std::map<short, FixedSortedList*>::const_iterator pair_it = edgeMap.begin();
             pair_it != edgeMap.end();
             pair_it++)
        {
            oss << "\t" << pair_it->first << ": ";
            oss << *(pair_it->second) << std::endl;
        }


        return oss.str();
    }

    friend std::ostream& operator<< (std::ostream& os, const FragmentEdgeMap& fem)
    {
        os << fem.toString();

        return os;
    }
};

#endif
