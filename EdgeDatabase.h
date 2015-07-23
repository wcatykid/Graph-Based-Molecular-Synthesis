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

#ifndef _EDGE_DATABASE_GUARD
#define _EDGE_DATABASE_GUARD 1

#include <utility>
#include <map>
#include <pthread.h>


#include "FragmentEdgeMap.h"
#include "Constants.h"

//
// A fragment has connection points.
// Considering a single fragment, we essentially have a hierarchical graph.
// Each unique fragment has unique connection ids.
// Connections are made between molecules via connection edges we store
// in this structure.
//
//
//        X ------- X                            X ------- X
//       /           \           edge           /           \
//      X   Frag Id   X(id:101) ------ (id:102)X   Frag ID   X
//       \     1     /                          \     2     /
//        X ------- X                            X ------- X
//
//  EdgeDatabase is keyed by the fragment ids.
//      MoleculeEdgeMap is keyed by 'from' in a (from, to) edge.
//          There are multiple lists in a fragments 
//
//

class EdgeDatabase
{
  private:
  
    // Size of the overall edge database AND
    // source of unique identifiers for the edges
    unsigned sz;
    pthread_mutex_t structure_lock;

  private:
    std::map<short, FragmentEdgeMap*> maps;
    
    //
    // Create the new container for the edge-list; add to the map
    //
    void addNewFragment(short fromFragmentID, short from, short to, short id)
    {
        FragmentEdgeMap* newTo = new FragmentEdgeMap();
        newTo->add(from, to, id);

        maps[fromFragmentID] = newTo;    
    }

    bool addEdge(short fromFragmentID, short from,
                 short to, short toFragmentID, unsigned id)
    {
        // Find the correct list in the map to search.
        std::map<short, FragmentEdgeMap*>::const_iterator vec_it = maps.find(fromFragmentID);
        
        // Vector (edge-list) not found; create and add.
        if (vec_it == maps.end())
        {
            addNewFragment(fromFragmentID, from, to, id);
        }
        else
        {
            // from-to; if duplicate, just indicate we did not add
            return vec_it->second->add(from, to, id);
        }
        
        return true;
    }

    // Returns the id of the edge
    short contains(short fromFragmentID, short from, short to) const
    {
        // Find the correct list in the map to search.
        std::map<short, FragmentEdgeMap*>::const_iterator it =  maps.find(fromFragmentID);
        
        // Vector (edge-list) not found
        if (it == maps.end()) return NOT_FOUND;

        // Otherwise, search the vector to see if the 'to' value is in the list.
        return it->second->contains(from, to);
    }
    
  public:
    EdgeDatabase() : sz(0)
    {
        pthread_mutex_init(&structure_lock, NULL);
    }
    
    //
    // Delete all entries
    //
    ~EdgeDatabase()
    {
        //  Traverse the map to remove all lists
        for (std::map<short, FragmentEdgeMap*>::iterator pair_it = maps.begin();
             pair_it != maps.end();
             pair_it++)
        {
            delete pair_it->second; 
        }
    }

    short add(short fromFragmentID, short from, short to, short toFragmentID)
    {
// std::cerr << "Attempting addition: (" << fromFragmentID << ", " << from << ") --- ("
//          << toFragmentID << ", " << to << ")" << std::endl;

        short ret_id = contains(fromFragmentID, from, to);
        if (ret_id != NOT_FOUND) return ret_id;
        
        pthread_mutex_lock(&structure_lock);

        // sz is the identifier generator
        short id = sz++;

        // Add both edges
        addEdge(fromFragmentID, from, to, toFragmentID, id);
        if (from != to) addEdge(toFragmentID, to, from, fromFragmentID, id);
        
        pthread_mutex_unlock(&structure_lock);

// std::cerr << "Added: (" << fromFragmentID << ", " << from << ") --- ("
//                         << toFragmentID << ", " << to << ")" << std::endl;

        return id;
    }

    std::string toString() const
    {
        std::ostringstream oss;

        //  Traverse the map to print all maps
        for (std::map<short, FragmentEdgeMap*>::const_iterator pair_it = maps.begin();
             pair_it != maps.end();
             pair_it++)
        {
            oss << pair_it->first << ":" << std::endl;
            oss << *(pair_it->second);
        }


        return oss.str();
    }

    friend std::ostream& operator<< (std::ostream& os, const EdgeDatabase& ed)
    {
        os << ed.toString();

        return os;
    }
};

#endif
