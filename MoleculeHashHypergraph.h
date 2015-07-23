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

#ifndef _MOLECULE_HASH_HYPERGRAPH_GUARD
#define _MOLECULE_HASH_HYPERGRAPH_GUARD 1


#include <pthread.h>


#include "LevelHashMap.h"


class MoleculeHashHypergraph
{
  private:
    LevelHashMap** level_maps;
    unsigned int numLevels;

    // disregards killed nodes counted
    unsigned int overall_sz;
    unsigned int num_killed;

    // Lock for the important hypergraph operations: these include only operations attributed
    // to structures interacting with the main structure.
    pthread_mutex_t graph_lock;

    // This vector links all levels together.
    std::vector<MinimalMolecule*> nodes;

  public:
    unsigned int currentSize() const { return overall_sz - num_killed; }
    unsigned int nonKilledSize() const { return overall_sz; }
    
    MoleculeHashHypergraph(unsigned int n) : overall_sz(0), num_killed(0)
    {
        numLevels = n;
        level_maps = new LevelHashMap*[numLevels];
		
        for (int i = 0; i < numLevels; i++)
	{
            level_maps[i] = new LevelHashMap();
        }
   
        pthread_mutex_init(&graph_lock, NULL);
    }
	
    //
    // Delete all hash
    //
    ~MoleculeHashHypergraph()
    {
        for (int i = 0; i < numLevels; i++)
        {
            // we may have already deleted a level; check before deletion.
            if (level_maps[i]) delete level_maps[i];
        }

        delete[] level_maps;
    }

    // Once we have completed processing of a particular, we can kill all memory
    // related to them. 
    void killLevel(unsigned int level)
    {
        num_killed += level_maps[level]->size();

        delete level_maps[level];
		
	level_maps[level] = 0;

        std::cerr << "Killed level " << level << std::endl;  
    }

    //
    // This function is provided ONLY for debugging purposes.
    // This should NOT be called at all to prevent duplicate isomorphism checks.
    //
    bool contains(MinimalMolecule* const that, unsigned int level) const
    {
        return level_maps[level]->contains(that);
    }
	
    // We assume a call to contains has already been made (as to limit the
    // number of graph isomorphisms)
    std::pair<unsigned int, bool> addNode(MinimalMolecule* const that, unsigned int level)
    {
        std::pair<unsigned int, bool> added = level_maps[level]->add(that);

        if (added.second)
        {
            overall_sz++;

            pthread_mutex_lock(&graph_lock);

            int index = nodes.size();
            nodes.push_back(that);
            that->uniqueIndexID = index;

            pthread_mutex_unlock(&graph_lock);

            return std::make_pair(index, true);
        }

        return added;
    }
};

#endif
