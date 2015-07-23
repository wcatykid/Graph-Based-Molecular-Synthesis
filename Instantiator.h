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

#ifndef _INSTANTIATOR_GUARD
#define _INSTANTIATOR_GUARD 1


#include <vector>
#include <map>
#include <queue>
#include <iostream>
#include <memory>
#include <pthread.h>


#include "Molecule.h"
#include "Rigid.h"
#include "Linker.h"
#include "MoleculeHashHypergraph.h"
#include "EdgeAnnotation.h"
#include "IdFactory.h"
#include "OBWriter.h"
#include "TimedHashMap.h"
#include "bloom_filter.hpp"



// threads require a struct to pass multiple arguments
struct Instantiator_ProcessLevel_Thread_Args
{
    int m; // level number
    MoleculeHashHypergraph* graph; // The single hypergraph used by all levels.
    void* this_pointer; // this pointer to calling class (Instantiator)
};

class Instantiator
{
  private:
    // Create necessary synthesis containers and init the linkers and rigids
    void InitializeBaseMolecules(const std::vector<Rigid*>& rigids,
                                 const std::vector<Linker*>& linkers,
                                 std::vector<Molecule*>& baseMolecules);

    void InitializeSynthesis(std::vector<Linker*>& linkers, std::vector<Rigid*>& rigids);

    // To generate unique molecular ids
    IdFactory moleculeIDFactory;

    // Contains all processed clauses and relationships amongst the clauses
    MoleculeHashHypergraph*  graph;

    // debug stream
    std::ostream& ds;

    void HandleNewMolecules(std::queue<Molecule*>& worklist,
                            pthread_mutex_t* wl_lock,
                            bloom_filter* const levelFilter,
                            std::vector<EdgeAggregator*>* newEdges);

    void SynthesizeWithMolecule(const Molecule* const currentMol, int level);
	
    void AddEdge(const std::vector<unsigned int>& antecedent,
                 unsigned int consequent,
                 EdgeAnnotationT* const annotation);

    std::pair<unsigned int, bool> AddNode(MinimalMolecule* const mol, unsigned int level);

    /*void ProcessLevel(std::vector<Molecule*>& baseMols,
                      std::queue<Molecule*>& inSet,
                      std::queue<Molecule*>& outSet,
                      pthread_mutex_t& in_lock,
                      pthread_mutex_t& out_lock,
                      bool* previousLevelComplete,
                      bool* thisLevelComplete);*/

    void InitOverallFilter();
    void InitLevelFilters();

    // Lock the hypergraph (for adding)
    pthread_mutex_t graph_lock;

    // A list of locks for all of the producer-consumer queues.
    pthread_mutex_t* queue_locks;

    // All of the hierarchical level threads.
    pthread_t* queue_threads;

    // Indicator that a level has completed processing.
    bool* completed_level;

    // The actual producer-consumer queue for each level.
    std::queue<Molecule*>* level_queues;

    // A bloom filter for each level beyond.
    std::vector<bloom_filter*> filters;

    // A bloom filter for each level beyond.
    bloom_filter* overall_filter;

    // array of args for each level thread
    Instantiator_ProcessLevel_Thread_Args *arg_pointer;

    // set of linkers and rigids (1-molecules)
    std::vector<Molecule*> baseMolecules;

    // Molecules per level (count) for debug
    int* moleculeLevelCount;
    unsigned long long overallMoleculeCount;

    // For output of molecules on the fly.
    OBWriter* const writer;

    // How many molecules were excluded using probabilistic techniques
    unsigned excluded;

    // The maximum number of molecules allowable in a queue.
    static const unsigned MAX_QUEUE_SIZES[22];

    // The expected number of molecules in a level, at maximum.
    static const unsigned long long LEVEL_SIZES[22]; 

    // On the fly validation of molecules synthesized;
    // Exits if the validation molecule was generated.
    void Validate(const std::string& syn_smi) const;

    // The molecule to validate
    std::string validation_smi;

  public:
    Instantiator(OBWriter*const obWriter, std::ostream& out = std::cout);

    ~Instantiator()
    {
        delete[] level_queues;
        delete[] moleculeLevelCount;
        delete graph;

        // Delete the Bloom filters.
        delete overall_filter;
        for (std::vector<bloom_filter*>::iterator it = filters.begin(); it != filters.end(); it++)
        {
            if (*it != 0) delete *it;
        }
        filters.clear();

        //
        // Threaded destruction
        //
        if (Options::THREADED)
        {
            delete[] queue_locks;
            delete[] queue_threads;
            delete[] completed_level;
            delete[] arg_pointer;
        }
    }

    // Main instantiation function for all linkers and rigidss; worklist technique to construct the graph
    MoleculeHashHypergraph* ThreadedInstantiate(std::vector<Linker*>& linkers,
                                                std::vector<Rigid*>& rigids);

    // Main instantiation function for all linkers and rigidss; worklist technique to construct the graph
    MoleculeHashHypergraph* SerialInstantiate(std::vector<Linker*>& linkers,
                                              std::vector<Rigid*>& rigids);

    // Recursive assistant for serial processing.
    void SerialInstantiateHelper(int level, unsigned& processedMols);

    unsigned getIncluded() const { return overallMoleculeCount; }
    unsigned getExcluded() const { return excluded; }

    // thread must be implemented as friend class
    friend void *ProcessLevel(void * args); // worker thread
};

#endif
