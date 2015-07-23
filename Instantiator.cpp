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

#include <vector>
#include <queue>
#include <iostream>
#include <memory>
#include <time.h>
#include <pthread.h>
#include <map>
#include <algorithm>


#include "Molecule.h"
#include "Rigid.h"
#include "Linker.h"


#include "MoleculeHashHypergraph.h"
#include "EdgeAggregator.h"


#include "Instantiator.h"
#include "OBWriter.h"
#include "Utilities.h"
#include "IdFactory.h"
#include "Constants.h"
#include "OBWriter.h"
#include "Options.h"
#include "bloom_filter.hpp"



// 0 indicates we let the queue size be limitless.
const unsigned Instantiator::MAX_QUEUE_SIZES[22] = { 0,   // Level 0
                                                     0,   //       1
                                                     300, //       2
                                                     10,  //       3
                                                     200, //       4
                                                     300, //       5
                                                     500, //       6
                                                     500, //       7
                                                     500, //       8
                                                     1000,//       9
                                                     1000,//       10
                                                     1000,//       11
                                                     1000,//       12
                                                     500, //       13
                                                     500, //       14
                                                     500, //       15
                                                     500, //       16
                                                     500, //       17
                                                     500, //       18
                                                     500, //       19
                                                     500, //       20
                                                     1    //       21
                                                   };

// The anticipated sizes of the level (at max). 0 indicates we are not using a Bloom filter.
const unsigned long long Instantiator::LEVEL_SIZES[22] = { 0,       // Level 0
                                                           0,       //       1
                                                           500,     //       2
                                                           10000,   //       3
                                                           300000,  //       4
                                                           1000000, //       5
                                                           5000000, //       6
                                                           15000000,//       7
                                                           30000000,//       8
                                                           30000000,//       9
                                                           30000000,//       10
                                                           15000000,//       11
                                                           5000000, //       12
                                                           2500000, //       13
                                                           1000000, //       14
                                                           500000,  //       15
                                                           100000,  //       16
                                                           50000,   //       17
                                                           25000,   //       18
                                                           10000,   //       19
                                                           5000,    //       20
                                                           1000     //       21
                                                         };



Instantiator::Instantiator(OBWriter*const obWriter, std::ostream& out) : writer(obWriter),
                                                                         ds(out),
                                                                         excluded(0),
                                                                         overallMoleculeCount(0)
{
    graph = new MoleculeHashHypergraph(HIERARCHICAL_LEVEL_BOUND + 1);

    // The hypergraph lock
    pthread_mutex_init(&graph_lock, NULL);

    // The threads and locks for the producer-consumer containers.
    level_queues = new std::queue<Molecule*>[HIERARCHICAL_LEVEL_BOUND + 1];
    moleculeLevelCount = new int[HIERARCHICAL_LEVEL_BOUND + 1];

    // Create the bloom filters

    if (Options::THREADED)
    {
        queue_locks = new pthread_mutex_t[HIERARCHICAL_LEVEL_BOUND + 1]; 
        queue_threads = new pthread_t[HIERARCHICAL_LEVEL_BOUND + 1];
        completed_level = new bool[HIERARCHICAL_LEVEL_BOUND + 1];
        arg_pointer = new Instantiator_ProcessLevel_Thread_Args[HIERARCHICAL_LEVEL_BOUND + 1];
    }
    for (int m = 1; m <= HIERARCHICAL_LEVEL_BOUND; m++)
    {
        if (Options::THREADED)
        {
            // Initialize those container locks.
            pthread_mutex_init(&queue_locks[m], NULL);
    
            // Initialize the fact that we have not computed this level. 
            completed_level[m] = false;

            // set up arg structs
            arg_pointer[m].m = m;
            arg_pointer[m].graph = graph;
            arg_pointer[m].this_pointer = this;
        }

        // We have create 0 molecules at this level, thus far.
        moleculeLevelCount[m] = 0;
    }

    InitOverallFilter();
    InitLevelFilters();
}

//
// Initialize the Bloom filter among all levels
//
void Instantiator::InitOverallFilter()
{
    bloom_parameters parameters;

    // How many elements roughly do we expect to insert?
    // Count the approximated level sizes
    unsigned long long approx_count = 0;
    for (int m = 0; m <= HIERARCHICAL_LEVEL_BOUND + 1; m++)
    {
        approx_count += LEVEL_SIZES[m];
    }
    parameters.projected_element_count = approx_count;

    // Maximum tolerable false positive probability? (0,1)
    parameters.false_positive_probability = 0.01; // 1%

    // Simple randomizer (optional)
    parameters.random_seed = 0xA5A5A5A5;

    if (!parameters)
    {
       std::cerr << "Error - Invalid set of bloom filter parameters!" << std::endl;
       return;
    }

    parameters.compute_optimal_parameters();

    // Create the Bloom filter.
    overall_filter = new bloom_filter(parameters);
}

//
// Initialize the Bloom filter at each level
//
void Instantiator::InitLevelFilters()
{
    bloom_parameters parameters;

    // Maximum tolerable false positive probability? (0,1)
    parameters.false_positive_probability = 0.001; // 1%

    // Simple randomizer (optional)
    parameters.random_seed = 0x5A5A5A5A;

    if (!parameters)
    {
       std::cerr << "Error - Invalid set of bloom filter parameters!" << std::endl;
    }

    //
    // Create the level filters
    //
    for (int m = 0; m <= HIERARCHICAL_LEVEL_BOUND + 1; m++)
    {
        if (LEVEL_SIZES[m] == 0)
        {
            filters.push_back(0);
        }
        else
        {
            // How many elements roughly do we expect to insert?
            parameters.projected_element_count = LEVEL_SIZES[m];

            parameters.compute_optimal_parameters();

            filters.push_back(new bloom_filter(parameters));
        }
    }
}



//
// Add the hyperedge to the hypergraph
//
void Instantiator::AddEdge(const std::vector<unsigned int>& antecedent,
                           unsigned int consequent,
                           EdgeAnnotationT* const annotation)
{
  /*
    pthread_mutex_lock(&graph_lock);

    graph->addEdge(antecedent, consequent, annotation);

    pthread_mutex_unlock(&graph_lock);
  */
}

//
// Add the hypernode to the hypergraph; success or failure is returned.
//
std::pair<unsigned int, bool> Instantiator::AddNode(MinimalMolecule* const mol, unsigned int sz)
{
    // We don't need to lock around the hypergraph since additions are level-based.
    // And each thread works on its own level.

    std::pair<int, bool> ret = graph->addNode(mol, sz);

    return ret;
}

//
// We first construct the base case of 2-Molecules.
// Then, we inductively start constructing 3-Molecules, 4-Molecules, etc.
//
#ifdef ZERO
MoleculeHashHypergraph* Instantiator::SerialInstantiate(std::vector<Linker*>& linkers,
                                                        std::vector<Rigid*>& rigids)
{
    //
    // Synthesizes level 2 molecules using SMI comparison.
    //
    InitializeSynthesis(linkers, rigids);

    // Indicate size of 1-M lists
    moleculeLevelCount[1] = baseMolecules.size();
    
    //
    // One level at a time:
    //    (a) Take the previous level molecules
    //    (b) Compose with the base molecules
    //    (c) Add those molecules to the next level queue
    //    (d) Kill the previous level's molecules   
    //
    for (int level = 2; level < HIERARCHICAL_LEVEL_BOUND; level++)
    {
        // Zip the smi file.
        //if (level == Options::SMI_LEVEL_BOUND + 1) writer->IndicateSMIwritingComplete();
        
        // Track the number of molecules at this level
        moleculeLevelCount[level] = level_queues[level].size();

        std::cerr << "Level " << level << " has " << moleculeLevelCount[level]
                  << " molecules to process." << std::endl; 

        int counter = 1;
        while(!level_queues[level].empty())
        {
            // Take a molecule from the in queue.
            Molecule* currentMol = level_queues[level].front();
            level_queues[level].pop();

            if (++counter % 500 == 0)
            {
                std::cerr << "Processing molecule " << counter
                          << " of " << moleculeLevelCount[level]
                          << " at level " << level << std::endl;
            }

            SynthesizeWithMolecule(currentMol, level);

            // Delete the current molecule; it has been processed completely.
            // Eliminate this code if we wish to kill an entire level, not molecule by molecule
            delete currentMol;
        }

        // Kill this level in the hypergraph
        graph->killLevel(level);

        // The Bloom Filter is no longer needed at this level.
        delete filters[level];
        filters[level] = 0;
    }

    std::cout << "Level\t" << "# Molecules" << std::endl;
    for (int m = 1; m <= HIERARCHICAL_LEVEL_BOUND; m++)
    {
       std::cout << m << "\t" << moleculeLevelCount[m] << std::endl;
    }

    // Tell the output engine we have completed synthesis.
    // This function then spins until the thread pool is complete.
    this->writer->IndicateSynthesisComplete();

    return graph;
}
#endif

//
// We first construct the base case of 2-Molecules.
// Then, we inductively start constructing 3-Molecules, 4-Molecules, etc.
//
MoleculeHashHypergraph* Instantiator::SerialInstantiate(std::vector<Linker*>& linkers,
                                                        std::vector<Rigid*>& rigids)
{
    //
    // Synthesizes level 2 molecules using SMI comparison.
    //
    InitializeSynthesis(linkers, rigids);

    // Indicate size of 1-M lists
    moleculeLevelCount[1] = baseMolecules.size();

    //
    // Using the level 2 molecules as a base case, process indicating non-completion.
    //
    unsigned molsProcessed = 0;
    while (!level_queues[2].empty())
    {
        SerialInstantiateHelper(2, molsProcessed);
    }

    //
    // Kill all levels
    //
    for (int m = 2; m <= HIERARCHICAL_LEVEL_BOUND; m++)
    {
        // Kill this level in the hypergraph
        graph->killLevel(m);

        // The Bloom Filter is no longer needed at this level.
        delete filters[m];
        filters[m] = 0;
    }

    std::cout << "Level\t" << "# Molecules" << std::endl;
    for (int m = 1; m <= HIERARCHICAL_LEVEL_BOUND; m++)
    {
       std::cout << m << "\t" << moleculeLevelCount[m] << std::endl;
    }

    // Tell the output engine we have completed synthesis.
    // This function then spins until the thread pool is complete.
    this->writer->IndicateSynthesisComplete();

    return graph;
}

//
// Given the current level, generate molecules in (level + 1) up to the capacity specified.
// When the capacity is exceeded, call this function recursively to process (level + 1)
// ...inductive completion. 
//
void Instantiator::SerialInstantiateHelper(int level, unsigned& processedMols)
{
    std::cerr << "Processing level " << level << std::endl;

    //
    // We max out at a specific level
    //
    if (level >= HIERARCHICAL_LEVEL_BOUND)
    {
        // Kill the contents of the queue
        while (!level_queues[level].empty())
        {
            Molecule* currentMol = level_queues[level].front();
            level_queues[level].pop();
            delete currentMol;
        }

        // Leave; no need to process.
        return;
    }

    //
    // Completely process all molecules in this level into level + 1
    //
    while (!level_queues[level].empty())
    {
        //
        // Adhere to capacities specified for each level
        //
        while (MAX_QUEUE_SIZES[level + 1] == 0 || level_queues[level + 1].size() < MAX_QUEUE_SIZES[level + 1])
        {
            //
            // Take a molecule from this level queue.
            //
            Molecule* currentMol = level_queues[level].front();
            level_queues[level].pop();

            moleculeLevelCount[level]++;

            if (++processedMols % 1000 == 0 || level <= 6)
            {
                std::cerr << "Processing molecule " << moleculeLevelCount[level]
                          << " at level " << level
                          << " queue contains (" << level_queues[level].size()
                          << "); Overall Processed Count: " 
                          << processedMols << std::endl;
            }

            // Dump the processed histogram of molecules
            if (processedMols % 1000000 == 0)
            {
                std::cerr << "Level\t" << "# Molecules" << std::endl;
                for (int m = 2; m <= HIERARCHICAL_LEVEL_BOUND; m++)
                {
                    std::cerr << m << "\t" << moleculeLevelCount[m] << std::endl;
                }
            }

            SynthesizeWithMolecule(currentMol, level);

            // Delete the current molecule; it has been processed completely.
            // Eliminate this code if we wish to kill an entire level, not molecule by molecule
            delete currentMol;

            // If nothing left to process at this level, quit and go to next levels above.
            if (level_queues[level].empty()) break;
        }

        //
        // Recursively process (level + 1)
        //
        SerialInstantiateHelper(level + 1, processedMols);
    }
}


//
// Creates the 2-molecules and initializes the fragments.
//
void Instantiator::InitializeSynthesis(std::vector<Linker*>& linkers,
                                       std::vector<Rigid*>& rigids)
{
    this->writer->IndicateSynthesisStarted();

    InitializeBaseMolecules(rigids, linkers, baseMolecules);

    // Add  all the base molecules to the hypergraph
    foreach_molecules(m_it, baseMolecules)
    {
        graph->addNode((*m_it)->ConstructMinimalMolecule(), 1);
    }

    //
    // Construct the set of 2-Molecules from the rigids and linkers.
    //
    for (int m1 = 0; m1 < baseMolecules.size(); m1++)
    {
        for (int m2 = m1; m2 < baseMolecules.size(); m2++)
        {
            std::vector<EdgeAggregator*>* newEdges =
                                          baseMolecules[m1]->Compose(*baseMolecules[m2]);

            HandleNewMolecules(level_queues[2], &queue_locks[2], filters[2], newEdges);
        }
    }

    std::cerr << "Done creating level 2" << std::endl;
}

//
// Add all new deduced clauses to the worklist if they have not been deduced before.
// If the given clause has been deduced before, update the hyperedges that were generated
// previously
//
// Forward Instantiation does not permit any cycles in the resultant graph.
//
void Instantiator::HandleNewMolecules(std::queue<Molecule*>& worklist,
                                      pthread_mutex_t* worklist_lock,
                                      bloom_filter* const levelFilter,
                                      std::vector<EdgeAggregator*>* newEdges)
{
    // Consider adding only if there are, in fact, new molecules
    if (newEdges->empty())
    {
        delete newEdges;
        return;
    }

    //
    // Since all molecules we have deduced are of the same size (using a level-based
    // construction), the size of the molecules are the same (equal num fragments)
    //
    unsigned level = (*newEdges->begin())->consequent->size();

    //
    // Add all molecules to the hypergraph
    //
    for (std::vector<EdgeAggregator*>:: const_iterator e_it = newEdges->begin();
         e_it != newEdges->end();
         e_it++)
    {
        // Did we generate this molecule previously? Or probability removal?
        bool killMolecule = false;

        // SMI for this molecule
        std::string smi = (*e_it)->consequent->ConstructSMI();

        // Add the consequent node to the graph directly.
        // std::pair<unsigned int, bool> addedResult = AddNode(minMol, level);

        // If we are validating the original molecule, check
        //if (VALIDATE) Validate();

        //
        // Check the memory-less dictionary for this level
        //
        static unsigned prob_excluded = 0;
        static unsigned overall_filtered = 0;
        if (levelFilter->contains(smi))
        {
            killMolecule = true;
        }
        //
        // Check the filter that applies to ALL molecules
        //
        else if (overall_filter->contains(smi))
        {
            killMolecule = true;

            if (++overall_filtered % 100 == 0)
            {
                std::cerr << "Overall filtered: " << overall_filtered <<  std::endl;
            }
        }
        //
        // Do we prune with probabilities?
        //
        else if (level >= Options::PROBABILITY_PRUNE_LEVEL_START)
        {
            if (Molecule::ProbabilisticExclusion((*e_it)->consequent))
            {
                killMolecule = true;

                if (++prob_excluded % 1000 == 0)
                {
                    std::cerr << "Probability excluding molecule: " << prob_excluded
                          << " (" << 100 * float(prob_excluded) / (overallMoleculeCount + prob_excluded)
                          << "\%)" << std::endl;
                }
            }
        }

        //
        // EXCLUDE
        //
        if (killMolecule)
        {
            delete (*e_it)->consequent;
        }

        //
        // INCLUDE
        //
        else
        {
            overallMoleculeCount++;

            // Add to the level bloom filters
            levelFilter->insert(smi);

            // Add to the overall bloom filter
            overall_filter->insert(smi);

            // Validation does not require output
            if (!VALIDATE) this->writer->OutputMoleculeAppendExternalSMI(smi);

            if (Options::THREADED) pthread_mutex_lock(worklist_lock);
            worklist.push((*e_it)->consequent);
            if (Options::THREADED) pthread_mutex_unlock(worklist_lock);
        }

        // Add the actual edge
        // AddEdge((*e_it)->antecedent, addedResult.first, (*e_it)->annotation);

        // We are done with this edge structure; delete it.
        delete (*e_it);
    }
    
    // Kill the edge list itself.
    delete newEdges;
}

//
// On the fly validation of molecules synthesized;
// Exits if the validation molecule was generated.
//
void Instantiator::Validate(const std::string& syn_smi) const
{
    // Convert
    if (syn_smi != validation_smi) return;

    std::cerr << "The give molecule has been synthesized: " << std::endl;
    std::cerr << "Validation: |" << validation_smi << "|" << std::endl; 
    std::cerr << "Synthesized: |" << syn_smi << "|" << std::endl;

    std::cerr << "Exiting..." << std::endl;
    exit(0);
}


//
// Initialize the linkers and rigids as required; the baseMolecules list will then be
// used as a reference container throughout synthesis.
//
void Instantiator::InitializeBaseMolecules(const std::vector<Rigid*>& rigids,
                                           const std::vector<Linker*>& linkers,
                                           std::vector<Molecule*>& baseMolecules)
{
    // Clear the list just in case.
    baseMolecules.clear();

    // Assign the linkers and rigids unique ids; these correspond EXACTLY to the indices of
    // the containers used for determing molecular (non)-isomorphism.
    foreach_rigids(r_it, rigids)
    {
        (*r_it)->setUniqueIndexID(moleculeIDFactory.getNextId());
        baseMolecules.push_back(*r_it);
    }

    foreach_linkers(l_it, linkers)
    {
        (*l_it)->setUniqueIndexID(moleculeIDFactory.getNextId());
        baseMolecules.push_back(*l_it);
    }

    // The set of base molecules is static in the synthesis process; therefore,
    // we set the (static) reference base set of molecules in the Molecule class
    // so the corresponding molecular fingerprint graph can be constructed and compared.
    Molecule::SetBaseMoleculeInfo(baseMolecules, rigids.size(), linkers.size());

    // Each molecule will contain a reference count of the number of each specific
    // linker / rigid in the particular molecule.
    foreach_molecules(m_it, baseMolecules)
    {
        (*m_it)->initFragmentDevices();
        (*m_it)->initGraphRepresentation();
    }
}

//
// Takes a single molecule and composes it with the base molecules to create the next level
// molecule.
//
void Instantiator::SynthesizeWithMolecule(const Molecule* const currentMol, int level)
{
    //
    // Compose with all of the base molecules
    //
    for (int m = 0; m < baseMolecules.size(); m++)
    {
        std::vector<EdgeAggregator*>* newEdges = currentMol->Compose(*baseMolecules[m]);

        //
        // Add the molecule to the next level queue; this depends on the level
        //
        HandleNewMolecules(level_queues[level + 1], 0, filters[level + 1], newEdges);
    }
}

//
//void Instantiator::ProcessLevel(std::vector<Molecule*>& baseMols,
//                                std::queue<Molecule*>& inSet,
//                                std::queue<Molecule*>& outSet,
//                                pthread_mutex_t& in_lock,
//                                pthread_mutex_t& out_lock,
//                                bool* previousLevelComplete,
//                                bool* thisLevelComplete)
//
void *ProcessLevel(void *ptr_void)
{
    //  unpacking arguments structure into mutiple local pointers
    Instantiator_ProcessLevel_Thread_Args * args = (Instantiator_ProcessLevel_Thread_Args *)ptr_void;
    int m = args->m; // level number
    Instantiator * This=(Instantiator *)args->this_pointer; // this pointer of calling class (Instantiator)

    //
    //  recast variables for local use (from the spawned thread record we were passed)
    //
    std::vector<Molecule*> *baseMols = &(This->baseMolecules);
    std::queue<Molecule*> *inSet = &(This->level_queues[m-1]);
    std::queue<Molecule*> *outSet = &(This->level_queues[m]);
    pthread_mutex_t *in_lock = &(This->queue_locks[m-1]);
    pthread_mutex_t *out_lock = &(This->queue_locks[m]);
    bool* previousLevelComplete = &(This->completed_level[m-1]);
    bool* thisLevelComplete = &(This->completed_level[m]);

    //
    // A structure for sleeping for 0.1 seconds
    //
    struct timespec sleepTime;
    struct timespec remTime;     // Remaining time
    sleepTime.tv_sec = 0;
    sleepTime.tv_nsec = 100000000L; // 0.1 seconds

    //
    // Keep consuming molecules as long as the previous level is incomplete or this
    // level queue contains molecules to process.
    //
    while (!(*previousLevelComplete) || !inSet->empty())
    {
        //
        // Nothing to process, currently, but the level is incomplete.
        //
        if (inSet->empty())
        {
            nanosleep(&sleepTime, &remTime);
        }
        else
        {
            //
            // If a greater level has elements in their queue, pause this thread for
            // a while.
            // Anything over level 13 should fly through.
            //
            bool process = false;

            if (m >= 13) process = true;
            else if (This->level_queues[m].size() < Instantiator::MAX_QUEUE_SIZES[m])
            {
                process = true;
            }
            //
            // Process a molecule in the queue
            //
            if (!process)
            {
                sleep(5);
            }
            else if (process)
            {
                //
                // Acquire a molecule to process.
                //
                pthread_mutex_lock(in_lock);
                Molecule* molToProcess = inSet->front();
                inSet->pop();
                pthread_mutex_unlock(in_lock);
                This->moleculeLevelCount[m-1]++;
                This->overallMoleculeCount++;

                if (This->overallMoleculeCount % 500 == 0 || m <= 6)
                {
                    std::cout << "Took molecule "
                              << This->moleculeLevelCount[m - 1]
                              << " off level " << m-1 << "; queue contains ("
                              << inSet->size() << "); Overall Count: "
                              << This->overallMoleculeCount << std::endl;
                }

                //
                // Process the molecule by composing it with all the base molecules.
                //
                int level = m - 1;
                for (int mol = 0; mol < Molecule::baseMolecules.size(); mol++)
                {
                    std::vector<EdgeAggregator*>* newEdges =
                                     molToProcess->Compose(*Molecule::baseMolecules[mol]);

                    //
                    // Add the molecule to the next level queue; this depends on the level
                    //
                    This->HandleNewMolecules(This->level_queues[level + 1],
                                             &This->queue_locks[level + 1],
                                             This->filters[level + 1],
                                             newEdges);
                }

                // We have successfully processed this molecule;
                // kill unneeded items in the molecule class.
                // Elements will persist in the MinimalMolecule representation
                // in the hypergraph.
                delete molToProcess;
            }
        }
    }

    // Indicate this level is complete.
    *thisLevelComplete = true;

    // Zip the smi file.
    //if (m == Options::SMI_LEVEL_BOUND + 1) This->writer->IndicateSMIwritingComplete();

    // We are done with this level so kill all references to it in the hypergraph.
    std::cerr << "Killing level " << (m - 1) << std::endl; 
    if (m > 2) args->graph->killLevel(m-1);

    std::cerr << "Level " << (m-1) << " created "
              << This->moleculeLevelCount[m-1] << " molecules." << std::endl; 

    std::cerr << "Level " << m << " complete." << std::endl; 
}


//
// Threaded construction of the hypergraph using a hierarchical list of threads and containers.
// We first construct the base case of 2-Molecules.
// Then, we inductively start constructing 3-Molecules, 4-Molecules, etc.
//
MoleculeHashHypergraph* Instantiator::ThreadedInstantiate(std::vector<Linker*>& linkers,
                                                          std::vector<Rigid*>& rigids)
{
    InitializeSynthesis(linkers, rigids);

    // 1-Molecules and 2-Molecules have been processed.
    completed_level[0] = true;
    completed_level[1] = true;
    completed_level[2] = true;

    // Indicate size of 1-M and 2-M lists
    moleculeLevelCount[1] = baseMolecules.size();

    //
    // For each level, start a thread and compose the elements with the base set of molecules.
    //
    for (int m = 3; m <= HIERARCHICAL_LEVEL_BOUND; m++)
    {
        if (~pthread_create(&queue_threads[m], NULL, ProcessLevel, (void*)&arg_pointer[m]))
            {if (g_debug_output) {std::cout << "Level " << m << " thread created" << std::endl;}}
        else
            {if (g_debug_output) {std::cout << "Level " << m << " creation failed" << std::endl;}}
    }
    for (int m = 3; m <= HIERARCHICAL_LEVEL_BOUND; m++)
    {
        (void) pthread_join(queue_threads[m], NULL);
        if (g_debug_output) std::cout << "Level " << m << " thread removed" << std::endl;
    }

    std::cout << "Level\t" << "# Molecules" << std::endl; 
    for (int m = 1; m <= HIERARCHICAL_LEVEL_BOUND; m++)
    {
       std::cout << m << "\t" << moleculeLevelCount[m] << std::endl; 
    }

    // Tell the output engine we have completed synthesis.
    // This function then spins until the thread pool is complete. 
    this->writer->IndicateSynthesisComplete();

    return graph;
}
