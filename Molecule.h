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

#ifndef _MOLECULE_GUARD
#define _MOLECULE_GUARD 1


#include <string>
#include <cstring> // for memset
#include <vector>
#include <memory>
#include <map>
#include <pthread.h>


#include <openbabel/mol.h>


#include "Bond.h"
#include "Atom.h"
#include "IdFactory.h"
#include "obgen.h"
#include "Constants.h"
#include "MinimalMolecule.h"
#include "SmiMinimalMolecule.h"
#include "EdgeDatabase.h"
#include "Utilities.h"


class EdgeAggregator;
class Rigid;
class Linker;
class SimpleFragmentGraph;

class Molecule
{
  public:
    Molecule();
    Molecule(OpenBabel::OBMol* mol, const std::string& theSMI);   //, MoleculeT t);
    ~Molecule();

    void setUniqueIndexID(unsigned int id) { uniqueIndexID = id; }
    unsigned int getUniqueIndexID() const { return uniqueIndexID; }

    //
    // Get Functions
    //
    virtual bool IsLinker() const { return false; }
    virtual bool IsComplex() const { return !this->IsLinker() && !this->IsRigid(); }
    virtual bool IsRigid() const { return false; }

    bool isLipinskiCompliant() const;
    double getMolWt() const { return MolWt; }
    double getHBD() const { return HBD; }
    double getHBA1() const { return HBA1; }
    double getlogP() const { return logP; }

    int getNumberOfAtoms() const { return this->atoms.size(); }
    int getNumberOfBonds() const { return this->bonds.size(); }

    // OpenBabel::OBMol* getOpenBabelMol() const { return obmol; }
    // void getSMI(std::string& s) const { s = smi; }
    // void setSMI(std::string that) { smi = that; }

    SimpleFragmentGraph* getFingerprint() const;

    bool addBond(int xID, int yID, unsigned int order); //, eTypeOfBondT bt, eStatusBitT s);

    std::string toString() const;
    friend std::ostream& operator<< (std::ostream& os, const Molecule& mol);
    virtual bool operator==(const Molecule& that) const;
    std::vector<EdgeAggregator*>* Compose(const Molecule&) const;

    // Acquire a summary of the linkers and rigids in this molecule.
    void GetNumLinkersRigids(int& numLinkers, int& numUniqueLinkers,
                             int& numRigids, int& numUniqueRigids) const;


    void openBabelPredictLipinski(OpenBabel::OBMol* obmol);
    static bool isOpenBabelLipinskiCompliant(OpenBabel::OBMol& mol);
    void estimateLipinski(const Molecule &mol1, const Molecule &mol2);
    static bool willExceedAdditiveThresholds(const Molecule &mol1, const Molecule &mol2);
    // Constructs a simple version of this molecule consisting of the fragment counts
    // and the fingerprint (fragment graph)
    MinimalMolecule* ConstructMinimalMolecule();
    SmiMinimalMolecule* ConstructSmiMinimalMolecule();

    std::string ConstructSMI() const;

    // The 'size' of a molecule is based on the number of total fragments.
    unsigned int size() const;

    // Initialize any containers to track fragments (linkers / rigids)
    void initFragmentDevices();

    // Calculate the number of linkers / rigids (copies and unique)
    //void calcFragmentInfo();

    // Initialize the fragment container
    void initFragmentInfo();

    // Initialize the graph-based representation of the fragment
    void initGraphRepresentation();

    // Collection of linkers and rigids for this synthesis.
    static std::vector<Molecule*> baseMolecules;
    static void SetBaseMoleculeInfo(const std::vector<Molecule*> baseMols,
                                    unsigned int numRigids, unsigned int numLinkers); 

    static unsigned int NUM_UNIQUE_FRAGMENTS;

    // Lock openbabel
    void init_openbabel_lock();
    static pthread_mutex_t openbabel_lock;

    void WriteToOpenBabelFormat(std::string&) const;

    static bool ProbabilisticExclusion(const Molecule* const);

  //
  /////////////////////////////////////////////////////////////////////////
  //
  protected:

    //
    // Instance Variables
    //

    // The unique identifier for this molecule
    unsigned int uniqueIndexID;

    // Local atoms and bonds
    std::vector<Atom*> atoms;
    std::vector<Bond> bonds;

    // Used for molecular comparison; the molecule represented as a graph
    SimpleFragmentGraph* fingerprint;

    // An array used to count the number of each specific linker /
    // rigid in this molecule
    unsigned short int* fragmentCounter;

    //
    // Lipinski Descriptors
    //
    double MolWt;
    double HBD;
    double HBA1;
    double logP;

    //
    // Statics
    //

    // Each linker / rigid has connection points;
    // we create unique ids for those connections.
    static IdFactory connectionIdMaker;

    //
    // Inline functions
    //
    virtual void parseAppendix(std::string& comment, int numAtoms = 0)
    {
        std::cerr << "Called Wrong parseAppendix::MOLECULE" << std::endl;
    }

  //
  /////////////////////////////////////////////////////////////////////////
  //
  private:
    void localizeOBMol(OpenBabel::OBMol* obmol);

    bool exceedsMaxEstimatedThresholds();
    bool ContainsLoops() const;
    bool satisfiesMoleculeSynthesisCriteria();
    /* Molecule* ComposeToNewOpenBabelMolecule(const Molecule& that,
                                            int thisAtomIndex,
                                            int thatAtomIndex) const;
    */
    Molecule* ComposeToNewLocalMolecule(const Molecule& that,
                                        int thisAtomIndex,
                                        int thatAtomIndex) const;


    static unsigned int RIGID_INDEX_START;
    static unsigned int RIGID_INDEX_END;
    static unsigned int LINKER_INDEX_START;
    static unsigned int LINKER_INDEX_END;
    static unsigned int FRAGMENT_END_INDEX;

    static EdgeDatabase edges;
};

#endif
