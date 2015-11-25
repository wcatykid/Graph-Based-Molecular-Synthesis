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

#include <cstring>
#include <vector>
#include <bitset>
#include <utility>
#include <iomanip>
#include <pthread.h>
#include <cmath>
#include <gsl/gsl_rng.h>


#include<openbabel/descriptor.h>
#include<openbabel/fingerprint.h>


#include "Molecule.h"
#include "Bond.h"
#include "Atom.h"
#include "obgen.h"
#include "Thread_Pool.h"
#include "Rigid.h"
#include "Linker.h"


#include "EdgeAggregator.h"
#include "EdgeAnnotation.h"
#include "IdFactory.h"
#include "Utilities.h"
#include "Constants.h"
#include "Options.h"
#include "OBWriter.h"
#include "MinimalMolecule.h"
#include "SmiMinimalMolecule.h"
#include "EdgeDatabase.h"


// global static lock for openbabel
pthread_mutex_t Molecule::openbabel_lock;

unsigned int Molecule::RIGID_INDEX_START = -1;
unsigned int Molecule::RIGID_INDEX_END = -1;
unsigned int Molecule::LINKER_INDEX_START = -1;
unsigned int Molecule::LINKER_INDEX_END = -1;
unsigned int Molecule::FRAGMENT_END_INDEX = -1;
unsigned int Molecule::NUM_UNIQUE_FRAGMENTS = -1;
EdgeDatabase Molecule::edges;

std::vector<Molecule*> Molecule::baseMolecules;
IdFactory Molecule::connectionIdMaker(100);
static const unsigned int NO_CONNECTION = -1;



Molecule::Molecule() : // obmol(0),
                       //type(COMPLEX),
                       fragmentCounter(0)
{
    init_openbabel_lock();
}

//
// Kill all items in the object except persistent information, which includes
// MinimalMolecule information:
//  (1) smi
//  (2) fragmentGraph
//
Molecule::~Molecule()
{
    // The OpenBabel molecule.
    // if (obmol) delete obmol;
    // obmol = 0;


    foreach_atoms(a_it, this->atoms)
    {
        delete (*a_it);
    }

/*
    foreach_bonds(b_it, this->bonds)
    {
        delete (*b_it);
    }
*/

    atoms.clear();
    bonds.clear();

    if (fragmentCounter) delete[] fragmentCounter;
    fragmentCounter = 0;
}

Molecule::Molecule(OpenBabel::OBMol* mol, const std::string& theSMI) : //, MoleculeT t) :
    uniqueIndexID(-1),
    //obmol(mol),
    //smi(theSMI),
    fingerprint(0),  
    //type(t),
    fragmentCounter(0)
{
    init_openbabel_lock();

    // Create the initial atom / bond data based on obmol.
    localizeOBMol(mol);

    // Acquire the OpenBabel Lipinski descriptor values
    openBabelPredictLipinski(mol);
}


void Molecule::init_openbabel_lock()
{
    //
    // initializing global openbabel lock (once)
    //
    static bool openbabel_lock_init = false;
    if (openbabel_lock_init)
    {
        pthread_mutex_init(&openbabel_lock, NULL);
        openbabel_lock_init = true;
    }
}

//
// Convert this Molecule to a simple representation defined by signature characteristics.
//
MinimalMolecule* Molecule::ConstructMinimalMolecule()
{
    return new MinimalMolecule(this->fingerprint,
                               this->fragmentCounter,
                               Molecule::NUM_UNIQUE_FRAGMENTS);
}

//
// Convert this Molecule to a simple representation defined by signature characteristics.
//
SmiMinimalMolecule* Molecule::ConstructSmiMinimalMolecule()
{
    std::string sdf;
    this->WriteToOpenBabelFormat(sdf);

    std::string smi;
    OBWriter::ConvertToSMI(sdf, smi);

    return new SmiMinimalMolecule(smi,
                                  this->fingerprint,
                                  this->fragmentCounter,
                                  Molecule::NUM_UNIQUE_FRAGMENTS);
}


//
// Convert this Molecule to a simple representation defined by signature characteristics.
//
std::string Molecule::ConstructSMI() const
{
    std::string sdf;
    this->WriteToOpenBabelFormat(sdf);

    std::string smi;
    OBWriter::ConvertToSMI(sdf, smi);

    return smi;
}


void Molecule::initFragmentDevices()
{
    initFragmentInfo();
//    calcFragmentInfo();


    // Indicate we are using this fragment
    fragmentCounter[uniqueIndexID] = 1;

    //
    // Create the connection identifiers for this linker / rigid
    //
    // Find the connections for this molecule and create ids for them.
    // These ids are unique to the linker and rigid.
    for (int a = 0; a < atoms.size(); a++)
    {
        unsigned id = NO_CONNECTION;

        if (atoms[a]->getMaxConnect() > 0)
        {
            id = connectionIdMaker.getNextId();

            // connectionIDs.push_back(id);
            atoms[a]->setConnectionID(id);

        }
    }
}

//
// Calculate the number of linkers / rigids (copies and unique)
//
void Molecule::initFragmentInfo()
{
    if (this->fragmentCounter == 0)
    {
        int sz = Molecule::NUM_UNIQUE_FRAGMENTS;

        // Create the reference count array
        fragmentCounter = new unsigned short int[sz];

        // init the counters to zero
        memset(fragmentCounter, 0, sz * sizeof(unsigned short int));
    }
}

//
// Fancy way to calculate the number of fragments in a molecule.
//
unsigned Molecule::size() const
{
    unsigned sz = 0;

    for (int r = 0; r < Molecule::NUM_UNIQUE_FRAGMENTS; r++)
    {
        sz += fragmentCounter[r];
    }

    return sz;
}


void Molecule::GetNumLinkersRigids(int& numLinkers, int& numUniqueLinkers,
                                   int& numRigids, int& numUniqueRigids) const
{
    numLinkers = 0;
    numUniqueLinkers = 0;
    numRigids = 0;
    numUniqueRigids = 0;

    for (int r = Molecule::RIGID_INDEX_START; r <= Molecule::RIGID_INDEX_END; r++)
    {
        if (fragmentCounter[r] != 0)
        {
            numUniqueRigids++;
            numRigids += fragmentCounter[r];
        }
    }

    for (int ell = Molecule::LINKER_INDEX_START; ell <= Molecule::LINKER_INDEX_END; ell++)
    {
        if (fragmentCounter[ell] != 0)
        {
            numUniqueLinkers++;
            numLinkers += fragmentCounter[ell];
        }
    }
}

void Molecule::initGraphRepresentation()
{
    // fingerprint = new SimpleFragmentGraph();
}

void Molecule::SetBaseMoleculeInfo(const std::vector<Molecule*> baseMols,
                                  unsigned int numRigids, unsigned int numLinkers)
{
    baseMolecules = baseMols;

    // Set the base molecule indices.
    Molecule::RIGID_INDEX_START = 0;
    Molecule::RIGID_INDEX_END = numRigids - 1;
    Molecule::LINKER_INDEX_START = numRigids;
    Molecule::LINKER_INDEX_END = numRigids + numLinkers - 1;
    Molecule::FRAGMENT_END_INDEX = Molecule::LINKER_INDEX_END;
    Molecule::NUM_UNIQUE_FRAGMENTS = numRigids + numLinkers;
}

void Molecule::openBabelPredictLipinski(OpenBabel::OBMol* obmol)
{
    pthread_mutex_lock(&Molecule::openbabel_lock);

    // calculate the molecular weight, H donors and acceptors and the plogp
    OpenBabel::OBDescriptor* pDesc1 = OpenBabel::OBDescriptor::FindType("HBD");
    OpenBabel::OBDescriptor* pDesc2 = OpenBabel::OBDescriptor::FindType("HBA1");
    OpenBabel::OBDescriptor* pDesc4 = OpenBabel::OBDescriptor::FindType("logP");

    if (!pDesc1) cerr << "HBD not found" << endl;
    if (!pDesc2) cerr << "HBA1 not found" << endl;
    if (!pDesc4) cerr << "logP not found" << endl;
    if (!pDesc1 || !pDesc2 || !pDesc4) return;

    MolWt = obmol->GetMolWt();
    HBD = pDesc1->Predict(obmol);
    HBA1 = pDesc2->Predict(obmol);
    logP = pDesc4->Predict(obmol);

    pthread_mutex_unlock(&Molecule::openbabel_lock);
}

//
// Near the end of the synthesis process, there is little benefit 
// to composing molecules if the two molecules will exceed the additive molecular weight.  
//
bool Molecule::willExceedAdditiveThresholds(const Molecule &mol1, const Molecule &mol2)
{
    // HBD 
    if (0.41189 + 0.4898 * (mol1.getHBD() + mol2.getHBD()) > HBD_UPPERBOUND) return true;

    // HBA1
    if (0.278 + 0.93778 * (mol1.getHBA1() + mol2.getHBA1()) > HBA1_UPPERBOUND) return true;

    // Molecular weight
    if (6.6746 + 0.95965 * (mol1.getMolWt() + mol2.getMolWt()) > MOLWT_UPPERBOUND) return true;

    return false;
}

void Molecule::estimateLipinski(const Molecule& mol1, const Molecule &mol2)
{
    double calc_MolWt = mol1.getMolWt() + mol2.getMolWt();
    double calc_HBD = mol1.getHBD() + mol2.getHBD();
    double calc_HBA1 = mol1.getHBA1() + mol2.getHBA1();
    double calc_logP = mol1.getlogP() + mol2.getlogP();

    this->MolWt = 6.6746 + 0.95965 * calc_MolWt;
    this->HBD = 0.41189 + 0.4898 * calc_HBD;
    this->HBA1 = 0.278 + 0.93778 * calc_HBA1;
    this->logP = 0.84121 + 0.59105 * calc_logP;
}

void Molecule::localizeOBMol(OpenBabel::OBMol* obmol)
{
    // Locking open babel since it is not thread-safe (at all)
    pthread_mutex_lock(&Molecule::openbabel_lock);

    int numOfAtoms = obmol->NumAtoms();
    int numOfBonds = obmol->NumBonds();

    //
    // Translate the OB atoms into our local atoms; create the space here and
    // then overwrite during parsing.
    //
/*
    for(int x = 0; x < numOfAtoms; x++)
    {
        this->atoms.push_back(Atom());
    }
*/

    //
    // Translate the OB Bonds into our local bonds.
    //
    for (int x = 0; x < numOfBonds; x++)
    {
        OpenBabel::OBBond* oneObBond = obmol->GetBondById(x);

        this->addBond((int)oneObBond->GetBeginAtom()->GetId(), 
                      (int)oneObBond->GetEndAtom()->GetId(),
                      oneObBond->GetBondOrder());
    }

    // Locking open babel since it is not thread-safe (at all)
    pthread_mutex_unlock(&Molecule::openbabel_lock);
}

//
// Molecular comparison is through the use of a local fingerprinting scheme.
// We construct a fingerprint by noting the connection anchors for fragments
// and constructing a graph based on those anchor points. A fingerprint equality
// check performs graph isomorphism.
//
bool Molecule::operator==(const Molecule& that) const
{
// std::cout << "Comparing: " << *this << " and " << that << std::endl;

    //
    // The fragment counter maintains the number of instances of each specific fragment;
    // if any of those counts differ, we have non-isomorphism.
    //
    for (int f = 0; f < Molecule::FRAGMENT_END_INDEX; f++)
    {
        if (this->fragmentCounter[f] != that.fragmentCounter[f])
        {
            return false;
        }
    }


    //
    // If we reach this point in the code, we can expect the two molecules to have the
    // same number of (1) linkers, (2) rigids, (3) unique rigids, (4) unique linkers,
    // (5) bonds, and (6) # atoms
    //
    //
    // Fingerprint verification.
    //
    // Fingerprint checking is last since it is slow; check other characteristics first.
    //
    return this->fingerprint->IsIsomorphicTo(that.getFingerprint());

    //
    // Compare molecules using the SMILES format: string comparison.
    //
    //return this->smi == that.smi;
}

// *****************************************************************************

//
// Exhaustively create any possible compositions between two Molecules
//
/*
    1. There should be no linker-to-linker connections.
       The only acceptable connections are rigid-to-rigid, and linker-to-rigid.

    2. The only loops that we should allow are loops that are inherent to the original linkers and rigids.
       In other words, if an original linker or rigid has a ring of atoms, then that loop is ok, but if you
       were to draw the molecule as a graph where each node is not an atom, but either a linker or a rigid,
       then that graph should not have any loops in it.

    3. The molecule should obey all 4 of Lipinski's rules, also known as the rule of 5.
       Here is a summary:
           a. A molecule should have a mass of less than 500 daltons, which can be calculated with
              the obmol function GetExactMass().
           b. No more than 5 hydrogen bond doners.
           c. No more than 10 hydrogen bond acceptors.
           d. An octanol-water partition coefficient log P not greater than 5.
*/

//
// (a) Check if the molecular weight is too heavy.
//
bool Molecule::exceedsMaxEstimatedThresholds()
{
    // We already checked if this molecule will exceed the upper bound for molecular weight.
    // No need to check it again.

    // (b) Hydrogen Bond donors
    if (HBD > HBD_UPPERBOUND)
    {
        return false;
    }

    // (c) Hydrogen Bond Acceptors
    if (HBA1 > HBA1_UPPERBOUND)
    {
        return false;
    }

    return MolWt > MOLWT_UPPERBOUND;
}

//
// Static function to check whether an OpenBabel Mol is Lipinski compliant.
//
// No need for locks since locks should go AROUND the function call.
bool Molecule::isOpenBabelLipinskiCompliant(OpenBabel::OBMol& mol)
{
    // calculate the molecular weight, H donors and acceptors and the plogp
    OpenBabel::OBDescriptor* pDesc1 = OpenBabel::OBDescriptor::FindType("HBD");
    OpenBabel::OBDescriptor* pDesc2 = OpenBabel::OBDescriptor::FindType("HBA1");
    OpenBabel::OBDescriptor* pDesc4 = OpenBabel::OBDescriptor::FindType("logP");

    if (!pDesc1) throw "HBD not found";
    if (!pDesc2) throw "HBA1 not found";
    if (!pDesc4) throw "logP not found";

    // (b) Hydrogen Bond donors
    if (pDesc1->Predict(&mol) > HBD_UPPERBOUND)
    {
        return false;
    }

    // (c) Hydrogen Bond Acceptors
    if (pDesc2->Predict(&mol) > HBA1_UPPERBOUND)
    {
        return false;
    }

    // Octanol-water partition coefficient log P not greater than 5
    if (pDesc4->Predict(&mol) > LOGP_UPPERBOUND)
    {
        return false;
    }

    return true;
}

bool Molecule::isLipinskiCompliant() const
{
    // (b) Hydrogen Bond donors
    if (HBD > HBD_UPPERBOUND) return false;

    // (c) Hydrogen Bond Acceptors
    if (HBA1 > HBA1_UPPERBOUND) return false;

    // Octanol-water partition coefficient log P not greater than 5
    if (logP > LOGP_UPPERBOUND) return false;

    return true;
}

/*
//
// Overall molecular criteria satisfaction.
//
bool Molecule::willSatisfiesMoleculeSynthesisCriteria()
{
    //
    // The only eliminating criteria is for molecule mass to be too large.
    //
    if (exceedsMaxEstimatedThresholds()) return false;

    // Do the linkers / rigids create a loop in the molecule?
    if (ContainsLoops()) return false;

    return true;
}
*/

std::vector<EdgeAggregator*>* Molecule::Compose(const Molecule& that) const
{
    std::vector<EdgeAggregator*>* newMolecules = new std::vector<EdgeAggregator*>();

    //
    // Pre-emptively check the Lipinski characteristics to see if there is a benefit
    // to composing these molecules; only perform this check if the user specified
    // Lipinski compliance.
    //
    if (Options::USE_LIPINSKI)
    {
        if (Molecule::willExceedAdditiveThresholds(*this, that)) return newMolecules;
    }

    //
    // Antecedent ids
    //
    std::vector<unsigned int> ante;
    ante.push_back(this->getUniqueIndexID());
    ante.push_back(that.getUniqueIndexID());

    //
    // For each atom in this molecule, does it connect to an atom in that molecule?
    //
    for (unsigned int thisA = 0; thisA < atoms.size(); thisA++)
    {
        for (unsigned int thatA = 0; thatA < that.atoms.size(); thatA++)
        {
            //
            // We've established the fact that these two particular atoms are connectable
            // Can we actually connect these two molecules at these two atoms?
            //
            if (atoms[thisA]->CanConnectTo(* that.atoms[thatA]))
            {

                if (g_debug_output)
                {
                    std::cerr << "Connection Possible: " << std::endl;
                    std::cerr << "\t" << atoms[thisA]->toString() << std::endl;
                    std::cerr << "\t" << that.atoms[thatA]->toString() << std::endl;
                }

                // Create a new molecule;
                // The indices are the new indices when the atoms and bonds are combined together.
                Molecule* newMol = 0;
                if (Options::OPENBABEL)
                {
                    // newMol = ComposeToNewOpenBabelMolecule(that, thisA + 1,
                    //                                 thatA + this->atoms.size() + 1);
                }
                else
                {
                    newMol = ComposeToNewLocalMolecule(that, thisA + 1,
                                                       thatA + this->atoms.size() + 1);
                }

                // Add the new molecule / edge to the list of new molecules
                newMolecules->push_back(new EdgeAggregator(ante, newMol, new EdgeAnnotationT()));
            }
        }
    } 

    return newMolecules;
}

// *****************************************************************************
//
// Create the local informations:
//    (1) Using the same style invoked by OpenBabel, we modify the indices by adding their
//        number of atoms in this to the index of the atoms in that.
//    (2) Bonds based on indices will be updated accordingly.
//

#ifdef ZERO

Molecule* Molecule::ComposeToNewOpenBabelMolecule(const Molecule& that,
                                                  int thisAtomIndex,
                                                  int thatAtomIndex) const
{
    static unsigned int num_blocked = 0;

    if (Options::THREADED)
    {
        num_blocked++;
        std::cerr << "Waiting on open babel lock: (" << num_blocked << ") are." << std::endl;

        pthread_mutex_lock(&Molecule::openbabel_lock);

        num_blocked--;
    }

    //
    // Combine the Open Babel representations.
    //
    OpenBabel::OBMol* newOBMol = new OpenBabel::OBMol(*this->obmol);

    *newOBMol += *that.obmol;

    // Add the new Open Babel bond.
    // This, along with the construction of the local molecule takes care of the new bond
    newOBMol->AddBond(thisAtomIndex, thatAtomIndex, 1); // Order of the bond is 1.

    // Remove the comment information as it is no longer relevant to this molecule.
    newOBMol->DeleteData("Comment");

    // Unlocking open babel; the next constructor call performs a relock
    std::string smi;
    OBWriter::ScrubAndConvertToSMIInternal(newOBMol, smi);

std::cerr << smi << std::endl;

    if (Options::THREADED)
    {
        pthread_mutex_unlock(&Molecule::openbabel_lock);
    }

    //
    //
    // Transfer the local data.
    //
    //
    // Create the new Molecule object; it will create the localized information.
    Molecule* newLocal = new Molecule(newOBMol, smi, COMPLEX);

    //
    // Copy the local atom information
    //
    int newAtomCount = 0;
    for (int a = 0; a < this->atoms.size(); a++, newAtomCount++)
    {
        newLocal->atoms[newAtomCount].SetBasedOn(this->atoms[a]);
    }

    int firstThatIndex = newAtomCount;
    for (int a = 0; a < that.atoms.size(); a++, newAtomCount++)
    {
        newLocal->atoms[newAtomCount].SetBasedOn(that.atoms[a]);
    }

    // Init the fragment counter container.
    newLocal->initFragmentInfo();

    // Combine all the linkers and rigids into this molecule.
    for (int f = 0; f <= FRAGMENT_END_INDEX; f++)
    {
        newLocal->fragmentCounter[f] = this->fragmentCounter[f] + that.fragmentCounter[f];
    }
    // Calculate all the fragment values: summary data.
    //newLocal->calcFragmentInfo();

    //
    // Add local information to the new molecule.
    // Bonds in open babel start indexing at 1.
    //
    newLocal->atoms[thisAtomIndex-1].addExternalConnection(); //thatAtomIndex-1);
    newLocal->atoms[thatAtomIndex-1].addExternalConnection(); //thisAtomIndex-1);

    //
    // Create the fingerprint graph for the new molecule by:
    //   (1) copying this fingerprint graph
/*
    newLocal->fingerprint = this->fingerprint->copy();

    // Add the new linker / rigid connection to the graph
    std::pair<unsigned int, unsigned int> toIndex;
    toIndex = newLocal->fingerprint->AddEdgeAndNode(
                          newLocal->atoms[thisAtomIndex - 1].getConnectionID(),
                          newLocal->atoms[thisAtomIndex - 1].getGraphNodeIndex(),
                          newLocal->atoms[thatAtomIndex - 1],
                          that);

std::cout << *this->fingerprint << std::endl << "+++++++++++" << std::endl;
std::cout << *that.fingerprint << std::endl << "===========" << std::endl;
std::cout << *newLocal->fingerprint << std::endl;

std::cout << "Graph node index: ("
          << newLocal->atoms[thisAtomIndex - 1].getGraphNodeIndex().first
          << ", " << newLocal->atoms[thisAtomIndex - 1].getGraphNodeIndex().second
          << ")" << std::endl;


    // Update the atoms of the new 'to' node to reflect the proper indices in the graph.
    for (int a = 0; a < that.atoms.size(); a++)
    {
        newLocal->atoms[firstThatIndex++].UpdateIndices(toIndex);
    }

*/

    // Estimate the Lipinski parameters.
    newLocal->estimateLipinski(*this, that);

std::string s;
newLocal->WriteToOpenBabelFormat(s);

    return newLocal;
}

#endif

// *****************************************************************************
//
// Create the local informations:
//    (1) Using the same style invoked by OpenBabel, we modify the indices by adding their
//        number of atoms in this to the index of the atoms in that.
//    (2) Bonds based on indices will be updated accordingly.
//
Molecule* Molecule::ComposeToNewLocalMolecule(const Molecule& that,
                                              int thisAtomIndex,
                                              int thatAtomIndex) const
{
    //
    //
    // Transfer the local data.
    //
    //
    // Create the new Molecule object; it will create the localized information.
    Molecule* newLocal = new Molecule();

    //
    // Copy the local atom information
    //
    foreach_atoms(a_it, this->atoms)
    {
        newLocal->atoms.push_back(Atom::ConstructAtom(**a_it));
    }

    foreach_atoms(a_it, that.atoms)
    {
        newLocal->atoms.push_back(Atom::ConstructAtom(**a_it));
    }

    //
    // Copy the local bond information
    //
    foreach_bonds(b_it, this->bonds)
    {
        newLocal->bonds.push_back(Bond(*b_it));
    }

    unsigned int offset = this->atoms.size();
    foreach_bonds(b_it, that.bonds)
    {
        newLocal->bonds.push_back(Bond(*b_it, offset));
    }
	
	// actual new bond (id, this-atom, that-atom, degree of bond)
    newLocal->bonds.push_back(Bond(thisAtomIndex - 1, thatAtomIndex - 1, 1));

    // Init the fragment counter container.
    newLocal->initFragmentInfo();

    // Combine all the linkers and rigids into this molecule.
    for (int f = 0; f <= FRAGMENT_END_INDEX; f++)
    {
        newLocal->fragmentCounter[f] = this->fragmentCounter[f] + that.fragmentCounter[f];
    }

    //
    // Add local information to the new molecule.
    // Bonds in open babel start indexing at 1.
    //
    newLocal->atoms[thisAtomIndex-1]->addExternalConnection(); // thatAtomIndex-1);
    newLocal->atoms[thatAtomIndex-1]->addExternalConnection(); // thisAtomIndex-1);

    // Create the fingerprint fragment graph for this new molecule.
    // unsigned short fromFragmentID = newLocal->atoms[thisAtomIndex - 1]->getOwnerFragment()->getUniqueIndexID();
    // unsigned short toFragmentID = newLocal->atoms[thatAtomIndex - 1]->getOwnerFragment()->getUniqueIndexID();
    // unsigned short fromConnID = newLocal->atoms[thisAtomIndex - 1]->getConnectionID();
    // unsigned short toConnID = newLocal->atoms[thatAtomIndex - 1]->getConnectionID();

    // short edgeID = edges.add(fromFragmentID, fromConnID, toConnID, toFragmentID);

    // newLocal->fingerprint = this->fingerprint->copyAndAppend(edgeID); 
    // std::cout << "Fingerprint: |" << *(newLocal->fingerprint) << "|" << std::endl;


    // Estimate the Lipinski parameters.
    newLocal->estimateLipinski(*this, that);

    return newLocal;
}


// *****************************************************************************

//
// On-Demand acquisition of the fingerprint one time.
//
SimpleFragmentGraph* Molecule::getFingerprint() const
{
    return this->fingerprint;
}

// *****************************************************************************

bool Molecule::addBond(int xID, int yID, unsigned int order) // , eTypeOfBondT bt, eStatusBitT s)
{
    this->bonds.push_back(Bond(/* this->bonds.size(), */ xID, yID, order));

    return true;
}

// *****************************************************************************

std::string Molecule::toString() const
{
    std::ostringstream oss;

    oss << "Molecule: " << uniqueIndexID << " ";

    if (IsLinker()) oss << " is a linker.";
    else if (IsRigid())  oss << " is a rigid.";

    else if (IsComplex())
    {
        int numLinkers;
        int numUniqueLinkers;
        int numRigids;
        int numUniqueRigids;

        GetNumLinkersRigids(numLinkers, numUniqueLinkers, numRigids, numUniqueRigids);

        oss << "There are ";
        oss << numRigids << " rigids, and ";
        oss << numLinkers;
        oss << " linkers." << std::endl;
    }

    oss << "There are ";
    oss << getNumberOfAtoms();
    oss << " atoms, and ";
    oss << getNumberOfBonds();
    oss << " bonds." << std::endl;

    if (getNumberOfAtoms() > 0)
    {
        oss << "Atoms:" << std::endl;

        foreach_atoms(a_it, this->atoms)
        {
            oss << "\t" << (*a_it)->toString() << std::endl;
        }
    }

    if (getNumberOfBonds() > 0)
    {
        oss << "Bonds:" << std::endl;

        foreach_bonds(b_it, this->bonds)
        {
            oss << "\t" << b_it->toString() << std::endl;
        }
    }

    return oss.str();
}


// *****************************************************************************

std::ostream& operator<< (std::ostream& os, const Molecule& mol)
{
    os << mol.toString() << std::endl;

    return os;
}

// *****************************************************************************

void Molecule::WriteToOpenBabelFormat(std::string& str) const
{
    std::ostringstream oss;

    //
    // Preamble
    //
    oss << "z_1.pdb" << std::endl;
    oss << " OpenBabel05191412593D" << std::endl; 

    oss << std::endl;

    //
    // Number of atoms, number of bonds
    //
    oss << std::setw(3) << this->atoms.size();
    oss << std::setw(3) << this->bonds.size();
    oss << std::setw(3) << 0;
    oss << std::setw(3) << 0;
    oss << std::setw(3) << 0;
    oss << std::setw(3) << 0;
    oss << std::setw(3) << 0;
    oss << std::setw(3) << 0;
    oss << std::setw(3) << 0;
    oss << std::setw(6) << "0999";
    oss << std::setw(6) << "V2000" << std::endl;

    foreach_atoms(a_it, this->atoms)
    {
        // Coordinates
        oss << std::setw(10) << "0.0000";
        oss << std::setw(10) << "0.0000";
        oss << std::setw(10) << "0.0000";

        oss << ' ';

        std::string atomtype = (*a_it)->getAtomType().getAtomType();
        oss << atomtype;
        if (atomtype.size() == 1) oss << ' ';

        for (int i = 1; i <= 12; i++)
        {
            oss << std::setw(3) << 0;
        }

        oss << std::endl;
    }

    foreach_bonds(b_it, this->bonds)
    {
        oss << std::setw(3) << b_it->getOriginAtomID() + 1;
        oss << std::setw(3) << b_it->getTargetAtomID() + 1;

        oss << std::setw(3) << b_it->getOrder();

        for (int i = 1; i <= 4; i++)
        {
            oss << std::setw(3) << 0;
        }

        oss << std::endl;
    }

    //
    // Post-amble
    //
    oss << "M  END";

// std::cerr << oss.str() << std::endl;

    // Assign for return
    str = oss.str();
}

// *****************************************************************************


//
// Probability-related code for inclusion / exclusion of a molecule
//
bool Molecule::ProbabilisticExclusion(const Molecule* const mol)
{
    static bool init_rng = false;
    static const gsl_rng_type* T;
    static gsl_rng* rec;
    if (!init_rng)
    {
        init_rng = true;

        gsl_rng_env_setup();

        T = gsl_rng_default;
        rec = gsl_rng_alloc (T);
    }

    int numLinkers;
    int numUniqueLinkers;
    int numRigids;
    int numUniqueRigids;

    mol->GetNumLinkersRigids(numLinkers, numUniqueLinkers, numRigids, numUniqueRigids);

    //
    // Acquire all of the probabilities associate with:
    //    (a) molecular weight
    //    (b) # rigid fragments
    //    (c) # linkers
    //    (d) log of ratio (linkers : rigids)
    //    (e) hydrogen binding donors
    //    (f) hydrogen binding acceptor 1
    //

    //    (a) molecular weight
    double mwProb = NormPdf(mol->getMolWt(), 428.366043, 91.124687);

    //    (b) # rigid fragments
    double numRigidProb = NormPdf(numRigids, 3.209722, 1.079512);

    //    (c) # linkers
    double numLinkerProb = LogisticPdf(numLinkers, 3.025175, 1.369960);

    //    (d) log of ratio (linkers : rigids)
    double log_ratio = log(((float)numLinkers) / ((float)numRigids));
    double ratioProb = LogisticPdf(log_ratio, -0.084292, 0.460030);

    //    (e) hydrogen binding donors
    double hbdProb = LogisticPdf(mol->getHBD(), 1.937285, 0.762586);

    //    (f) hydrogen binding acceptor 1
    double hbaProb = LogisticPdf(mol->getHBA1(), 6.056996, 1.312437);

    // Acquire the (cumulative) join probability distribution
    double cumProb = mwProb * numRigidProb * numLinkerProb * ratioProb * hbdProb * hbaProb;

    // Generate a random number between 0 and 1.
    double randJointProb = 1;
    for (int i = 0; i < 6; i++)
    {
        randJointProb *= gsl_rng_uniform(rec);
    }

// std::cerr << cumProb << " < " << randJointProb << " : " << (cumProb > randJointProb) << std::endl;

    // return false;

    return cumProb > randJointProb; 
    
    // gsl_rng_free(rec);
}

