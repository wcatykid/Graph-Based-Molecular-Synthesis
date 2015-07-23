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
#include <fstream>


#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/fingerprint.h>


//#include "HyperGraph.h"
#include "Validator.h"
#include "Options.h"
#include "Utilities.h"
#include "Constants.h"

//
// Validate a single molecule.
//
bool Validator::Validate(OpenBabel::OBMol& validationMol)
{
    std::vector<unsigned int> validationFP;
    bool returnval = false;

    if (g_debug_output)
    {
        std::cerr << "Atoms: " << validationMol.NumAtoms() << std::endl;
        std::cerr << "Bonds: " << validationMol.NumBonds() << std::endl;
    }

    // Create the fingerprint for the validation molecule
    OpenBabel::OBFingerprint* fpType = OpenBabel::OBFingerprint::FindFingerprint("");

    // Acquire the fingerprint of the validation molecule so we can use it for
    // Tanimoto comparison.
    fpType->GetFingerprint(&validationMol, validationFP);

    if (g_debug_output)
    {
        std::cerr << "Validation: " << std::endl;
        foreach_uints(u_it, validationFP)
        {
            std::cerr << *u_it << "|";
        }
        std::cerr << std::endl;
    }

    //
    // Check this validation fingerprint against all of fingerprints in the list of
    // valid molecules.
    // 
    // Tracks the largest tanimoto coefficient of the synthesized molecule (and index?)
    //
    double maxTanimoto = -1;
    int maxIndex = -1;
    int molIndex = 0;
    foreach_obmol_points(m_it, this->molecules)
    {
        // Acquire the fingerprint of the hypergraph molecule.
        std::vector<unsigned int> hgFP;
        OpenBabel::OBFingerprint* fpType = OpenBabel::OBFingerprint::FindFingerprint("");
        fpType->GetFingerprint(*m_it, hgFP);

        // Debug output of the fingerprint.
        if (g_debug_output)
        {
            std::cerr << "Hypergraph: " << std::endl;
            foreach_uints(u_it, hgFP)
            {
                std::cerr << *u_it << "|";
            }
            std::cerr << std::endl;
        }

        // Acquire the tanimoto coefficient for the validation molecule and
        // the synthesized molecules in the hypergraph.
        double tanimoto = OpenBabel::OBFingerprint::Tanimoto(validationFP, hgFP);

        if (tanimoto > maxTanimoto)
        {
            maxTanimoto = tanimoto;
            maxIndex = molIndex;
        }

        if (g_debug_output) std::cerr << "Tanimoto: " << tanimoto << std::endl;

        if (tanimoto > (double)Options::TANIMOTO) 
        {
            returnval = true;
            break;
        }

        molIndex++;
    }


    std::ofstream logfile("Validation_logfile.txt", std::ofstream::out | std::ofstream::app); // append
    logfile << "Validation Molecule: " << validationMol.GetTitle() << " with ";
    logfile << "Synth Molecule: " << this->molecules[maxIndex]->GetTitle() << "\n";
    logfile << maxIndex << ": maxTanimoto = " << maxTanimoto;
    // if (returnval) logfile << " - Validated\n";
    // else logfile << " - Failed to Validate\n";
    logfile << std::endl;
    logfile.close();

    return returnval;
}


//
// Validate a list of molecules.
//
void Validator::Validate(std::vector<OpenBabel::OBMol>& molsToValidate)
{
    int counter = 1;

    foreach_obmols(m_it, molsToValidate)
    {
        if (!Validate(*m_it))
        {
            std::cerr << "Failed to validate: " << m_it->GetTitle() << std::endl;
        }
        else
        {
            std::cerr << "Validated molecule #" << counter
                      << ": " << m_it->GetTitle() << std::endl;
        }
        counter++;
    }
}


//
// For validation purposes, we use the original files in MOL format (before
// being stripped for linkers and rigids).
//
// (1) Parse the input files for all molecules to validate. 
// (2) Acquire all of the hypergraph molecules and perform obgen to re-acquire hydrogen bonds.
// (3) Compare the molecules using the Tanimoto similarity. (Get fingerprints of both
//     molecules and compare).
//
void Validator::Validate(const std::string& fileName)
{
    if (fileName == "")
    {
        std::cerr << "Validation file not specified; will not validate." << std::endl;
        return;
    }

    //
    // Input parser conversion functionality for Open babel
    //
    OpenBabel::OBConversion obConversion;
    obConversion.SetInFormat("MOL2");

    // The molecules to validate.
    std::vector<OpenBabel::OBMol> molsToValidate;
   
    std::cerr << "Reading Validation file " << fileName << std::endl;
 
    //
    // Read all of the OBMol objects using Open Babel; using their sample code style for reading.
    //
    OpenBabel::OBMol* mol = new OpenBabel::OBMol();
    bool notAtEnd = obConversion.ReadFile(mol, fileName.c_str());
    molsToValidate.push_back(*mol);

    while(notAtEnd)
    {
        // Create and parse using Open Babel
        OpenBabel::OBMol* mol = new OpenBabel::OBMol();
        notAtEnd = obConversion.Read(mol);
        if (notAtEnd) molsToValidate.push_back(*mol);
    }

    // Validate all molecules in the given file.
    Validate(molsToValidate);
}


