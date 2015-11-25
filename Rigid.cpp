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
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <cctype>


#include "Rigid.h"
#include "Utilities.h"
#include "RigidConnectableAtom.h"


Rigid::Rigid(OpenBabel::OBMol* obmol, const std::string& name) : Molecule(obmol, name)
{
    //
    // Acquire the comment data, make a copy, parse that comment.
    //
    OpenBabel::OBCommentData* comment =
                           static_cast<OpenBabel::OBCommentData*>(obmol->GetData("Comment"));

    std::string commentStr = comment->GetData();

    parseAppendix(commentStr, obmol->NumAtoms());

    // if (Options::OPENBABEL) OBWriter::ScrubAndConvertToSMIInternal(obmol, this->smi); 
}

void Rigid::parseAppendix(std::string& suffix, int numAtoms)
{
    // Use a string stream instead of manipulatiing the string
    std::stringstream suffStream(suffix);

    //
    // Read until we get "> <"
    //
    std::string line = "";
    while(line.find("> <") == std::string::npos)
    {
        getline(suffStream, line);
    }

    //
    // Read the Atom Types into a temporary
    //
    std::string atomType;
    std::vector<std::string> atomTypes;
    for(int x = 0; x < numAtoms; x++)
    {
        suffStream >> atomType;

std::cerr << atomType << std::endl;

        atomTypes.push_back(atomType);
    }

    // Get the next line.
    getline(suffStream, line);

    //
    // Read until we get "> <"
    //
    while(line.find("> <") == std::string::npos)
    {
        getline(suffStream, line);
    }

    //
    // Read Branches
    //

    // Parallels the atom arrays
    std::vector<std::string>* conns = new std::vector<std::string>[atomTypes.size()];    

    int atomId = -1;
    while (!isspace(suffStream.peek()))    
    {
        suffStream >> atomId;

        while (suffStream.peek() != '\n' && suffStream.peek() != '\r')
        {
            suffStream >> atomType;

            conns[atomId - 1].push_back(atomType);

            eatWhiteToNewLineOrChar(suffStream);
        }

        // Get the newline
        suffStream.get();
    }

    //
    // Read through the $$$$
    //
    while(line.find("$$$$") == std::string::npos)
    {
        getline(suffStream, line);
    }

    //
    // Actually create the atoms for this molecule.
    //
    for (int a = 0; a < numAtoms; a++)
    {
        if (conns[a].empty()) this->atoms.push_back(new Atom(atomTypes[a]));
        else
        {
            this->atoms.push_back(new RigidConnectableAtom(atomTypes[a], this, conns[a]));
            conns[a].clear();
        }
    }

    delete[] conns;
}
