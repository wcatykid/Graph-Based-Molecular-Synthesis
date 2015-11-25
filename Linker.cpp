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

#include <vector>
#include <iostream>
#include <string>
#include <map>
#include <sstream>


#include <openbabel/mol.h>
#include <openbabel/bond.h>


#include "Linker.h"
#include "LinkerConnectableAtom.h"


Linker::Linker(OpenBabel::OBMol* obmol, const std::string& name) : Molecule(obmol, name)
{
    //
    // Acquire the comment data, make a copy, parse that comment.
    //
    OpenBabel::OBCommentData* comment = static_cast<OpenBabel::OBCommentData*>(obmol->GetData("Comment"));

    std::string commentStr = comment->GetData();

    parseAppendix(commentStr, obmol->NumAtoms());

    // if (Options::OPENBABEL) OBWriter::ScrubAndConvertToSMIInternal(obmol, this->smi); 
}

//
// Parse suffix to add max connection for each atom.
//
void Linker::parseAppendix(std::string& suffix, int numAtoms)
{
    // std::cerr << "Linker::parseAppendix: " << suffix << "|" << std::endl;

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
    // Now, read the MAX Connections for each atom.
    //
    int maxConnections = -1;
    std::string atomType;

    for(int x = 0; x < numAtoms; x++)
    {
        suffStream >> maxConnections;
        suffStream >> atomType;

        // A linker can link to any atom.
        this->atoms.push_back(new LinkerConnectableAtom(maxConnections, atomType, this));
    }
}
