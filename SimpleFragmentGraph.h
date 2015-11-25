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

#ifndef _SIMPLE_FRAGMENT_GRAPH_GUARD
#define _SIMPLE_FRAGMENT_GRAPH_GUARD 1


#include <cstring>
#include <sstream>


//
// Wrapper around a simple list of integers.
// Those integers represent the edges that define the connections in a molecule.
//
class SimpleFragmentGraph
{
  public:
    SimpleFragmentGraph() : edgeIndices(0), sz(0) {}

    ~SimpleFragmentGraph()
    {
        if (edgeIndices) delete[] edgeIndices;
    }

    SimpleFragmentGraph* copyAndAppend(short edgeID) const
    {
        //
        // Create the new Simple representation
        //
        SimpleFragmentGraph* theCopy = new SimpleFragmentGraph();

        theCopy->edgeIndices = new short[this->sz + 1];

        //
        // Copy the old edges to the new edges; insert new edge into the ordered list
        //
        memcpy(theCopy->edgeIndices, this->edgeIndices, this->sz * sizeof(short)); 

        int index;
        for (index = sz; index > 0 && theCopy->edgeIndices[index - 1] > edgeID; index--)
        {
            theCopy->edgeIndices[index] = theCopy->edgeIndices[index - 1];
        }

        theCopy->edgeIndices[index] = edgeID;

        // Update the fact that we have 1 more edge in this list.
        theCopy->sz = this->sz + 1;

        return theCopy;
    }

    //
    // Graph Isomorphism using a simple comparison of the edges of the fragment graphs.
    //
    bool IsIsomorphicTo(SimpleFragmentGraph* const that) const
    {
        if (this->sz != that->sz) return false;

        short* first = this->edgeIndices;
        short* second = that->edgeIndices;

        int counter;
        for (counter = 0; counter < this->sz && *first == *second; counter++)
        {
            first++;
            second++;
        }

        return counter == this->sz;
    }


    std::string toString() const
    {
        std::ostringstream oss;

        oss << "(#" << sz << "): ";
        for (int index = 0; index < sz; index++)
        {
            oss << edgeIndices[index] << ' ';
        }

        return oss.str();
    }

    friend std::ostream& operator<< (std::ostream& os, const SimpleFragmentGraph& fg)
    {
        os << fg.toString();

        return os;
    }

  private:
    short* edgeIndices;
    unsigned short int sz;
};
	
#endif
