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

#ifndef _EDGE_AGGREGATOR_GUARD
#define _EDGE_AGGREGATOR_GUARD 1

#include <vector>

#include "Molecule.h"
#include "EdgeAnnotation.h"

//
// An aggregation class of information to pass back from instantiation.
//
class EdgeAggregator
{
  public:
    std::vector<unsigned int> antecedent;
    Molecule* consequent;
    EdgeAnnotationT* annotation;
        
    EdgeAggregator(const std::vector<unsigned int>& ante, Molecule* c, EdgeAnnotationT* ann)
                  : antecedent(ante), consequent(c), annotation(ann)
    {
    }

    ~EdgeAggregator()
    {
        // This annotation should persist into the hypergraph.
        // For now we delete; it is not-needed.
        delete annotation;
    }
};

#endif
