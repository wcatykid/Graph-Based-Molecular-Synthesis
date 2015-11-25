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

#ifndef _UTILITIES_GUARD
#define _UTILITIES_GUARD 1

#include <cmath>
#include <vector>
#include <iostream>
#include <string>


using namespace std;


// log base 2 of the input value.
double log2(double value);

// How many bits in binary to represent this number?
unsigned int numBinaryBits(unsigned int value);

string MakeString(const char[], int);
std::string MakeString(const char s1[], const char s2[]);
std::string MakeString(const char s1[], std::string s2);

template<class T>
bool Contains(std::vector<T> list, const T& val)
{
    for (int i = 0; i < list.size(); i++)
    {
        if (list[i].Equals(val)) return true;
    }
    return false;
}

//
// To simulate instanceof in Java.
//
template<typename CheckType, typename InstanceType>
bool instanceOf(InstanceType &Instance)
{
  return dynamic_cast<CheckType*>(Instance) != NULL;
}

void eatWhiteToNewLineOrChar(std::istream&);
void eatWhiteLines(std::istream&);


void MakeBoolVector(vector<bool>& vec, int size);
bool ContainsFalse(const vector<bool>& vec);

//
// Directory related functionality
//
void MakeDirectory(const std::string& theDir);
bool DoesDirectoryExist(const std::string& theDir);
void CleanDirectory(const std::string& theDir);

//
// Probability distributions
//
float NormPdf(float x, float loc, float scale);
float CauchyPdf(float x, float loc, float scale);
float LogisticPdf(float x, float loc, float scale);
float WaldPdf(float x, float loc, float scale);
float LaplacePdf(float x, float loc, float scale);


//
// Macros for simplifying code a bit.
//
#define foreach_rigids(it, expr) \
  for(std::vector<Rigid*>::const_iterator it = (expr).begin(); \
      it != (expr).end(); it++)

#define foreach_linkers(it, expr) \
  for(std::vector<Linker*>::const_iterator it = (expr).begin(); \
      it != (expr).end(); it++)

#define foreach_molecules(it, expr) \
  for(std::vector<Molecule*>::const_iterator it = (expr).begin(); \
      it != (expr).end(); it++)

#define foreach_bonds(it, expr) \
  for(std::vector<Bond>::const_iterator it = (expr).begin(); \
      it != (expr).end(); it++)

#define foreach_atoms(it, expr) \
  for(std::vector<Atom*>::const_iterator it = (expr).begin(); \
      it != (expr).end(); it++)

#define foreach_nodes(it, expr) \
  for(std::vector<FragmentGraphNode*>::const_iterator it = (expr).begin(); \
      it != (expr).end(); it++)

#define foreach_subnodes(it, expr) \
  for(std::vector<FragmentSubNode*>::const_iterator it = (expr).begin(); \
      it != (expr).end(); it++)

#define foreach_uints(it, expr) \
  for(std::vector<unsigned int>::const_iterator it = (expr).begin(); \
      it != (expr).end(); it++)

#define foreach_obmols(it, expr) \
  for(std::vector<OpenBabel::OBMol>::iterator it = (expr).begin(); \
      it != (expr).end(); it++)

#define foreach_obmol_points(it, expr) \
  for(std::vector<OpenBabel::OBMol*>::iterator it = (expr).begin(); \
      it != (expr).end(); it++)

#endif
