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

#include <string>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <cctype>
#include <algorithm>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>
#include "errno.h"


#include "Utilities.h"



//
// log base 2 of the input value.
//
double log2(double value)
{
    if (value <= 0) return 0;

    return log10(value) / log10(2.0);
}

//
// How many bits in binary to represent this number?
//
unsigned int numBinaryBits(unsigned int value)
{
    return floor(log2(value)) + 1;
}

std::string MakeString(const char s[], int val)
{
    char buff[strlen(s) + 32];

    sprintf(buff, "%s: %d", s, val);
    std::string temp = buff;

    return temp;
}

std::string MakeString(const char s1[], const char s2[])
{
    char buff[strlen(s1) + strlen(s2)];

    sprintf(buff, "%s: %s", s1, s2);
    std::string temp = buff;

    return temp;
}

std::string MakeString(const char s1[], std::string s2)
{
    char buff[strlen(s1) + s2.size()];

    sprintf(buff, "%s: %s", s1, s2.c_str());
    std::string temp = buff;

    return temp;
}

void eatWhiteLines(std::istream& in)
{
    while (isspace(in.peek()))
    {
        while (in.get() != '\n')
        {
        }
    }
}

void eatWhiteToNewLineOrChar(std::istream& in)
{
    for (char c = in.peek(); c != '\n' && isalnum(c); c = in.peek())
    {
        in.get();
    }
}

void MakeBoolVector(vector<bool>& vec, int size)
{
    for (int i = 0; i < size; i++)
    {
        vec.push_back(false);
    }
}

bool ContainsFalse(const vector<bool>& vec)
{
    return std::find(vec.begin(), vec.end(), false) != vec.end();
}


float
NormPdf(float x, float loc, float scale){

  float norm_para, prob, pdf_val;

  norm_para = 1/(scale * sqrt(2 * M_PI));
  prob = exp (0.f - (x - loc) * (x - loc) / (2 * scale * scale));
  
  pdf_val = norm_para * prob;
  
  return pdf_val;
}


float
CauchyPdf(float x, float loc, float scale){
  float norm_para, prob, pdf_val;
  
  norm_para = 1/(M_PI * scale);
  prob = 1/(1 + ( (x-loc)/scale) * ((x-loc)/scale) );
  
  pdf_val = norm_para * prob;
  
  return pdf_val;
}



float
LogisticPdf(float x, float loc, float scale){
  float norm_para, e_power, prob, pdf_val;
  
  norm_para = 1/scale;
  e_power = exp( -(x - loc) / scale);
  prob = e_power / powf(1+e_power, 2.0);
  
  pdf_val = norm_para * prob;
  
  return pdf_val;
}


float
WaldPdf(float x, float loc, float scale){
  float norm_para, prob, pdf_val;
  
  float normed_x = (x - loc)/scale;

  norm_para = 1/(sqrt( 2 * M_PI * powf(normed_x, 3.0) ) * scale);
  prob = exp(-pow(normed_x-1, 2)/(2*normed_x));

  if (normed_x < 0) 
    pdf_val = 0.00000001f;
  else
    pdf_val = norm_para * prob;

  return pdf_val;
}


float 
LaplacePdf(float x, float loc, float scale){
  float normed_x, pdf_val;

  normed_x = fabs(x-loc) / scale;

  pdf_val = (1/(2 * scale)) * exp(- normed_x);
  
  return pdf_val;
}


bool DoesDirectoryExist(const std::string& theDir)
{
    struct stat s;
    int err = stat(theDir.c_str(), &s);

    if (err == -1)
    {
        if (ENOENT == errno)
        {
            // Does not exist
            return false;
        }
        else
        {
            std::cerr << "Error in Stat: Directory Existence." << std::endl;
            exit(1);
        }
    }
    else
    {
        if (S_ISDIR(s.st_mode))
        {
            // It's a dir
            return true;
        }
        else
        {
            // Exists but is not a directory 
            std::cerr << "Output directory " << theDir
                      << " exists, but is not a directory" << std::endl;
        }
    }

    return false;
}

void CleanDirectory(const std::string& theDir)
{
    //
    // Remove the contents of the directory
    //
    if (DoesDirectoryExist(theDir))
    {
        std::cout << "Cleaning output directory: " << theDir << std::endl;

        struct dirent *next_file;
        DIR *theFolder;

        char filepath[256];

        theFolder = opendir(theDir.c_str());

        while (next_file = readdir(theFolder))
        {
            // build the full path for each file in the folder
            sprintf(filepath, "%s/%s", theDir.c_str(), next_file->d_name);
            remove(filepath);
        }

        closedir(theFolder);
    }
}

void MakeDirectory(const std::string& theDir)
{
    if (!DoesDirectoryExist(theDir))
    {
        //
        // The directory does not exist;
        // we do need to create it.
        //
        mode_t mode;
        if (mkdir(theDir.c_str(), mode) == -1)
        {
            std::cout << "Creating output directory: " << theDir << std::endl;
        }
    }

    //
    // Change the permissions of the directory to be readable and writable.
    //
    char permission[] = "0777";
    int i = strtol(permission, 0, 8);

    if (chmod (theDir.c_str(), i) < 0)
    {
        fprintf(stderr, "Error in chmod(%s, %s)\n",
                theDir.c_str(), permission);
        exit(1);
    }
}
