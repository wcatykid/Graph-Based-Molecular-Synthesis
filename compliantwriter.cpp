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

#include <fstream>
#include <string>


#include <unistd.h>


#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/descriptor.h>


static const double MOLWT_UPPERBOUND = 570;
static const double HBD_UPPERBOUND = 5;
static const double HBA1_UPPERBOUND = 10;
static const double LOGP_UPPERBOUND = 7.2;

static char* outputFile;

bool isOpenBabelLipinskiCompliant(OpenBabel::OBMol* const mol)
{
    // calculate the molecular weight, H donors and acceptors and the plogp
    OpenBabel::OBDescriptor* pDesc1 = OpenBabel::OBDescriptor::FindType("HBD");
    OpenBabel::OBDescriptor* pDesc2 = OpenBabel::OBDescriptor::FindType("HBA1");
    OpenBabel::OBDescriptor* pDesc4 = OpenBabel::OBDescriptor::FindType("logP");

    if (!pDesc1) throw "HBD not found";
    if (!pDesc2) throw "HBA1 not found";
    if (!pDesc4) throw "logP not found";

    // (b) Hydrogen Bond donors
    if (pDesc1->Predict(mol) > HBD_UPPERBOUND)
    {
        // std::cerr << "Failed HBD" << std::endl;
        return false;
    }

    // (c) Hydrogen Bond Acceptors
    if (pDesc2->Predict(mol) > HBA1_UPPERBOUND)
    {
        // std::cerr << "Failed HBA1" << std::endl;
        return false;
    }

    // Octanol-water partition coefficient log P not greater than 5
    if (pDesc4->Predict(mol) > LOGP_UPPERBOUND)
    {
        // std::cerr << "Failed PLOGP" << std::endl;
        return false;
    }

    return true;
}

void Output(std::string smi)
{
    // Open the file for writing.
    std::ofstream out;

    out.open(outputFile, std::ofstream::out | std::ofstream::app); // append

    // Loop until we get access to the output file.
    while (out.fail())
    {
        sleep(1);

        out.open(outputFile, std::ofstream::out | std::ofstream::app); // append
    }

    // Write the molecule to the file.
    out << smi << std::endl;

    out.close();
}


void OutputMolecule(OpenBabel::OBMol* const mol)
{
    // (3) Check Lipinski Compliance
    // The molecule must be Lipinski compliant (using Open Babel)
    if (!isOpenBabelLipinskiCompliant(mol))
    {
        std::cerr << "Not OpenBabel Lipinski compliant." << std::endl;
        return;
    }

    // Create the SMI; output
    OpenBabel::OBConversion SMI_conv(&std::cin, &std::cout);

    // set conversion type(s) and verify it worked
    if(!SMI_conv.SetInAndOutFormats("SMI","SMI")) throw "SetInAndOutFormats failed!";

    // Pre-emptive extra run of OBGen, seems to stop segmentation fault
    //OBGen::fast_obgen(mol);

    std::string smi = SMI_conv.WriteString(mol); // convert and write to string

    // remove garbage and acquire only the SMI
    smi = smi.substr(0, smi.find('\t'));

    Output(smi);
}

int main(int argc,char **argv)
{
    if (argc != 3)
    {
        std::cerr << "Usage: compliantwriter <SDF-FILE>.sdf <output-file>" << std::endl; 
        return 1;
    }

    char* program_name = argv[0];
    std::string infilename = argv[1];
    outputFile = argv[2];

    // SMI reader
    OpenBabel::OBConversion ob_conv(&std::cin, &std::cout);

    // Check if the input format is SDF or SMI based on the file extension 
    if (infilename.rfind(".smi") != std::string::npos)
    {
        // set conversion type(s) and verify it worked
        if(!ob_conv.SetInAndOutFormats("SMI","SMI"))
        {
            std::cerr << "SetInAndOutFormats failed!" << std::endl;
            return 1;
        }
    }
    else if (infilename.rfind(".sdf") != std::string::npos)
    {
        // set conversion type(s) and verify it worked
        if(!ob_conv.SetInAndOutFormats("SDF","SMI"))
        {
            std::cerr << "SetInAndOutFormats failed!" << std::endl;
            return 1;
        }
    }
    else
    {
        std::cerr << ".sdf or .smi extension for input file not found." << std::endl;
        return 1;
    }

    std::ifstream ifs;

    // Read the file
    ifs.open(infilename.c_str());
    if (!ifs)
    {
        std::cerr << program_name << ": cannot read input file!" << std::endl;
        exit (-1);
    }

    OpenBabel::OBMol mol;
    if (ob_conv.Read(&mol, &ifs)) OutputMolecule(&mol);

    return 0;
}
