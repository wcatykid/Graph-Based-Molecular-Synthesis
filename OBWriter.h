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

#ifndef _OB_WRITER_GUARD
#define _OB_WRITER_GUARD 1


#include <vector>
#include <string>
#include <iostream>
#include <queue>
#include <pthread.h>


#include <openbabel/mol.h>


#include "Molecule.h"
#include "Thread_Pool.h"
#include "IdFactory.h"


//
// A class to dump all of the molecules to a file.
//
class OBWriter
{
  public:
     OBWriter(unsigned int numThreads);
    ~OBWriter();

    // static void InitializeFile(const char* fileName);
    void OutputMoleculeInternal(unsigned int, unsigned int, Molecule&);
    void OutputMoleculeExternalSMI(Molecule&);
    void OutputMoleculeExternalSDF(Molecule&);
    void OutputMoleculeAppendExternalSMI(const std::string& smi);
    void OutputMoleculeAppendExternalSDF(Molecule&);
    static int OutputSingleMolecule(std::string smiMol);
    static std::vector<OpenBabel::OBMol*> compliantMols;

    void IndicateSynthesisStarted();
    void IndicateSynthesisComplete();
    void InitiateOutputThreadPool();
    void IndicateSMIwritingComplete() const;

    static void InitializeFile(const std::string& outFile);

    void write(std::vector<Molecule> molecules);

    static void TurnValidationOff() { performValidation = false; }
    static void SetPool(Thread_Pool<std::string, int>* thePool) { staticPool = thePool; }
    static unsigned InputPoolSize() { return staticPool->in_q_size(); }
    static unsigned OutputPoolSize() { return staticPool->out_q_size(); }
    static unsigned NumCompliantMolecules() { return numCompliant; }

    static void ScrubAndConvertToSMIInternal(OpenBabel::OBMol* mol, std::string& smi);
    static void ScrubAndConvertToSMIExternal(OpenBabel::OBMol* mol, std::string& smi);
    static void ConvertToSMI(const std::string& sdf, std::string& smi);

  private:
    unsigned int mCounter; 
    unsigned int mFailCounter; 
    bool writing_complete;
    bool writing_started;

    static bool synthesis_complete;
    static bool performValidation;
    static pthread_mutex_t smi_output_file_lock;
    static pthread_mutex_t sdf_output_file_lock;
    static pthread_mutex_t id_lock;
    static pthread_mutex_t popen_lock;
    static pthread_mutex_t smi_popen_lock;
    static pthread_mutex_t writer_popen_lock;
    static pthread_mutex_t valid_molecule_lock;
    static std::ofstream out;
    static std::string outFileName;
    static unsigned numCompliant;
    static OpenBabel::OBConversion SDF_to_SMI_conv;


    Thread_Pool<std::string, int>* pool;  

    // The same static version of the pool.
    static Thread_Pool<std::string, int>* staticPool;

    void Initialize();

    void ScrubAndExportSMI(std::vector<Molecule>& molecules);
    void CallsBeforeWriting(std::vector<Molecule>& molecules);

    unsigned molCounter;
    std::string prefix;
    std::string sdfSuffix;
    std::string smiSuffix;
    unsigned UPPERBOUND;
    std::string infix;
    std::string outputDir;
    std::string sdfOutfileName;
    std::string smiOutfileName;
};

#endif
