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
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <pthread.h>
#include <unistd.h>
#include <cstdio>
#include <sys/stat.h>
#include <dirent.h>



#include <openbabel/mol.h>
#include <openbabel/obconversion.h>



#include "OBWriter.h"
#include "Molecule.h"
#include "Rigid.h"
#include "Linker.h"
#include "obgen.h"
#include "Thread_Pool.h"
#include "Constants.h"
#include "Utilities.h"
#include "IdFactory.h"
#include "Options.h"
#include "zpipe.h"



// Static Definitions
pthread_mutex_t OBWriter::valid_molecule_lock;
pthread_mutex_t OBWriter::sdf_output_file_lock;
pthread_mutex_t OBWriter::smi_output_file_lock;
pthread_mutex_t OBWriter::id_lock;
pthread_mutex_t OBWriter::popen_lock;
pthread_mutex_t OBWriter::smi_popen_lock;
pthread_mutex_t OBWriter::writer_popen_lock;
//IdFactory OBWriter::molIDmaker(1000);
std::ofstream OBWriter::out;
std::string OBWriter::outFileName;
std::vector<OpenBabel::OBMol*> OBWriter::compliantMols;
bool OBWriter::synthesis_complete = false;
bool OBWriter::performValidation = true;
Thread_Pool<std::string, int>* OBWriter::staticPool = 0;
unsigned OBWriter::numCompliant = 0;
OpenBabel::OBConversion OBWriter::SDF_to_SMI_conv;


// ****************************************************************************

OBWriter::OBWriter(unsigned int threadCount) : mCounter(0),
                                               mFailCounter(0),
                                               writing_complete(false),
                                               writing_started(false)          
{
    // Create the thread pool
    if (!Options::SMI_ONLY)
    {
        pool = new Thread_Pool<std::string, int>(threadCount, OBWriter::OutputSingleMolecule);
    }

    molCounter = 0;
    prefix = "molecules";
    sdfSuffix = ".sdf";
    smiSuffix = ".smi";
    UPPERBOUND = 250000;
    outputDir = DEFAULT_OUTPUT_DIR;
    sdfOutfileName = outputDir + "/" + prefix + "-1-10000" + sdfSuffix;
    smiOutfileName = outputDir + "/" + prefix + "-1-250000" + smiSuffix;
}

// ****************************************************************************

OBWriter::~OBWriter()
{
    // Killing the thread pool to force all threads to join.
    if (!Options::SMI_ONLY) delete pool;
}

// ****************************************************************************

void OBWriter::InitializeFile(const std::string& outFile)
{
    outFileName = outFile;

    /*

    out.open(outFile.c_str());

    if (out.fail()) throw "Output stream opening failed.";
    */
}


// ****************************************************************************

void OBWriter::Initialize()
{
    pthread_mutex_init(&OBWriter::valid_molecule_lock, NULL);
    pthread_mutex_init(&OBWriter::smi_output_file_lock, NULL);
    pthread_mutex_init(&OBWriter::sdf_output_file_lock, NULL);
    pthread_mutex_init(&OBWriter::popen_lock, NULL);
}

// ****************************************************************************

void OBWriter::IndicateSynthesisStarted()
{
    OBWriter::SetPool(this->pool);

    //
    // Construct the actual directory
    //
    std::string theDir = DEFAULT_OUTPUT_DIR;
    if (Options::OUTPUT_DIR_SUFFIX != "")
    {
        theDir += "_";
        theDir += Options::OUTPUT_DIR_SUFFIX;
    }

    std::cout << "Will output to directory: " << theDir << std::endl;

    bool overwrite = true;
    if (DoesDirectoryExist(theDir))
    {
        std::cout << "The output directory " << theDir << " exists." << std::endl;
        std::cout << "Do you wish to overwrite? (Y / N)" << std::endl; 

        // Read in the entire input string and observe only the first character
        std::string decision;
        std::cin >> decision;

        if (decision[0] == 'y' || decision[0] == 'Y') overwrite = true;
        else overwrite = false;
    }
    else
    {
        MakeDirectory(theDir);
    }

    if (!overwrite)
    {
        std::cerr << "Re-run ./synth with different output directory." << std::endl;
        exit(1);
    }

//exit(1);

    // Remove all files in the directory
    CleanDirectory(theDir);

    // Set specific output file information.
    outputDir = theDir;
    sdfOutfileName = outputDir + "/" + prefix + "-1-10000" + sdfSuffix;
    smiOutfileName = outputDir + "/" + prefix + "-1-250000" + smiSuffix;    
}

// ****************************************************************************

void OBWriter::IndicateSynthesisComplete()
{
    synthesis_complete = true;

    if (Options::SMI_ONLY)
    {
        OBWriter::out.close();

        zlib_compress(smiOutfileName, smiOutfileName + ".zlib");
        remove(smiOutfileName.c_str());
    }
    else
    {
        std::cerr << "Synthesis is complete; writing continues." << std::endl;
        std::cerr << "Input pool contains " << pool->in_q_size()
                  << " molecules to process with obgen." << std::endl;

        // Spin until writing is complete.
        while (writing_started && mCounter > pool->out_q_size())
        {
            // Sleep 5 seconds; obgen takes a while.
            sleep(5);
        }
        std::cerr << "Writing of the molecules with obgen is complete." << std::endl;
    }
}

// ****************************************************************************

void OBWriter::OutputMoleculeInternal(unsigned int non_killed_sz, unsigned int current_sz, Molecule& mol)
{
    // If this is the first call to output, save the fact we are writing
    if (!writing_started)
    {
        if (!Options::SMI_ONLY) OBWriter::Initialize();
        writing_started = true;
    }

    //
    // Output molecule
    //

    //
    // Process the molecule for output
    //
    // (1) Lock around open babel
static unsigned int num_blocked = 0;
num_blocked++;

std::cerr << "Num Blocked: " << num_blocked << std::endl;

    pthread_mutex_lock(&Molecule::openbabel_lock);  
num_blocked--;
    // (2) make a copy with the copy constructor.
    OpenBabel::OBMol theMol; // Needs fixing *(mol.getOpenBabelMol());

    // (3) Check Lipinski Compliance
    // The molecule must be Lipinski compliant (using Open Babel)
    // We use the copy as not to disrupt the approximations we use during synthesis.
    if (!Molecule::isOpenBabelLipinskiCompliant(theMol))
    {
        this->mFailCounter++;
        pthread_mutex_unlock(&Molecule::openbabel_lock);
        return;
    }
    else this->mCounter++;

    pthread_mutex_unlock(&Molecule::openbabel_lock);

    // (4) Export to SMI
    std::string smiMol = "Needs Fixing.";
    // mol.getSMI(smiMol);

    // (5) Add the molecule to the queue for processing.
    if (Options::SMI_ONLY)
    {
        if (numCompliant % 10 == 0)
        {
            std::cerr << "Writing Lipinski compliant molecule " << numCompliant 
                      << "; hypergraph contains (" << non_killed_sz << ", "
                      << current_sz << ")" << std::endl; 
        }
        // output only the SMI version of the information to the output file.
        OBWriter::out << smiMol << std::endl;
    }
    else
    {
        pool->push(smiMol);
    }

    // Maintain a count
    OBWriter::numCompliant++;
}

// ****************************************************************************

void OBWriter::OutputMoleculeExternalSMI(Molecule& mol)
{
    // (a) create temp file
    char tmpFileNameBuff[32]; // temp file name
    int filedes = -1; // file descriptor

    // memset the buffer to 0 (thread-safe)
    memset(tmpFileNameBuff, 0, sizeof(tmpFileNameBuff));

    // Copy the dir/template the buffer// (thread-safe)
    std::string shmpath = Options::shmPath;
    shmpath += "chemTmpFile-XXXXXX.smi";
    strncpy(tmpFileNameBuff, shmpath.c_str(), shmpath.size());

    // generate unique non existing name from template, afterwards, buffer
    // contains actual file name(thread-safe)
    filedes = mkstemps(tmpFileNameBuff, 4);
    if(filedes < 1)
    {
        std::cerr << "Creation of temp file failed" << std::endl;
        return;
    }

    // (b) write SMI into file
    FILE *tmp_file; // file pointer
    tmp_file = fdopen(filedes, "w"); // 'front end' file pointer made from file descriptor
    std::string smi = "Needs fixing.";
    // mol.getSMI(smi);
    fputs(smi.c_str(), tmp_file); // write to tmp file (fputs requires file pointer)
    fclose(tmp_file); // no longer need front end, get rid of it

    //
    // Call external instance of compliance writer
    //
    // (a) set up function call, we suppress all stderr messages with 2> and
    // keep stdout with output.
    std::string call = Options::writerPath;
    call += COMPLIANT_EXE + " ";
    call += string(tmpFileNameBuff) + " molecules.smi";
    std::cerr << "Calling: " << call << std::endl;

    // (b) Spawn obgen process to work on the temp file. resulting SDF in string 'result'
    pthread_mutex_lock(&writer_popen_lock);
    FILE* pipe = popen(&call[0], "r"); // not clear whether this is thread safe
    pthread_mutex_unlock(&writer_popen_lock);

    if (!pipe)
    {
        std::cerr << "Creation of obgen caller pipe failed" << std::endl;
        return;
    }

    char data_buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
        if(fgets(data_buffer, 128, pipe) != NULL)
                result += data_buffer;
    }

    pthread_mutex_lock(&smi_popen_lock);
    pclose(pipe); // not clear wether is thread safe
    pthread_mutex_unlock(&smi_popen_lock);

    // (c) close and unlink temporary file
    close(filedes);
    unlink(tmpFileNameBuff); // even if called, file won't actually delete until file is closed
    // Removal of the temporary file.
    std::cerr << "Removing: " << tmpFileNameBuff << std::endl;
    remove(tmpFileNameBuff);
}

// ****************************************************************************

void OBWriter::OutputMoleculeAppendExternalSMI(const std::string& smi)
{
    //
    // Append an SMI version of the molecule to the output file.
    //
    pthread_mutex_lock(&smi_output_file_lock);

    molCounter++;

    //
    // Update the file we are writing to.
    //
    if (molCounter % UPPERBOUND == 0)
    {
        // zip the file we just created and remove it.
        zlib_compress(smiOutfileName, smiOutfileName + ".zlib");
        remove(smiOutfileName.c_str());

        //
        // Create the new output file name.
        //
        std::ostringstream oss;

        oss << outputDir << "/" << prefix << "-" << molCounter
            << "-" << (molCounter + UPPERBOUND) << smiSuffix;

        smiOutfileName = oss.str();
    }


    std::ofstream outfile(smiOutfileName.c_str(),
                          std::ofstream::out | std::ofstream::app); // append

    outfile << smi << std::endl;

    outfile.close();

    pthread_mutex_unlock(&smi_output_file_lock);
}

// ****************************************************************************

void OBWriter::IndicateSMIwritingComplete() const
{
    std::cout << "Compressing the smi file: " << smiOutfileName << std::endl;

    // zip the file we just created and remove it.
    //zlib_compress(smiOutfileName, smiOutfileName + ".zlib");
    std::string zip = "gzip " + smiOutfileName; 
    system(zip.c_str());
}

// ****************************************************************************

void OBWriter::OutputMoleculeAppendExternalSDF(Molecule& mol)
{
    //
    // Create the SDF format
    //
    std::string sdf;
    mol.WriteToOpenBabelFormat(sdf);
   
    //
    // Append an SDF version of the molecule to the output file.
    //
    pthread_mutex_lock(&sdf_output_file_lock);

    molCounter++;

    //
    // Update the file we are writing to.
    //
    if (molCounter % UPPERBOUND == 0 && molCounter > UPPERBOUND)
    {
        // zip the file we just created and remove it.
        zlib_compress(sdfOutfileName, sdfOutfileName + ".zlib");
        remove(sdfOutfileName.c_str());

        //
        // Create the new output file name.
        //
        std::ostringstream oss;

        oss << outputDir << "/" << prefix << "-" << (molCounter - UPPERBOUND)
            << "-" << molCounter << sdfSuffix;

        sdfOutfileName = oss.str();
    }

    std::ofstream outfile(sdfOutfileName.c_str(),
                          std::ofstream::out | std::ofstream::app); // append

    outfile << sdf << std::endl << std::endl;
    outfile << "$$$$" << std::endl << std::endl;

    outfile.close();

    pthread_mutex_unlock(&sdf_output_file_lock);
}

// ****************************************************************************

void OBWriter::OutputMoleculeExternalSDF(Molecule& mol)
{

    // (a) create temp file
    char tmpFileNameBuff[32]; // temp file name
    int filedes = -1; // file descriptor

    // memset the buffer to 0 (thread-safe)
    memset(tmpFileNameBuff, 0, sizeof(tmpFileNameBuff));

    pthread_mutex_lock(&sdf_output_file_lock);

    // Copy the dir/template the buffer// (thread-safe)
    std::string shmpath = Options::shmPath;
    shmpath += "chemTmpFile-XXXXXX.sdf";
    strncpy(tmpFileNameBuff, shmpath.c_str(), shmpath.size());

    // generate unique non existing name from template, afterwards, buffer
    // contains actual file name(thread-safe)
    filedes = mkstemps(tmpFileNameBuff, 4);

    if(filedes < 1)
    {
        std::cerr << "Creation of temp file failed" << std::endl;
        return;
    }


    // (b) write SDF into file
    FILE *tmp_file; // file pointer
    tmp_file = fdopen(filedes, "w"); // 'front end' file pointer made from file descriptor
    std::string sdf;
    mol.WriteToOpenBabelFormat(sdf);

    if (tmp_file == NULL) 
    {
        std:: cerr << "|" << tmpFileNameBuff << "|" << std::endl;
        throw "Problem.";
    }
    fputs(sdf.c_str(), tmp_file); // write to tmp file (fputs requires file pointer)
    fclose(tmp_file); // no longer need front end, get rid of it

    //
    // Call external instance of compliance writer
    //
    // (a) set up function call, we suppress all stderr messages with 2> and
    // keep stdout with output.
    std::string call = Options::writerPath;
    call += COMPLIANT_EXE + " ";
    call += string(tmpFileNameBuff) + " molecules.smi";
    std::cerr << "Calling: " << call << std::endl;

    // (b) Spawn obgen process to work on the temp file. resulting SDF in string 'result'
    pthread_mutex_lock(&writer_popen_lock);
    FILE* pipe = popen(&call[0], "r"); // not clear whether this is thread safe
    pthread_mutex_unlock(&writer_popen_lock);

    if (!pipe)
    {
        std::cerr << "Creation of obgen caller pipe failed" << std::endl;
        return;
    }

    char data_buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
        if(fgets(data_buffer, 128, pipe) != NULL)
                result += data_buffer;
    }

    pthread_mutex_lock(&smi_popen_lock);
    pclose(pipe); // not clear wether is thread safe
    pthread_mutex_unlock(&smi_popen_lock);

    // (c) close and unlink temporary file
    close(filedes);
    unlink(tmpFileNameBuff); // even if called, file won't actually delete until file is closed

    // Removal of the temporary file.
    std::cerr << "Removing: " << tmpFileNameBuff << std::endl; 
    remove(tmpFileNameBuff);

    pthread_mutex_unlock(&sdf_output_file_lock);
}


// ****************************************************************************

int OBWriter::OutputSingleMolecule(std::string smiMol)
{
    // Output debugging information / progress bar
    unsigned inPoolSize = OBWriter::InputPoolSize();
    unsigned outPoolSize = OBWriter::OutputPoolSize();

    std::cerr << "Pool IN queue (" << inPoolSize
              << "); OUT queue (" << outPoolSize << ")" << std::endl;

    // 
    // Write SMI to a temp file
    //

    // (a) create temp file
    char tmpFileNameBuff[32]; // temp file name
    int filedes = -1; // file descriptor

    // memset the buffer to 0 (thread-safe)
    memset(tmpFileNameBuff, 0, sizeof(tmpFileNameBuff));

    // Copy the dir/template the buffer// (thread-safe)
    std::string shmpath = Options::shmPath;
    shmpath += "chemTmpFile-XXXXXX.smi";
    strncpy(tmpFileNameBuff, shmpath.c_str(), shmpath.size());

    // generate unique non existing name from template, afterwards, buffer
    // contains actual file name(thread-safe)
    filedes = mkstemps(tmpFileNameBuff, 4);
    if(filedes < 1)
    {
        std::cerr << "Creation of temp file failed" << std::endl;
        return -1;
    }

    // (b) write SMI into file
    FILE *tmp_file; // file pointer
    tmp_file = fdopen(filedes, "w"); // 'front end' file pointer made from file descriptor
    fputs(smiMol.c_str(), tmp_file); // write to tmp file (fputs requires file pointer)
    fclose(tmp_file); // no longer need front end, get rid of it

    //
    // Call external instance of obgen
    //
    // (a) set up function call, we suppress all stderr messages with 2> and
    // keep stdout with output.
    std::string obgenCall = "./synthobgen " + string(tmpFileNameBuff) + " 2> /dev/null";
    std::cerr << "Calling: " << obgenCall << std::endl;

    // (b) Spawn obgen process to work on the temp file. resulting SDF in string 'result'
    pthread_mutex_lock(&popen_lock);
    FILE* pipe = popen(&obgenCall[0], "r"); // not clear wether is thread safe
    pthread_mutex_unlock(&popen_lock);

    if (!pipe)
    {
        std::cerr << "Creation of obgen caller pipe failed" << std::endl;
        return -1;
    }

    char data_buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
    	if(fgets(data_buffer, 128, pipe) != NULL)
    		result += data_buffer;
    }

    pthread_mutex_lock(&popen_lock);
    pclose(pipe); // not clear wether is thread safe
    pthread_mutex_unlock(&popen_lock);

    // (c) close and unlink temporary file
    close(filedes);
    unlink(tmpFileNameBuff); // even if called, file won't actually delete until file is closed

    //
    // (d) Append output to a total output file.
    //
    pthread_mutex_lock(&sdf_output_file_lock);
    OBWriter::out << result;
    pthread_mutex_unlock(&sdf_output_file_lock);

    //
    // (e) We keep the SDF version of the synthesized molecule; only for validation purposes
    //
    if (OBWriter::performValidation)
    {
        // Begin open babel usage
        pthread_mutex_lock(& Molecule::openbabel_lock);

        // store result into a new molecule
        OpenBabel::OBMol* mol = new OpenBabel::OBMol();
        OpenBabel::OBConversion SDF_conv;
        SDF_conv.SetInAndOutFormats("SDF", "SDF");
        SDF_conv.ReadString(mol, result);

        // End open babel usage
        pthread_mutex_unlock(& Molecule::openbabel_lock);

        // Save the valid molecule for validation purposes
        pthread_mutex_lock(& OBWriter::valid_molecule_lock);
        OBWriter::compliantMols.push_back(mol);
        pthread_mutex_unlock(& OBWriter::valid_molecule_lock);
    }

    // Maintain a count
    OBWriter::numCompliant++;

    std::cerr << tmpFileNameBuff << " completed." << std::endl;

    return 0;
}

///////////////////////////////////////////////////////////////////

void OBWriter::ConvertToSMI(const std::string& sdf, std::string& smi)
{
    // Begin open babel usage
    pthread_mutex_lock(& Molecule::openbabel_lock);

    // store sdf into an open babel molecule
    OpenBabel::OBMol mol;

    // Use the static converter in OBWriter
    SDF_to_SMI_conv.SetInAndOutFormats("SDF", "SMI");

    // Convert to SDF
    SDF_to_SMI_conv.ReadString(&mol, sdf);

    // Convert to SMI
    smi = SDF_to_SMI_conv.WriteString(&mol);

    // End open babel usage
    pthread_mutex_unlock(& Molecule::openbabel_lock);

    // Clean up the smi value; ensures only molecule values
    smi = smi.substr(0, smi.find('\t'));
}


// ****************************************************************************

void OBWriter::ScrubAndConvertToSMIInternal(OpenBabel::OBMol* mol, std::string& smi)
{
    // (2) make a copy with the copy constructor.
    OpenBabel::OBMol theMol = *(mol);

    OpenBabel::OBConversion SMI_conv(&std::cin, &std::cout);

    // set conversion type(s) and verify it worked
    if(!SMI_conv.SetInAndOutFormats("SMI","SMI"))
        throw "SetInAndOutFormats failed!";

    // Pre-emptive extra run of OBGen, seems to stop segmentation fault
    OBGen::fast_obgen(&theMol);

    smi = SMI_conv.WriteString(&theMol); // convert and write to string

    // remove garbage and acquire only the SMI
    smi = smi.substr(0, smi.find('\t'));
}


// ****************************************************************************

void OBWriter::ScrubAndConvertToSMIExternal(OpenBabel::OBMol* mol, std::string& smi)
{

smi = "TBD";

    //
    // Write SDF to a temp file
    //
/*
    // (a) create temp file
    char tmpFileNameBuff[32]; // temp file name
    int filedes = -1; // file descriptor

    // memset the buffer to 0 (thread-safe)
    memset(tmpFileNameBuff, 0, sizeof(tmpFileNameBuff));

    pthread_mutex_lock(&sdf_output_file_lock);

    // Copy the dir/template the buffer// (thread-safe)
    strncpy(tmpFileNameBuff, "/run/shm/chemTmpFile-XXXXXX.smi", 31);

    // generate unique non existing name from template, afterwards, buffer
    // contains actual file name(thread-safe)
    filedes = mkstemps(tmpFileNameBuff, 4);
    if(filedes < 1)
    {
        std::cerr << "Creation of temp file failed" << std::endl;
        return;
    }

    // (b) write SMI into file
    FILE *tmp_file; // file pointer
    tmp_file = fdopen(filedes, "w"); // 'front end' file pointer made from file descriptor
if (tmp_file == NULL) throw "Problem.";
std:: cerr << "|" << smiMol << "|" << std::endl;
    fputs(smiMol.c_str(), tmp_file); // write to tmp file (fputs requires file pointer)
    pthread_mutex_unlock(&sdf_output_file_lock);
    fclose(tmp_file); // no longer need front end, get rid of it

    //
    // Call external instance of obgen
    //
    // (a) set up function call, we suppress all stderr messages with 2> and
    // keep stdout with output.
    std::string call = "./convertToSMI " + string(tmpFileNameBuff) + " 2> /dev/null";
    std::cerr << "Calling: " << call << std::endl;

    // (b) Spawn obgen process to work on the temp file. resulting SDF in string 'result'
    pthread_mutex_lock(&smi_popen_lock);
    FILE* pipe = popen(&obgenCall[0], "r"); // not clear wether is thread safe
    pthread_mutex_unlock(&smi_popen_lock);

    if (!pipe)
    {
        std::cerr << "Creation of obgen caller pipe failed" << std::endl;
        return;
    }

    char data_buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
        if(fgets(data_buffer, 128, pipe) != NULL)
                result += data_buffer;
    }

    pthread_mutex_lock(&smi_popen_lock);
    pclose(pipe); // not clear wether is thread safe
    pthread_mutex_unlock(&smi_popen_lock);

    // (c) close and unlink temporary file
    close(filedes);
    unlink(tmpFileNameBuff); // even if called, file won't actually delete until file is closed

    //
    // (d) Append output to a total output file.
    //
    pthread_mutex_lock(&sdf_output_file_lock);
    OBWriter::out << result;
    pthread_mutex_unlock(&sdf_output_file_lock);

    //
    // (e) We keep the SDF version of the synthesized molecule; only for validation purposes
    //
    if (OBWriter::performValidation)
    {
        // Begin open babel usage
        pthread_mutex_lock(& Molecule::openbabel_lock);

        // store result into a new molecule
        OpenBabel::OBMol* mol = new OpenBabel::OBMol();
        OpenBabel::OBConversion SDF_conv;
        SDF_conv.SetInAndOutFormats("SDF", "SDF");
        SDF_conv.ReadString(mol, result);

        // End open babel usage
        pthread_mutex_unlock(& Molecule::openbabel_lock);

        // Save the valid molecule for validation purposes
        pthread_mutex_lock(& OBWriter::valid_molecule_lock);
        OBWriter::compliantMols.push_back(mol);
        pthread_mutex_unlock(& OBWriter::valid_molecule_lock);
    }

    // Maintain a count
    OBWriter::numCompliant++;

    std::cerr << tmpFileNameBuff << " completed." << std::endl;

    return 0;
*/
}

// ****************************************************************************

void OBWriter::ScrubAndExportSMI(std::vector<Molecule>& molecules)
{
    OpenBabel::OBConversion SDF_conv(&std::cin, &std::cout); // conversion to/from SDF (has xyz coords)
    OpenBabel::OBConversion SMI_conv(&std::cin, &std::cout); // conversion to/from SMI (no xyz coords)

    // set conversion type(s) and verify it worked
    if(!SMI_conv.SetInAndOutFormats("SMI","SMI") || !SDF_conv.SetInAndOutFormats("SDF","SDF"))
    {
        std::cerr << "SetInAndOutFormats failed!" << std::endl;
        return;    
    }

    std::string s; // temporary buffer
    int i = 1;     // counter (for debugging output)

    std::ofstream logfile("ScrubAndExportSMI_logfile.txt", std::ofstream::out);

    if (logfile.is_open())
    {
        std::cout  << "ScrubAndExportSMI: logfile.is_open";
    }
    else
    {
        std::cerr << "ScrubAndExportSMI: Unable to open logfile";
    }

    //
    // Iterate through all molecules.
    //
    for (std::vector<Molecule>::iterator it = molecules.begin(); it != molecules.end(); it++)
    {
        // Set up data struct and display basic information
        if (g_debug_output) std::cout << "ScrubAndExportSMI: getOpenBabelMol..." << std::endl;

        OpenBabel::OBMol* mol = 0; //it->getOpenBabelMol();  // OBMol mol;

        std::cout << "ScrubAndExportSMI: #" << i 
                  << ": NumAtoms=" << mol->NumAtoms()
                  << ", NumBonds=" << mol->NumBonds() << std::endl;
        logfile << "ScrubAndExportSMI: #" << i
                << ": NumAtoms=" << mol->NumAtoms()
                << ", NumBonds=" << mol->NumBonds() << std::endl;
        //std::cout << it->toString() << std::endl; // more detailed data

        // Pre-emptive extra run of OBGen, seems to stop segmentation fault
        if (g_debug_output) std::cout << "ScrubAndExportSMI: OBGen(fast)..." << std::endl;
        OBGen::fast_obgen(mol);

        // log "before" molecule
        //logfile << "molecule #" << i << ", SDF format, before:" << std::endl;
        //s=SDF_conv.WriteString(mol);
        //logfile << s;

        // Write to then read from SMI; should remove xyz coords
        if (g_debug_output) std::cout << "ScrubAndExportSMI: WriteString:" << std::endl;

        s = SMI_conv.WriteString(mol);

        logfile << "molecule #" << i << ", SMI: " << s; // log "SMI" version

        if (SMI_conv.ReadString(mol, s))
        {
            if (g_debug_output) std::cout << "ScrubAndExportSMI: ReadString - successful" << std::endl;

            // log "after" molecule
            if (g_debug_output) std::cout << "molecule #" << i << ", SDF format, after:" << std::endl;
                
            //logfile << "molecule #" << i << ", SDF format, after:" << std::endl;
            s = SDF_conv.WriteString(mol);

            if (g_debug_output) std::cout << s;
            //logfile << s;
        }
        else if (g_debug_output) std::cout << "ScrubAndExportSMI: ReadString - failed" << std::endl;

        if (g_debug_output)
        {
            std::cout << "-----------------------------------------------------------------" << std::endl;
        }

        logfile << "-----------------------------------------------------------------" << std::endl;
        i++;
    }

    logfile.close();
}

// ****************************************************************************

void OBWriter::CallsBeforeWriting(std::vector<Molecule>& molecules)
{
    int counter = 0;
    for (std::vector<Molecule>::iterator it = molecules.begin(); it != molecules.end(); it++)
    {
        OpenBabel::OBMol* obmol = 0; //it->getOpenBabelMol();

        // if (g_debug_output)
        std::cerr << "Converting molecule " << counter++ << " with obgen" << std::endl; 

        OBGen::obgen(obmol);
    }
}


// ****************************************************************************

void OBWriter::write(std::vector<Molecule> molecules)
{
    //
    // Refine the list of molecules to those that are lipinski compliant.
    //
    std::vector<Molecule> synthMolecules;
    for (std::vector<Molecule>::iterator it = molecules.begin(); it != molecules.end(); it++)
    {
        //it->openBabelPredictLipinski();

        if (it->IsComplex() && it->isLipinskiCompliant())
        {
            synthMolecules.push_back(*it);
        }
    }

    //
    // Convert to SMI and then use obgen to regenerate 3D molecules.
    //
    std::cout << "OBWriter::write: ScrubAndExportSMI(synthMolecules)..." << std::endl;
    ScrubAndExportSMI(synthMolecules);

    std::cout << "OBWriter::write: CallsBeforeWriting(synthMolecules)..." << std::endl;
    CallsBeforeWriting(synthMolecules);

    // Converter to output the synthesized molecules
    OpenBabel::OBConversion toSDF(&std::cin, &this->out);
    toSDF.SetOutFormat("SDF");

    //
    // Process all refined molecules
    //
    for (std::vector<Molecule>::iterator it = synthMolecules.begin();
         it != synthMolecules.end();
         it++)
    {
        //
        // Print the molecule number
        //
        this->out << "#### ";
        this->out << mCounter++;
        this->out << " ####";

        //
        // Print the names of all the linkers / rigids.
        //
    /*
        std::vector<Rigid*> rigids;
        it->getRigids(rigids);
        std::vector<Linker*> linkers;
        it->getLinkers(linkers);

        foreach_rigids(r_it, rigids)
        {
            this->out << (*r_it)->getName() << std::endl;
        }

        foreach_linkers(l_it, linkers)
        {
            this->out << (*l_it)->getName() << std::endl;
        }
    */

        //
        // Take the SMI and convert it back to SDF by populating the coordinates.
        //
        OpenBabel::OBMol* obmol = 0; // it->getOpenBabelMol();
        toSDF.Write(obmol);
    }
}

// ****************************************************************************
