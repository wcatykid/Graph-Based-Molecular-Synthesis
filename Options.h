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

#ifndef _OPTIONS_GUARD
#define _OPTIONS_GUARD 1

#include <string>
#include <vector>

//
// A aggregation class for options specified in the command-line
//
class Options
{
  public:
    Options(int, char**);
    bool parseCommandLine();

    bool AnalyzeEnvironment();
    static std::string writerPath;
    static std::string shmPath;

    std::string outFile;
    std::string outFileSMI;
    std::string validationFile;
    std::vector<std::string> inFiles;

    static double TANIMOTO;
    static bool THREADED;
    static bool SERIAL;
    static bool OPENBABEL;
    static bool USE_LIPINSKI;
    //static unsigned SMI_LEVEL_BOUND;
    static unsigned PROBABILITY_PRUNE_LEVEL_START;
    static unsigned int OBGEN_THREAD_POOL_SIZE;
    static bool SMI_ONLY;
    static std::string OUTPUT_DIR_SUFFIX;

  private:
    int argc;
    char** argv;

    bool handleOption(int& index);

    bool acquireEnvironmentVariable(const std::string& variable,
                                    const std::string& suffix,
                                    std::string& value);
};

#endif
