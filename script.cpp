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

#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <fstream>
using namespace std;

int main(int argc, char** argv)
{
	string sysCommand, scenarioFileName, scenarioName, moleculeFileName;
	ifstream scenarioFile, moleculesFile;
	ofstream logFile;
	string synthCommand = "./synth";

	if (argc <= 1)
    {
        cout << "Please include the protein name to test" << endl;
        return 1;
    }

	/* get name of protein */
	string nameOfProtein = string(argv[1]);

	/* make new directory in outputFiles */
	sysCommand = "mkdir outputFiles/" + nameOfProtein + "\n";
	cout << sysCommand << endl;
	system(sysCommand.c_str());

	/* unpack protein tar into chemistry folder */
	sysCommand = "tar -xf moleculeLib/output-" + nameOfProtein + ".tar \n";
	cout << sysCommand << endl;
	system(sysCommand.c_str());

	/* get list of scenarios */
	sysCommand = "ls output-" + nameOfProtein + "/ >> deleteme.txt\n";
	cout << sysCommand << endl;
	system(sysCommand.c_str());
	scenarioFile.open("deleteme.txt");

	while(scenarioFile >> scenarioFileName)
	{ /* for each scenario */

		/* get scenario name */
		scenarioName = scenarioFileName.substr(0, scenarioFileName.find("."));

		/* unpack scenario tar into chemistry folder */
		sysCommand = "tar -xf output-" + nameOfProtein + "/" + scenarioFileName + "\n";
		cout << sysCommand << endl;
		system(sysCommand.c_str());

		/* remove all molecule files from chemistry folder */
		//system("rm *.sdf\n");

		/* move linker files into chemistry folder */
		sysCommand = "mv " + scenarioName + "/linkers/*.sdf* .\n";
		cout << sysCommand << endl;
		system(sysCommand.c_str());

		/* move rigid files into chemistry folder */
		sysCommand = "mv " + scenarioName + "/rigids/*.sdf* .\n";
		cout << sysCommand << endl;
		system(sysCommand.c_str());

		/* generate synth command */
		
		system("ls *.sdf* >> moleculeFileList.txt\n");
		moleculesFile.open("moleculeFileList.txt");
		cout << endl << endl;
		while(moleculesFile >> moleculeFileName)
		{
			synthCommand += " " + moleculeFileName;
		}

		synthCommand += " -o output-" + scenarioName + ".sdf\n";

		moleculesFile.close();
		system("rm moleculeFileList.txt");

		/* run synth */
		cout << synthCommand << endl;
		system(synthCommand.c_str());

		/* move output file into outputFiles/[nameOfProtein] */
		sysCommand = "mv output-" + scenarioName + ".sdf outputFiles/" + nameOfProtein + "/output-" + scenarioName + ".sdf\n";
		cout << sysCommand << endl;
		system(sysCommand.c_str());

		/* remove linkers and rigids from chemistry folder */
		system("rm *.sdf");
		

		/* remove scenario folder */
		sysCommand = "rm -rf " + scenarioName + "\n";
		cout << sysCommand << endl;
		system(sysCommand.c_str());
		
		/*sysCommand = "rmdir " + scenarioName + "\n";
		system(sysCommand.c_str());*/
		/*sysCommand = "rm " + nameOfProtein + "/" + scenarioFileName + ".tar\n";
		system(sysCommand.c_str());*/

	} /* for each scenario */

	scenarioFile.close();
	system("rm deleteme.txt");

	/* remove protein folder */
	sysCommand = "rm -rf output-" + nameOfProtein + "\n";
	cout << sysCommand << endl;
	system(sysCommand.c_str());
	

	/* remove protein tar */
	//sysCommand = "rm moleculeLib/output-" + nameOfProtein + ".tar \n";

	/* make tar file from outputFiles/[nameOfProtein] */
	sysCommand = "tar -czf outputFiles/" + nameOfProtein + ".tar.gz outputFiles/" + nameOfProtein + "\n";
	cout << sysCommand << endl;
	system(sysCommand.c_str());

	logFile.open("outputLog.txt");
	logFile << synthCommand << endl;

	cout << "Done." << endl;

    return 0;
}
