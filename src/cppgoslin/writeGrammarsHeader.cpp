/*
MIT License

Copyright (c) the authors (listed in global LICENSE file)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


#include <fstream>
#include <iostream>
#include <stdlib.h>

using namespace std;

void replaceAll(std::string& str, const std::string& from, const std::string& to) {
    if(from.empty())
        return;
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }
}



void addingGrammar(ofstream& offile, string grammarName, string grammarFilename){

    
    ifstream infile(grammarFilename.c_str());
    offile << "const string " + grammarName + " = \"";
    if (!infile.good()){
        cout << "Error: file '" + grammarFilename + "' not found." << endl;
        exit(-1);
    }
    
    string line;
    while (getline(infile, line)){
        replaceAll(line, "\\", "\\\\");
        replaceAll(line, "\"", "\\\"");
        line += " \\n\\";
        offile << line << endl;
    }
    offile << "\";" << endl;
}



void writeGrammarHeader(string ofFileName){
    ofstream offile(ofFileName.c_str());
    
    offile << "/*" << endl;
    offile << "MIT License" << endl;
    offile << endl;
    offile << "Copyright (c) the authors (listed in global LICENSE file)" << endl;
    offile << endl;
    offile << "Permission is hereby granted, free of charge, to any person obtaining a copy" << endl;
    offile << "of this software and associated documentation files (the \"Software\"), to deal" << endl;
    offile << "in the Software without restriction, including without limitation the rights" << endl;
    offile << "to use, copy, modify, merge, publish, distribute, sublicense, and/or sell" << endl;
    offile << "copies of the Software, and to permit persons to whom the Software is" << endl;
    offile << "furnished to do so, subject to the following conditions:" << endl;
    offile << endl;
    offile << "The above copyright notice and this permission notice shall be included in all" << endl;
    offile << "copies or substantial portions of the Software." << endl;
    offile << endl;
    offile << "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR" << endl;
    offile << "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY," << endl;
    offile << "FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE" << endl;
    offile << "AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER" << endl;
    offile << "LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM," << endl;
    offile << "OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE" << endl;
    offile << "SOFTWARE." << endl;
    offile << "*/" << endl;
    offile << endl;
    offile << endl;
    offile << endl;
    
    offile << "#ifndef KNOWN_GRAMMARS_H" << endl;
    offile << "#define KNOWN_GRAMMARS_H" << endl << endl;

    offile << "#include <string>" << endl;

    offile << "using namespace std;" << endl;
    
    addingGrammar(offile, "shorthand_grammar", "data/goslin/Shorthand2020.g4");
    offile << endl << endl << endl;
    
    addingGrammar(offile, "goslin_grammar", "data/goslin/Goslin.g4");
    offile << endl << endl << endl;
    
    addingGrammar(offile, "lipid_maps_grammar", "data/goslin/LipidMaps.g4");
    offile << endl << endl << endl;
    
    addingGrammar(offile, "swiss_lipids_grammar", "data/goslin/SwissLipids.g4");
    offile << endl << endl << endl;
    
    addingGrammar(offile, "hmdb_grammar", "data/goslin/HMDB.g4");
    offile << endl << endl << endl;
    
    addingGrammar(offile, "sum_formula_grammar", "data/goslin/SumFormula.g4");
    offile << endl << endl << endl;
    
    addingGrammar(offile, "fatty_acid_grammar", "data/goslin/FattyAcids.g4");
    offile << endl << endl << endl;
    
    addingGrammar(offile, "systematic_grammar", "data/goslin/Systematic.g4");
    offile << endl << endl << endl;
    
    offile << "#endif /* KNOWN_GRAMMARS_H */" << endl;
}







int main(int argc, char** argv){
    if (argc < 2){
        cout << "Error, specify grammar output filename." << endl;
        exit(-1);
    }

    writeGrammarHeader(argv[1]);
    return 0;
}
