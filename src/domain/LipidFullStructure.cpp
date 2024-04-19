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


#include "cppgoslin/domain/LipidFullStructure.h"

LipidFullStructure::LipidFullStructure(Headgroup* _headgroup, vector<FattyAcid*> *_fa) : LipidStructureDefined(_headgroup, _fa) {            
    info->level = FULL_STRUCTURE;
}


LipidLevel LipidFullStructure::get_lipid_level(){
    return FULL_STRUCTURE;
}


string LipidFullStructure::get_lipid_string(LipidLevel level){
    switch(level){
        case NO_LEVEL:
        case FULL_STRUCTURE:
            return LipidMolecularSpecies::build_lipid_subspecies_name(FULL_STRUCTURE);
            
        case STRUCTURE_DEFINED:
        case SN_POSITION:
        case MOLECULAR_SPECIES:
        case SPECIES:
        case CATEGORY:
        case CLASS:
            return LipidStructureDefined::get_lipid_string(level);
    
        default:
            throw IllegalArgumentException("LipidFullStructure does not know how to create a lipid string for level " + std::to_string(level));
    }
}
