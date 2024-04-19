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


#include "cppgoslin/cppgoslin.h"
#include <set>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <cassert>
#include <fstream>
#include <math.h>
#include <cstdint>

using namespace std;
using namespace goslin;

int main(int argc, char** argv){
    char PARSER_QUOTE = '\'';
    
    
    LipidAdduct* lipid;
    string lipid_name;
    
    GoslinParser goslin_parser;
    LipidParser lipid_parser;
    LipidMapsParser lipid_maps_parser;
    //GoslinFragmentParser goslin_fragment_parser;
    SwissLipidsParser swiss_lipids_parser;
    HmdbParser hmdb_parser;
    int num_hydroxyl = 0;
    
    // test bitfield
    srand(time(0));
    int num_numbers = 1000 + (rand() % 1000);
    int max_number = num_numbers * 100;
    set<int> check;
    Bitfield b(max_number);
    
    for (int i = 0; i < num_numbers; ++i){
        int num = rand() % max_number;
        b.insert(num);
        check.insert(num);
    }
    
    uint32_t cnt = 0;
    for (int i : b){
        assert (check.find(i) != check.end());
        ++cnt;
    }
    assert (cnt == check.size());
    
    
    
    
    // check lipids maps parser for headless FA
    lipid_name = "15:6(2Z,4E,6Z,8E,12E,14)(6Me,8Me,10Me[S],13Me)";
    lipid = lipid_maps_parser.parse(lipid_name);
    assert (lipid_maps_parser.parser_event_handler->word_in_grammar);
    assert (lipid);
    delete lipid;
    
    
    
    // Pure Parser test
    GoslinParserEventHandler goslin_parser_event_handler;
    Parser<LipidAdduct*> goslin_parser_pure(&goslin_parser_event_handler, "data/goslin/Goslin.g4", PARSER_QUOTE);
    
    // check goslin pure parser with illegal lipid name
    string failLipid = "TAG 16::1-18:1-24:0";
    try {
        lipid = goslin_parser_pure.parse(failLipid);
        assert (false);
    }
    catch (LipidException &e){ }

    // glycerophospholipid
    for (auto the_lipid_name : { "TAG 16:1-18:1-24:0", "PE 16:1/12:0", "DAG 16:1-12:0", "12-HETE", "15S-HETE-d8", "HexCer 18:1;2/16:0"}){
        LipidAdduct* lipid = goslin_parser_pure.parse(the_lipid_name);
        assert (lipid);
        delete lipid;
    }
    
    
    
    
    // check lipid maps and swiss lipids parser with illegal lipid name
    string failLipidSL = "TG(16::1_18:1_24:0)";
    try {
        lipid = lipid_maps_parser.parse(failLipidSL);
        assert (false);
    }
    catch (LipidException &e){ }
    try {
        lipid = swiss_lipids_parser.parse(failLipidSL);
        assert (false);
    }
    catch (LipidException &e){ }
    
    
    /*
    failLipidSL = "TG 16::1-18:1-24:0";
    // check goslin fragment parser with illegal lipid name 
    try {
        lipid = goslin_fragment_parser.parse(failLipidSL);
        assert (false);
    }
    catch (LipidException &e){ }
    
    
    // check if goslin fragment parser parses correctly lipid name with fragment
    lipid_name = "PE 16:1-12:0 - -(H20)";
    lipid = goslin_fragment_parser.parse(lipid_name);
    assert(lipid);
    assert(lipid->fragment);
    assert(lipid->fragment->name == "-(H20)");
    delete lipid;
    */
    
    
    
    
    
    lipid = swiss_lipids_parser.parse("Cer(d18:1(4E)/24:0-2OH)");
    assert (lipid);
    assert(lipid->get_lipid_string() == "Cer 18:1(4);(OH)2/24:0;OH");
    assert(lipid->get_sum_formula() == "C42H83NO4");
    assert (abs(lipid->get_mass() - 665.632209) < 1e-3);
    delete lipid;
        
        
    lipid = swiss_lipids_parser.parse("Cer(d18:1(4E)/24:0(2OH))");
    assert (lipid);
    assert(lipid->get_lipid_string() == "Cer 18:1(4);(OH)2/24:0;OH");
    assert(lipid->get_sum_formula() == "C42H83NO4");
    assert (abs(lipid->get_mass() - 665.632209) < 1e-3);
    delete lipid;
    
        
        
    lipid = lipid_maps_parser.parse("Cer(d18:1(4E)/24:0(2OH))");
    assert (lipid);
    assert(lipid->get_lipid_string() == "Cer 18:1(4E);1OH,3OH/24:0;2OH");
    assert(lipid->get_sum_formula() == "C42H83NO4");
    assert (abs(lipid->get_mass() - 665.632209) < 1e-3);
    delete lipid;
        
    
    lipid = goslin_parser.parse("Cer 18:1(4E);2/24:0;1");
    assert (lipid);
    assert(lipid->get_lipid_string() == "Cer 18:1(4);(OH)2/24:0;OH");
    assert(lipid->get_sum_formula() == "C42H83NO4");
    assert (abs(lipid->get_mass() - 665.632209) < 1e-3);
    delete lipid;
        
    lipid = hmdb_parser.parse("SM(d18:1(9Z)/16:1(9Z)(OH))");
    assert (lipid);
    assert(lipid->get_lipid_string() == "SM 18:1(9);OH/16:1(9);OH");
    assert(lipid->get_sum_formula() == "C39H77N2O7P");
    assert (abs(lipid->get_mass() - 716.546841) < 1e-3);
    delete lipid;
    
    
    
    
    
    
    lipid_name = "PG(22:1(5Z)/12:0)";
    lipid = swiss_lipids_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->lipid->info->level == FULL_STRUCTURE);
    assert(lipid->get_lipid_string() == "PG 22:1(5Z)/12:0");
    delete lipid;
    
    
    lipid_name = "PG(22:1/12:0)";
    lipid = swiss_lipids_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->lipid->info->level == SN_POSITION);
    assert(lipid->get_lipid_string() == "PG 22:1/12:0");
    delete lipid;
    
    
    lipid_name = "PG(22:1_12:0)";
    lipid = swiss_lipids_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->lipid->info->level == MOLECULAR_SPECIES);
    assert(lipid->get_lipid_string() == "PG 22:1_12:0");
    delete lipid;
    
    lipid_name = "LPG(O-22:1)";
    lipid = swiss_lipids_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->lipid->info->level == SPECIES);
    assert(lipid->get_lipid_string() == "LPG O-22:1");
    delete lipid;
    
    
    lipid_name = "LPC O-16:1a";
    lipid = lipid_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "LPC O-16:1");
    delete lipid;
    
    
    
    lipid_name = "2-tetracosyl-3-hydroxy-pentatriaconta-18Z-enoic acid";
    lipid = lipid_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "FA 35:1(18Z);2(24:0);3OH");
    assert(lipid_parser.lastSuccessfulParser);
    assert(lipid_parser.lastSuccessfulParser->grammar_name == "FattyAcids");
    delete lipid;
    
    
    
    
    
    lipid_name = "SM(1) 18:1(10Z);3OH/14:0";
    lipid = lipid_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "SM(1) 18:1(10Z);3OH/14:0");
    delete lipid;
    
    
    
    lipid_name = "LPE O-16:4p";
    lipid = lipid_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "LPE P-16:3");
    delete lipid;
    
    
    lipid_name = "LPE O-16:2p";
    lipid = lipid_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string(CLASS) == "LPE");
    delete lipid;
    
    
    lipid_name = "fail";
    try {
        lipid = lipid_parser.parse(lipid_name);
        assert (false);
    }
    catch(LipidException &e){ }
    
    
    lipid_name = "LP 19:1p";
    try {
        lipid = lipid_parser.parse(lipid_name);
        assert (false);
    }
    catch(LipidException &e){
    }
    
    
    
    lipid_name = "PA 19:2(12E)/12:0";
    try {
        lipid = lipid_parser.parse(lipid_name);
        string l = lipid->get_lipid_string();
        assert (false);
    }
    catch(LipidException &e){
    }
    
    
    
    // test SwissLipids parser
    lipid_name = "TG(O-16:0/18:3(6Z,9Z,12Z)/18:1(11Z))";
    lipid = swiss_lipids_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "TG O-16:0/18:3(6Z,9Z,12Z)/18:1(11Z)");
    delete lipid;
    
    lipid_name = "PIP2[4,5](21:0/24:4(9Z,12Z,15Z,18Z))";
    lipid = swiss_lipids_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "PIP2(4',5') 21:0/24:4(9Z,12Z,15Z,18Z)");
    delete lipid;
    
    lipid_name = "GalGb3Cer(d18:0/14:0)";
    lipid = swiss_lipids_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "GalGb3Cer 18:0;OH/14:0");
    delete lipid;
    
    lipid_name = "PI(34:5(19Z,22Z,25Z,28Z,31Z)/18:1(6Z))";
    lipid = swiss_lipids_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "PI 34:5(19Z,22Z,25Z,28Z,31Z)/18:1(6Z)");
    delete lipid;
    
    lipid_name = "PIP(22:6/20:4)";
    lipid = swiss_lipids_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "PIP 22:6/20:4");
    delete lipid;
    
    lipid_name = "NAPE (12:0/30:4(15Z,18Z,21Z,24Z)/12:0)";
    lipid = swiss_lipids_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "PE-N(FA 12:0) 12:0/30:4(15Z,18Z,21Z,24Z)");
    delete lipid;
    
    lipid_name = "NAPE (12:0/0:0/12:0)";
    lipid = swiss_lipids_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "LPE-N(FA 12:0) 12:0/0:0");
    delete lipid;
    
    /*
    vector<FattyAcid*>* fa_list = new vector<FattyAcid*>;
    LipidSpecies* ls = new LipidStructuralSubspecies("PA", fa_list);
    try {
        cout << ls->get_lipid_string() << endl;
    }
    catch (LipidException &e){
        cout << "Exception:" << endl;
        cout << e.what() << endl;
    }
    
    //delete ls;
    //delete fa_list;
    */
    
    // test lipid parser
    lipid_name = "PE 16:1-12:0";
    lipid = lipid_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "PE 16:1_12:0");
    delete lipid;
    
    
    lipid_name = "PE O-16:1p/12:0";
    lipid = lipid_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "PE P-16:0/12:0");
    // check that all FAs have been initialized properly
    assert(lipid->lipid->get_fa_list().size() == 2);
    int faCnt = 1;
    for(FattyAcid* fa:lipid->lipid->get_fa_list()) {
        // cout << fa->name << " (lcb: " << (fa->lcb?"true":"false") << "): pos=" << fa->position << ", #carbon=" << fa->num_carbon << ", #hydroxyl=" << fa->num_hydroxyl << ", #double-bonds=" << fa->double_bonds->get_num() << endl;
        assert (fa->position == faCnt);
        switch(faCnt) {
        case 1:
        assert(fa->lipid_FA_bond_type != LCB_EXCEPTION && fa->lipid_FA_bond_type != LCB_REGULAR);
        assert (fa->num_carbon == 16);
        num_hydroxyl = contains_val_p(fa->functional_groups, "OH") ? fa->functional_groups->at("OH").size() : 0;
        assert (num_hydroxyl == 0);
        assert (FattyAcid::get_prefix(fa->lipid_FA_bond_type) == "P-");
        assert (fa->double_bonds->get_num() == 0);
        break;
        case 2:
        assert(fa->lipid_FA_bond_type != LCB_EXCEPTION && fa->lipid_FA_bond_type != LCB_REGULAR);
        assert (fa->num_carbon == 12);
        num_hydroxyl = contains_val_p(fa->functional_groups, "OH") ? fa->functional_groups->at("OH").size() : 0;
        assert (num_hydroxyl == 0);
        assert (FattyAcid::get_prefix(fa->lipid_FA_bond_type) == "");
        assert (fa->double_bonds->get_num() == 0);
        break;
        }
        ++faCnt;
    }
    delete lipid;
    
    
    lipid_name = "PAT16 16:1/12:0/14:1/8:0";
    lipid = goslin_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "PAT16 16:1/12:0/14:1/8:0");
    delete lipid;
    
    lipid_name = "SLBPA 16:1/12:0/14:1";
    lipid = lipid_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "SLBPA 16:1_12:0_14:1");
    delete lipid;
    
    lipid_name = "MLCL 16:1/12:0/14:1";
    lipid = lipid_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "LCL 16:1_12:0_14:1");
    delete lipid;
    
    
    lipid_name = "DLCL 14:1/8:0";
    lipid = lipid_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "DLCL 14:1_8:0");
    delete lipid;
    
    
    
    lipid_name = "PIP[3'] 17:0/20:4(5Z,8Z,11Z,14Z)";
    lipid = lipid_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "PIP(3') 17:0/20:4(5Z,8Z,11Z,14Z)");
    delete lipid;
    
    lipid_name = "PIMIP 12:0/14:1";
    lipid = lipid_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "PIMIP 12:0/14:1");
    delete lipid;
    
    lipid_name = "LCDPDAG 14:1";
    lipid = lipid_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "LCDPDAG 14:1");
    delete lipid;
    
    
    lipid_name = "CPA 8:0";
    lipid = lipid_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "CPA 8:0");
    delete lipid;
    
    
    lipid_name = "LPIN 20:4(5Z,8Z,11Z,14Z)";
    lipid = lipid_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "LPIN 20:4");
    delete lipid;

    lipid = lipid_parser.parse("Cer 36:1;2");
    int ohCount = lipid->lipid->info->get_total_functional_group_count("OH");
    assert(2==ohCount);
    int asdCount = lipid->lipid->info->get_total_functional_group_count("ASD");
    assert(0==asdCount);
    lipid = lipid_parser.parse("Cer d36:1");
    ohCount = lipid->lipid->info->get_total_functional_group_count("OH");
    assert(2==ohCount);
    lipid = lipid_parser.parse("Cer 18:1;2/18:0");
    ohCount = lipid->lipid->info->get_total_functional_group_count("OH");
    assert(2==ohCount);
    lipid = lipid_parser.parse("Cer d18:1/18:0");
    ohCount = lipid->lipid->info->get_total_functional_group_count("OH");
    assert(2==ohCount);
    lipid = lipid_parser.parse("Cer 18:1;(OH)2/18:0");
    ohCount = lipid->lipid->info->get_total_functional_group_count("OH");
    assert(2==ohCount);
    
    // testing lipid maps parser
    vector< vector<string> > lmp_data{{"PA(16:1/12:0)", "PA 16:1/12:0"},
                        {"PA(4:0/12:0)", "PA 4:0/12:0"},
                        {"PC(O-14:0/0:0)", "LPC O-14:0/0:0"},
                        {"SQMG(16:1(11Z)/0:0)", "SQMG 16:1(11Z)/0:0"},
                        {"TG(13:0/22:3(10Z,13Z,16Z)/22:5(7Z,10Z,13Z,16Z,19Z))[iso6]", "TG 13:0/22:3(10Z,13Z,16Z)/22:5(7Z,10Z,13Z,16Z,19Z)"},
                        {"13R-HODE", "13R-HODE"},
                        {"CL(1'-[20:0/20:0],3'-[20:4(5Z,8Z,11Z,14Z)/18:2(9Z,12Z)])", "CL 20:0/20:0/20:4(5Z,8Z,11Z,14Z)/18:2(9Z,12Z)"},
                        {"PA(P-20:0/18:3(6Z,9Z,12Z))", "PA P-20:0/18:3(6Z,9Z,12Z)"},
                        {"M(IP)2C(t18:0/20:0(2OH))", "M(IP)2C(1) 18:0;3OH,4OH/20:0;2OH"},
                        {"Cer(d16:2(4E,6E)/22:0(2OH))", "Cer 16:2(4E,6E);1OH,3OH/22:0;2OH"},
                        {"MG(18:1(11E)/0:0/0:0)[rac]", "MG 18:1(11E)/0:0/0:0"},
                        {"PAT18(24:1(2E)(2Me,4Me[S],6Me[S])/25:1(2E)(2Me,4Me[S],6Me[S])/26:1(2E)(2Me,4Me[S],6Me[S])/24:1(2E)(2Me,4Me[S],6Me[S]))", "PAT18 24:1(2E);2Me,4Me,6Me/25:1(2E);2Me,4Me,6Me/26:1(2E);2Me,4Me,6Me/24:1(2E);2Me,4Me,6Me"},
                        {"(3'-sulfo)Galbeta-Cer(d18:1/20:0)", "SHexCer 18:1;O2/20:0"},
                        {"GlcCer(d15:2(4E,6E)/22:0(2OH))", "GlcCer(1) 15:2(4E,6E);3OH/22:0;2OH"}};
    
    for (auto &lmp : lmp_data){
        lipid = lipid_maps_parser.parse(lmp.at(0));
        assert (lipid);
        assert(lipid->get_lipid_string() == lmp.at(1));
        delete lipid;
    }

    
    
    // testing lyso lipids
    lipid = lipid_parser.parse("LPA O-16:1a");
    assert(lipid);
    assert(lipid->get_lipid_string() == "LPA O-16:1");
    delete lipid;
    
    lipid = lipid_parser.parse("LPC O-16:1a");
    assert(lipid);
    assert(lipid->get_lipid_string() == "LPC O-16:1");
    delete lipid;
    
    lipid = lipid_parser.parse("LPE O-16:1p");
    assert(lipid);
    assert(lipid->get_lipid_string() == "LPE P-16:0");
    delete lipid;
    
    lipid = lipid_parser.parse("LPS O-16:1p");
    assert(lipid);
    assert(lipid->get_lipid_string() == "LPS P-16:0");
    delete lipid;
    
    lipid = lipid_parser.parse("LCB 18:1;2");
    assert(lipid);
    assert(lipid->get_lipid_string() == "SPB 18:1;O2");
    delete lipid;
    
    
    try {
        lipid = lipid_parser.parse("LPE O-16:1p/12:0");
        assert(false);
    }
    catch (LipidException &e){}
    
    
    
    
    
    // testing some failing situations
    for (string test_lipid_name : {"PE 18:3:1-16:2", "CL 18:1/18:1-18:1-18:1", "MLL 18:1-18:1-18:1", "M(IP )2C 18:0;3/26:0;1", "PET 16:0-18:1", "CDPDAG 18:1-18:1-18:1", "BeMP 18:1-18:1", "PE O18:0a-22:6", "GB3 18:1;2-24:1;12", "DGDGDG 16:0-16:1", "DG  16:0-16:1", "LPC 14:1a", "LPA14:1a", "LCB 17:1;2/6:0", "LP 19:1p", "xx2-34", "Hex2Cer 18:1;2/12:0/2:0", "Cer 18:1;2\\16:0", "MAG 16:0-12:0", "PC 21:0+22:6", "PE 16:2-18:3;1 - GP(153e)", "PE 16:2-18:3;1 - FA 18:3;1(+HO)", "CL 14:1-16:1-18:1-20:1 - FA 18:1(-H)-FA 20:1(+HO) [-2]", "MLCL 14:0-16:1-18:2 - FA 14:0(+O) +FA 16:1(+O) + HG(GP)", "SHexCer 18:0;3/26:0;1 - LCB 18:1;3(-CH3O)", "PEt 16:0-18:1 - -(-NH4)", "CDPDAG 18:1-18:1 - HG(C,2,08)", "BMP 18:1-18:1 - FA 18:2(+O)+ HG(GP,155)", "PE O-18:0a/22:6-GP(136)", "Hex3Cer 18:1;2/24:1 - LCB 24:1(-CH3O2)", "DGDG 16:0-16:1 - HG(MGDG,379)", "LCB 17:1;2 -LCB(60)", "Hex2Cer 18:1;2/12:0 -  -(+CH3COO)", "Cer 18:1;2/16:0 - FA 18:1;2(-C2H3NO)", "MAG 16:0 - FA 17:0", "PIP2 21:0-22:6 - -(2 *H3PO4+H2O,214)", "LPC(+[13]C8) O-12:1a - (C5H13NO,104)(+[13]O)", "LPC(+[18]O7) O-12:1a - (C5H13NO,104)"}){
        try {
            lipid = lipid_maps_parser.parse(test_lipid_name);
            assert(false);
        }
        catch (LipidException &e){}
    }
    
    
    
    
    
    // testing dedicated lipid species
    lipid_name = "Cer 28:1;2";
    lipid = goslin_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "Cer 28:1;O2");
    delete lipid;
    
    lipid_name = "DAG 38:1";
    lipid = goslin_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string() == "DG 38:1");
    delete lipid;
    
    
    
    
    // check if goslin parser fails correctly on parsing lipid name with fragment
    lipid_name = "PE 16:1-12:0 - -(H20)";
    try {
        lipid = goslin_parser.parse(lipid_name);
        assert(false);
    }
    catch (LipidException &e){}
    
    
    
    
        
    // test the down leveling of lipid names
    // glycerophospholipid;
    lipid_name = "PE 16:1(2E)/12:0";
    lipid = goslin_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string(FULL_STRUCTURE) == "PE 16:1(2E)/12:0");
    assert(lipid->get_lipid_string(STRUCTURE_DEFINED) == "PE 16:1(2)/12:0");
    assert(lipid->get_lipid_string(SN_POSITION) == "PE 16:1/12:0");
    assert(lipid->get_lipid_string(MOLECULAR_SPECIES) == "PE 16:1_12:0");
    assert(lipid->get_lipid_string(SPECIES) == "PE 28:1");
    assert(lipid->get_lipid_string(CLASS) == "PE");
    assert(lipid->get_lipid_string(CATEGORY) == "GP");
    delete lipid;
    
    
    // sphingolipid;
    lipid_name = "Cer 16:1(14E);2/12:0";
    lipid = goslin_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string(STRUCTURE_DEFINED) == "Cer 16:1(14);(OH)2/12:0");
    assert(lipid->get_lipid_string(SN_POSITION) == "Cer 16:1;O2/12:0");
    assert(lipid->get_lipid_string(MOLECULAR_SPECIES) == "Cer 16:1;O2/12:0");
    assert(lipid->get_lipid_string(SPECIES) == "Cer 28:1;O2");
    assert(lipid->get_lipid_string(CLASS) == "Cer");
    assert(lipid->get_lipid_string(CATEGORY) == "SP");
    // check that all FAs have been initialized properly
    assert(lipid->lipid->get_fa_list().size() == 2);
    faCnt = 1;
    for(FattyAcid* fa : lipid->lipid->get_fa_list()) {
        switch(faCnt) {
        case 1:
            assert(fa->lipid_FA_bond_type == LCB_EXCEPTION || fa->lipid_FA_bond_type == LCB_REGULAR);
            assert (fa->num_carbon == 16);
            num_hydroxyl = contains_val_p(fa->functional_groups, "OH") ? fa->functional_groups->at("OH").size() : 0;
            assert (num_hydroxyl == 1);
            assert (FattyAcid::get_prefix(fa->lipid_FA_bond_type) == "");
            assert (fa->double_bonds->get_num() == 1);
            break;
        case 2:
            assert(fa->lipid_FA_bond_type != LCB_EXCEPTION && fa->lipid_FA_bond_type != LCB_REGULAR);
            assert (fa->num_carbon == 12);
            num_hydroxyl = contains_val_p(fa->functional_groups, "OH") ? fa->functional_groups->at("OH").size() : 0;
            assert (num_hydroxyl == 0);
            assert (FattyAcid::get_prefix(fa->lipid_FA_bond_type) == "");
            assert (fa->double_bonds->get_num() == 0);
            break;
        }
        ++faCnt;
    }
    delete lipid;
    
    // glycerolipid;
    
    lipid_name = "TAG 16:1(5E)/18:0/20:2(3Z,6Z)";
    lipid = goslin_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string(STRUCTURE_DEFINED) == "TG 16:1(5)/18:0/20:2(3,6)");
    assert(lipid->get_lipid_string(SN_POSITION) == "TG 16:1/18:0/20:2");
    assert(lipid->get_lipid_string(MOLECULAR_SPECIES) == "TG 16:1_18:0_20:2");
    assert(lipid->get_lipid_string(SPECIES) == "TG 54:3");
    assert(lipid->get_lipid_string(CLASS) == "TG");
    assert(lipid->get_lipid_string(CATEGORY) == "GL");
    assert(lipid->lipid);
    
    
    // try to retrieve LipidSpeciesInfo for summary information
    LipidSpeciesInfo *lsi = lipid->lipid->info;
    assert(lsi->lipid_FA_bond_type != LCB_EXCEPTION && lsi->lipid_FA_bond_type != LCB_REGULAR);
    assert (lsi->level == FULL_STRUCTURE);
    assert (lsi->lipid_FA_bond_type == ESTER);
    assert (lsi->num_carbon == 54);
    assert (lsi->double_bonds->get_num() == 3);
    num_hydroxyl = contains_val_p(lsi->functional_groups, "OH") ? lsi->functional_groups->at("OH").size() : 0;
    assert (num_hydroxyl == 0);
    assert (lsi->position == 0);
    // check that all FAs have been initialized properly
    assert(lipid->lipid->get_fa_list().size() == 3);
    faCnt = 1;
    for(FattyAcid* fa : lipid->lipid->get_fa_list()) {
        // cout << fa->name << " (lcb: " << (fa->lcb?"true":"false") << "): pos=" << fa->position << ", #carbon=" << fa->num_carbon << ", #hydroxyl=" << fa->num_hydroxyl << ", #double-bonds=" << fa->double_bonds->get_num() << endl;
        assert (fa->position == faCnt);
        switch(faCnt) {
        case 1:
            assert(fa->lipid_FA_bond_type != LCB_EXCEPTION && fa->lipid_FA_bond_type != LCB_REGULAR);
            assert (fa->num_carbon == 16);
            num_hydroxyl = contains_val_p(fa->functional_groups, "OH") ? fa->functional_groups->at("OH").size() : 0;
            assert (num_hydroxyl == 0);
            assert (FattyAcid::get_prefix(fa->lipid_FA_bond_type) == "");
            assert (fa->double_bonds->get_num() == 1);
            break;
        case 2:
            assert(fa->lipid_FA_bond_type != LCB_EXCEPTION && fa->lipid_FA_bond_type != LCB_REGULAR);
            assert (fa->num_carbon == 18);
            num_hydroxyl = contains_val_p(fa->functional_groups, "OH") ? fa->functional_groups->at("OH").size() : 0;
            assert (num_hydroxyl == 0);
            assert (FattyAcid::get_prefix(fa->lipid_FA_bond_type) == "");
            assert (fa->double_bonds->get_num() == 0);
            break;
        case 3:
            assert(fa->lipid_FA_bond_type != LCB_EXCEPTION && fa->lipid_FA_bond_type != LCB_REGULAR);
            assert (fa->num_carbon == 20);
            num_hydroxyl = contains_val_p(fa->functional_groups, "OH") ? fa->functional_groups->at("OH").size() : 0;
            assert (num_hydroxyl == 0);
            assert (FattyAcid::get_prefix(fa->lipid_FA_bond_type) == "");
            assert (fa->double_bonds->get_num() == 2);
            break;
        }
        ++faCnt;
    }
    delete lipid;
    
    lipid_name = "TAG 16:1(8E)/12:0/20:2(4Z,8Z)";
    lipid = goslin_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string(STRUCTURE_DEFINED) == "TG 16:1(8)/12:0/20:2(4,8)");
    assert(lipid->get_lipid_string(SN_POSITION) == "TG 16:1/12:0/20:2");
    assert(lipid->get_lipid_string(MOLECULAR_SPECIES) == "TG 16:1_12:0_20:2");
    assert(lipid->get_lipid_string(SPECIES) == "TG 48:3");
    assert(lipid->get_lipid_string(CLASS) == "TG");
    assert(lipid->get_lipid_string(CATEGORY) == "GL");
    assert(lipid->lipid);
    
    
    // try to retrieve LipidSpeciesInfo for summary information
    lsi = lipid->lipid->info;
    assert(lsi->lipid_FA_bond_type != LCB_EXCEPTION && lsi->lipid_FA_bond_type != LCB_REGULAR);
    assert (lsi->level == FULL_STRUCTURE);
    assert (lsi->lipid_FA_bond_type == ESTER);
    assert (lsi->num_carbon == 48);
    assert (lsi->double_bonds->get_num() == 3);
    num_hydroxyl = contains_val_p(lsi->functional_groups, "OH") ? lsi->functional_groups->at("OH").size() : 0;
    assert (num_hydroxyl == 0);
    assert (lsi->position == 0);
    // check that all FAs have been initialized properly
    assert(lipid->lipid->get_fa_list().size() == 3);
    faCnt = 1;
    for(FattyAcid* fa:lipid->lipid->get_fa_list()) {
        // cout << fa->name << " (lcb: " << (fa->lcb?"true":"false") << "): pos=" << fa->position << ", #carbon=" << fa->num_carbon << ", #hydroxyl=" << fa->num_hydroxyl << ", #double-bonds=" << fa->double_bonds->get_num() << endl;
        assert (fa->position == faCnt);
        switch(faCnt) {
        case 1:
            assert(fa->lipid_FA_bond_type != LCB_EXCEPTION && fa->lipid_FA_bond_type != LCB_REGULAR);
            assert (fa->num_carbon == 16);
            num_hydroxyl = contains_val_p(fa->functional_groups, "OH") ? fa->functional_groups->at("OH").size() : 0;
            assert (num_hydroxyl == 0);
            assert (FattyAcid::get_prefix(fa->lipid_FA_bond_type) == "");
            assert (fa->double_bonds->get_num() == 1);
            break;
        case 2:
            assert(fa->lipid_FA_bond_type != LCB_EXCEPTION && fa->lipid_FA_bond_type != LCB_REGULAR);
            assert (fa->num_carbon == 12);
            num_hydroxyl = contains_val_p(fa->functional_groups, "OH") ? fa->functional_groups->at("OH").size() : 0;
            assert (num_hydroxyl == 0);
            assert (FattyAcid::get_prefix(fa->lipid_FA_bond_type) == "");
            assert (fa->double_bonds->get_num() == 0);
            break;
        case 3:
            assert(fa->lipid_FA_bond_type != LCB_EXCEPTION && fa->lipid_FA_bond_type != LCB_REGULAR);
            assert (fa->num_carbon == 20);
            num_hydroxyl = contains_val_p(fa->functional_groups, "OH") ? fa->functional_groups->at("OH").size() : 0;
            assert (num_hydroxyl == 0);
            assert (FattyAcid::get_prefix(fa->lipid_FA_bond_type) == "");
            assert (fa->double_bonds->get_num() == 2);
            break;
        }
        ++faCnt;
    }
    delete lipid;
    
    
    
    // sterol;
    lipid_name = "ChE 16:1(12E)";
    lipid = goslin_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string(STRUCTURE_DEFINED) == "SE 27:1/16:1(12)");
    assert(lipid->get_lipid_string(SN_POSITION) == "SE 27:1/16:1");
    assert(lipid->get_lipid_string(MOLECULAR_SPECIES) == "SE 27:1/16:1");
    assert(lipid->get_lipid_string(SPECIES) == "SE 27:1/16:1");
    assert(lipid->get_lipid_string(CLASS) == "SE 27:1");
    assert(lipid->get_lipid_string(CATEGORY) == "ST");
    delete lipid;
    
    
    
    
    // sterol;
    lipid_name = "ChE 16:1(12E)";
    lipid = goslin_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string(STRUCTURE_DEFINED) == "SE 27:1/16:1(12)");
    assert(lipid->get_lipid_string(SN_POSITION) == "SE 27:1/16:1");
    assert(lipid->get_lipid_string(MOLECULAR_SPECIES) == "SE 27:1/16:1");
    assert(lipid->get_lipid_string(SPECIES) == "SE 27:1/16:1");
    assert(lipid->get_lipid_string(CLASS) == "SE 27:1");
    assert(lipid->get_lipid_string(CATEGORY) == "ST");
    delete lipid;
    

    
    // PC;
    lipid_name = "PC O-16:1(13E)a/12:0";
    lipid = goslin_parser.parse(lipid_name);
    assert (lipid);
    assert(lipid->get_lipid_string(STRUCTURE_DEFINED) == "PC O-16:1(13)/12:0");
    assert(lipid->get_lipid_string(SN_POSITION) == "PC O-16:1/12:0");
    assert(lipid->get_lipid_string(MOLECULAR_SPECIES) == "PC O-16:1_12:0");
    assert(lipid->get_lipid_string(SPECIES) == "PC O-28:1");
    assert(lipid->get_lipid_string(CLASS) == "PC");
    assert(lipid->get_lipid_string(CATEGORY) == "GP");
    delete lipid;
    
    
    // testing adducts
    lipid_name = "PE 16:1/12:0[M+H]1+";
    lipid = goslin_parser.parse(lipid_name);
    assert(lipid);
    assert(lipid->get_lipid_string() == "PE 16:1/12:0[M+H]1+");
    delete lipid;

    
    

    // test several more lipid names
    vector<string> lipidnames;
    ifstream infile("data/goslin/testfiles/lipidnames.csv");
    assert(infile.good());
    string line;
    while (getline(infile, line)){
        line = strip(line, ' ');
        if (line.length() < 2) continue;
        if (line[0] == '#') continue;
        vector<string> *v = split_string(line, ',', '"');
        lipidnames.push_back(strip(v->at(0), '"'));
        delete v;
    }
    infile.close();
    
    for (auto &test_lipid_name : lipidnames){
        lipid = lipid_parser.parse(test_lipid_name);
        if (lipid == NULL){
            cout << "Fail: '" << test_lipid_name << "'" << endl;
        }
        assert(lipid);
        delete lipid;
    }
        

    cout << "All tests passed without any problem" << endl;
    
    return EXIT_SUCCESS;
}
