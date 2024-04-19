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


#include "cppgoslin/parser/FattyAcidParserEventHandler.h"

#define reg(x, y) BaseParserEventHandler<LipidAdduct*>::registered_events->insert({x, bind(&FattyAcidParserEventHandler::y, this, placeholders::_1)})
#define FA_I ("fa" + std::to_string(fatty_acyl_stack.size()))

FattyAcidParserEventHandler::FattyAcidParserEventHandler() : BaseParserEventHandler<LipidAdduct*>() {
    
        
    reg("lipid_pre_event", reset_lipid);
    reg("lipid_post_event", build_lipid);
    reg("fatty_acid_post_event", set_fatty_acid);
    reg("fatty_acid_recursion_post_event", set_fatty_acid);
    
    reg("acid_single_type_pre_event", set_fatty_acyl_type);
    reg("ol_ending_pre_event", set_fatty_acyl_type);
    reg("double_bond_position_pre_event", set_double_bond_information);
    reg("double_bond_position_post_event", add_double_bond_information);
    reg("db_number_post_event", set_double_bond_position);
    reg("cistrans_post_event", set_cistrans);
    reg("acid_type_double_post_event", check_db);
    reg("db_length_pre_event", open_db_length);
    reg("db_length_post_event", close_db_length);
    
    // lengths
    reg("functional_length_pre_event", reset_length);
    reg("fatty_length_pre_event", reset_length);
    reg("functional_length_post_event", set_functional_length);
    reg("fatty_length_post_event", set_fatty_length);
    
    // numbers
    reg("notation_specials_pre_event", special_number);
    reg("notation_last_digit_pre_event", last_number);
    reg("notation_second_digit_pre_event", second_number);
    
    // furan
    reg("tetrahydrofuran_pre_event", set_tetrahydrofuran);
    reg("furan_pre_event", set_furan);
    
    // functional groups
    reg("functional_group_pre_event", set_functional_group);
    reg("functional_group_post_event", add_functional_group);
    reg("functional_pos_pre_event", set_functional_pos);
    reg("functional_position_pre_event", set_functional_position);
    reg("functional_group_type_pre_event", set_functional_type);
    
    // cyclo / epoxy
    reg("cyclo_position_pre_event", set_functional_group);
    reg("cyclo_position_post_event", rearrange_cycle);
    reg("epoxy_pre_event", set_functional_group);
    reg("epoxy_post_event", add_epoxy);
    reg("cycle_pre_event", set_cycle);
    reg("methylene_post_event", set_methylene);

    // dioic
    reg("dioic_pre_event", set_functional_group);
    reg("dioic_post_event", set_dioic);
    reg("dioic_acid_pre_event", set_fatty_acyl_type);
    reg("dial_post_event", set_dial);

    
    // prosta
    reg("prosta_pre_event", set_prosta);
    reg("prosta_post_event", add_cyclo);
    reg("reduction_pre_event", set_functional_group);
    reg("reduction_post_event", reduction);
    reg("homo_post_event", homo);

    
    // recursion
    reg("recursion_description_pre_event", set_recursion);
    reg("recursion_description_post_event", add_recursion);
    reg("recursion_pos_pre_event", set_recursion_pos);
    reg("yl_ending_pre_event", set_yl_ending);
    reg("acetic_acid_post_event", set_acetic_acid);
    reg("acetic_recursion_pre_event", set_recursion);
    reg("acetic_recursion_post_event", add_recursion);
    reg("hydroxyl_number_pre_event", add_hydroxyl);
    reg("ol_pre_event", setup_hydroxyl);
    reg("ol_post_event", add_hydroxyls);
    reg("ol_pos_post_event", set_yl_ending);
    
    
    // wax esters
    reg("wax_ester_pre_event", set_recursion);
    reg("wax_ester_post_event", add_wax_ester);
    reg("ate_post_event", set_ate);
    reg("isoprop_post_event", set_iso);
    reg("isobut_post_event", set_iso);
    
    // CoA
    reg("coa_post_event", set_coa);
    reg("methyl_pre_event", set_methyl);
    
    // CAR
    reg("car_pre_event", set_car);
    reg("car_post_event", add_car);
    
    // amine
    reg("ethanolamine_post_event", add_ethanolamine);
    reg("amine_n_pre_event", set_recursion);
    reg("amine_n_post_event", add_amine);
    reg("amine_post_event", add_amine_name);
    
    // functional group position summary
    reg("fg_pos_summary_pre_event", set_functional_group);
    reg("fg_pos_summary_post_event", add_summary);
    reg("func_stereo_pre_event", add_func_stereo);
    
    debug = "";
}



const map<string, int> FattyAcidParserEventHandler::last_numbers{{"un", 1}, {"hen", 1}, {"do", 2}, {"di", 2}, {"tri", 3}, {"buta", 4}, {"but", 4}, {"tetra", 4}, {"penta", 5}, {"pent", 5}, {"hexa", 6}, {"hex", 6}, {"hepta", 7}, {"hept", 7}, {"octa", 8}, {"oct", 8}, {"nona", 9}, {"non", 9}};


const map<string, int> FattyAcidParserEventHandler::second_numbers {{"deca", 10}, {"dec", 10}, {"eicosa", 20}, {"eicos", 20 }, {"cosa", 20}, {"cos", 20}, {"triaconta", 30}, {"triacont", 30}, {"tetraconta", 40}, {"tetracont", 40}, {"pentaconta", 50}, {"pentacont", 50}, {"hexaconta", 60}, {"hexacont", 60}, {"heptaconta", 70}, {"heptacont", 70}, {"octaconta", 80}, {"octacont", 80}, {"nonaconta", 90}, {"nonacont", 90}};

const map<string, string> FattyAcidParserEventHandler::func_groups {{"keto", "oxo"}, {"ethyl", "Et"}, {"hydroxy", "OH"}, {"phospho", "Ph"}, {"oxo", "oxo"}, {"bromo", "Br"}, {"methyl", "Me"}, {"hydroperoxy", "OOH"}, {"homo", ""}, {"Epoxy", "Ep"}, {"fluro", "F"}, {"fluoro", "F"}, {"chloro", "Cl"}, {"methylene", "My"}, {"sulfooxy", "Su"}, {"amino", "NH2"}, {"sulfanyl", "SH"}, {"methoxy", "OMe"}, {"iodo", "I"}, {"cyano", "CN"}, {"nitro", "NO2"}, {"OH", "OH"}, {"thio", "SH"}, {"mercapto", "SH"}, {"carboxy", "COOH"}, {"acetoxy", "Ac"}, {"cysteinyl", "Cys"}, {"phenyl", "Phe"}, {"s-glutathionyl", "SGlu"}, {"s-cysteinyl", "SCys"}, {"butylperoxy", "BOO"}, {"dimethylarsinoyl", "MMAs"}, {"methylsulfanyl", "SMe"}, {"imino", "NH"}, {"s-cysteinylglycinyl", "SCG"}};

const map<string, int> FattyAcidParserEventHandler::ate {{"formate", 1}, {"acetate", 2}, {"butyrate", 4}, {"propionate", 3}, {"valerate", 5}, {"isobutyrate", 4}};

const map<string, int> FattyAcidParserEventHandler::special_numbers {{"meth", 1}, {"etha", 2}, {"eth", 2}, {"propa", 3}, {"isoprop", 3}, {"prop", 3}, {"propi", 3}, {"propio", 3}, {"buta", 4}, {"but", 4}, {"butr", 4}, {"furan", 5}, {"valer", 5}, {"eicosa", 20}, {"eicos", 20}, {"icosa", 20}, {"icos", 20}, {"prosta", 20}, {"prost", 20}, {"prostan", 20}};


FattyAcidParserEventHandler::~FattyAcidParserEventHandler(){
}


void FattyAcidParserEventHandler::set_lipid_level(LipidLevel _level){
    level = min(level, _level);
}


void FattyAcidParserEventHandler::reset_lipid(TreeNode *node) {
    BaseParserEventHandler<LipidAdduct*>::content = 0;
    level = FULL_STRUCTURE;
    headgroup = "";
    fatty_acyl_stack.clear();
    fatty_acyl_stack.push_back(new FattyAcid("FA"));
    tmp.remove_all();
    tmp.set_dictionary("fa1", new GenericDictionary());
}


void FattyAcidParserEventHandler::build_lipid(TreeNode *node) {
    
    if (tmp.contains_key("cyclo_yl")) {
        tmp.set_list("fg_pos", new GenericList());
        tmp.get_list("fg_pos")->add_list(new GenericList());
        tmp.get_list("fg_pos")->get_list(0)->add_int(1);
        tmp.get_list("fg_pos")->get_list(0)->add_string("");
        tmp.get_list("fg_pos")->add_list(new GenericList());
        tmp.get_list("fg_pos")->get_list(1)->add_int(tmp.get_int("cyclo_len"));
        tmp.get_list("fg_pos")->get_list(1)->add_string("");
        add_cyclo(node);
        tmp.remove("cyclo_yl");
        tmp.remove("cyclo_len");
    }
    
            
    
    if (tmp.contains_key("post_adding")){
        FattyAcid *curr_fa = fatty_acyl_stack.back();
        int s = tmp.get_list("post_adding")->list.size();
        curr_fa->num_carbon += s;
        for (int i = 0; i < s; ++i){
            int pos = tmp.get_list("post_adding")->get_int(i);
            curr_fa->add_position(pos);
            DoubleBonds* db = new DoubleBonds(curr_fa->double_bonds->num_double_bonds);
            for (auto &kv : curr_fa->double_bonds->double_bond_positions){
                db->double_bond_positions.insert({kv.first + (kv.first >= pos), kv.second});
            }
            db->num_double_bonds = db->double_bond_positions.size();
            delete curr_fa->double_bonds;
            curr_fa->double_bonds = db;
        }
    }
    
    FattyAcid *curr_fa = fatty_acyl_stack.back();
    if (!curr_fa->double_bonds->double_bond_positions.empty()){
        int db_right = 0;
        for (auto &kv : curr_fa->double_bonds->double_bond_positions) db_right += kv.second.length() > 0;
        if (db_right != (int)curr_fa->double_bonds->double_bond_positions.size()){
            set_lipid_level(STRUCTURE_DEFINED);
        }
    }
    
    Headgroup *head_group = new Headgroup(headgroup);
    
    lipid = new LipidAdduct();
    
    switch(level){
        case COMPLETE_STRUCTURE:
            lipid->lipid = new LipidCompleteStructure(head_group, &fatty_acyl_stack);
            break;
        
        case FULL_STRUCTURE:
            lipid->lipid = new LipidFullStructure(head_group, &fatty_acyl_stack);
            break;
            
        case STRUCTURE_DEFINED:
            lipid->lipid = new LipidStructureDefined(head_group, &fatty_acyl_stack);
            break;
            
        case SN_POSITION:
            lipid->lipid = new LipidSnPosition(head_group, &fatty_acyl_stack);
            break;
            
        case MOLECULAR_SPECIES:
            lipid->lipid = new LipidMolecularSpecies(head_group, &fatty_acyl_stack);
            break;
            
        case SPECIES:
            lipid->lipid = new LipidSpecies(head_group, &fatty_acyl_stack);
            break;
            
        default:
            break;
    }
    BaseParserEventHandler<LipidAdduct*>::content = lipid;
}


void FattyAcidParserEventHandler::switch_position(FunctionalGroup* func_group, int switch_num){
    func_group->position = switch_num - func_group->position;
    for (auto &kv : *(func_group->functional_groups)){
        for (auto &fg : kv.second){
            switch_position(fg, switch_num);
        }
    }
}


void FattyAcidParserEventHandler::set_fatty_acid(TreeNode *node) {
    FattyAcid* curr_fa = fatty_acyl_stack.back();
    
    if (tmp.contains_key("length_pattern")) {
        
        string length_pattern = tmp.get_string("length_pattern");
        const int n_len = tmp.get_list("length_tokens")->list.size();
        int* num = new int[n_len];
        for (int i = 0; i < (int)tmp.get_list("length_tokens")->list.size(); ++i) num[i] = tmp.get_list("length_tokens")->get_int(i);
        
        int l = 0, d = 0;
        if (length_pattern == "L" || length_pattern == "S"){
            l += num[0];
        }
            
        else if (length_pattern == "LS"){
            l += num[0] + num[1];
        }
        
        else if (length_pattern == "LL" || length_pattern == "SL" || length_pattern == "SS"){
            l += num[0];
            d += num[1];
        }
            
        else if (length_pattern == "LSL" || length_pattern == "LSS"){
            l += num[0] + num[1];
            d += num[2];
        }
            
        else if (length_pattern == "LSLS"){
            l += num[0] + num[1];
            d += num[2] + num[3];
        }
            
        else if (length_pattern == "SLS"){
            l += num[0];
            d += num[1] + num[2];
        }
            
        else if (length_pattern.length() > 0 && length_pattern[0] == 'X'){
            l += num[0];
            for (int i = 1; i < (int)tmp.get_list("length_tokens")->list.size(); ++i) d += num[i];
        }
        
        else if (length_pattern == "LLS"){ // false
            throw RuntimeException("Cannot determine fatty acid and double bond length in '" + node->get_text() + "'");
        }
        
        curr_fa->num_carbon += l;
        if (curr_fa->double_bonds->double_bond_positions.size() == 0 && d > 0) curr_fa->double_bonds->num_double_bonds = d;
        delete []num;
    }
    
    
    if (contains_val_p(curr_fa->functional_groups, "noyloxy")){
        if (headgroup == "FA") headgroup = "FAHFA";
        
        while (!curr_fa->functional_groups->at("noyloxy").empty()){
            FattyAcid* fa = (FattyAcid*)curr_fa->functional_groups->at("noyloxy").back();
            curr_fa->functional_groups->at("noyloxy").pop_back();
        
            AcylAlkylGroup* acyl = new AcylAlkylGroup(fa);
            acyl->position = fa->position;
            
            if (uncontains_val_p(curr_fa->functional_groups, "acyl")) curr_fa->functional_groups->insert({"acyl", vector<FunctionalGroup*>()});
            curr_fa->functional_groups->at("acyl").push_back(acyl);
        }
        curr_fa->functional_groups->erase("noyloxy");
    }
        
    else if (contains_val_p(curr_fa->functional_groups, "nyloxy") || contains_val_p(curr_fa->functional_groups, "yloxy")){
        string yloxy = contains_val_p(curr_fa->functional_groups, "nyloxy") ? "nyloxy" : "yloxy";
        while (!curr_fa->functional_groups->at(yloxy).empty()){
            FattyAcid* fa = (FattyAcid*)curr_fa->functional_groups->at(yloxy).back();
            curr_fa->functional_groups->at(yloxy).pop_back();        
            
            AcylAlkylGroup* alkyl = new AcylAlkylGroup(fa, -1, 1, true);
            alkyl->position = fa->position;
            
            if (uncontains_val_p(curr_fa->functional_groups, "alkyl")) curr_fa->functional_groups->insert({"alkyl", vector<FunctionalGroup*>()});
            curr_fa->functional_groups->at("alkyl").push_back(alkyl);
        }
        curr_fa->functional_groups->erase(yloxy);
    }
            
    else {
        bool has_yl = false;
        for (auto &kv : *(curr_fa->functional_groups)){
            if (endswith(kv.first, "yl")){
                has_yl = true;
                break;
            }
        }
        if (has_yl){
            while (true){
                string yl = "";
                for (auto &kv : *(curr_fa->functional_groups)){
                    if (endswith(kv.first, "yl")){
                        yl = kv.first;
                        break;
                    }
                }
                if (yl.length() == 0) {
                    break;
                }
            
                while (!curr_fa->functional_groups->at(yl).empty()){
                    FattyAcid* fa = (FattyAcid*)curr_fa->functional_groups->at(yl).back();
                    curr_fa->functional_groups->at(yl).pop_back();
                    
                    if (tmp.contains_key("cyclo")){
                        int cyclo_len = curr_fa->num_carbon;
                        tmp.set_int("cyclo_len", cyclo_len);
                        if (fa->position != cyclo_len && !tmp.contains_key("furan")) {
                            switch_position(curr_fa, 2 + cyclo_len);
                        }
                        fa->shift_positions(cyclo_len);
                        if (tmp.contains_key("furan")) curr_fa->shift_positions(-1);
                        
                        for (auto &kv : *(fa->functional_groups)){
                            if (uncontains_val_p(curr_fa->functional_groups, kv.first)) {
                                curr_fa->functional_groups->insert({kv.first, vector<FunctionalGroup*>()});
                            }
                            for (auto &func_group : kv.second) curr_fa->functional_groups->at(kv.first).push_back(func_group);
                        }
                            
                        curr_fa->num_carbon = cyclo_len + fa->num_carbon;
                        
                        for (auto &kv : fa->double_bonds->double_bond_positions){
                            curr_fa->double_bonds->double_bond_positions.insert({kv.first + cyclo_len, kv.second});
                        }
                        curr_fa->double_bonds->num_double_bonds = curr_fa->double_bonds->double_bond_positions.size();
                        if (!tmp.contains_key("tetrahydrofuran") && tmp.contains_key("furan")){
                            curr_fa->double_bonds->num_double_bonds += 2;
                            if (uncontains_val(curr_fa->double_bonds->double_bond_positions, 1)) curr_fa->double_bonds->double_bond_positions.insert({1, "E"});
                            if (uncontains_val(curr_fa->double_bonds->double_bond_positions, 3)) curr_fa->double_bonds->double_bond_positions.insert({3, "E"});
                        }
                            
                        tmp.set_int("cyclo_yl", 1);
                    }
                    else {
                        // add carbon chains here here
                        // special chains: i.e. ethyl, methyl
                        string fg_name = "";
                        if (fa->double_bonds->get_num() == 0 && fa->functional_groups->empty()){
                            FunctionalGroup *fg = 0;
                            if (fa->num_carbon == 1){
                                fg_name = "Me";
                                fg = KnownFunctionalGroups::get_functional_group(fg_name);
                            }
                            else if (fa->num_carbon == 2){
                                fg_name = "Et";
                                fg = KnownFunctionalGroups::get_functional_group(fg_name);
                            }
                            if (fg_name.length() > 0){
                                fg->position = fa->position;
                                if (uncontains_val_p(curr_fa->functional_groups, fg_name)) curr_fa->functional_groups->insert({fg_name, vector<FunctionalGroup*>()});
                                curr_fa->functional_groups->at(fg_name).push_back(fg);
                            }
                        }
                        if (fg_name.length() == 0){
                            CarbonChain *cc = new CarbonChain(fa, fa->position);
                            if (uncontains_val_p(curr_fa->functional_groups, "cc")) curr_fa->functional_groups->insert({"cc", vector<FunctionalGroup*>()});
                            curr_fa->functional_groups->at("cc").push_back(cc);
                        }
                    }
                }
                if (tmp.contains_key("cyclo")) tmp.remove("cyclo");
                curr_fa->functional_groups->erase(yl);
            }
        }
    }
        
    if (contains_val_p(curr_fa->functional_groups, "cyclo")){
        FattyAcid *fa = (FattyAcid*)curr_fa->functional_groups->at("cyclo").front();
        curr_fa->functional_groups->erase("cyclo");
        if (!tmp.contains_key("cyclo_len")) tmp.set_int("cyclo_len", 5);
        int start_pos = curr_fa->num_carbon + 1;
        int end_pos = curr_fa->num_carbon + tmp.get_int("cyclo_len");
        fa->shift_positions(start_pos - 1);
        
        if (contains_val_p(curr_fa->functional_groups, "cy")){
            for (auto &cy : curr_fa->functional_groups->at("cy")){
                cy->shift_positions(start_pos - 1);
            }
        }
        for (auto &kv : *(fa->functional_groups)){
            if (uncontains_val_p(curr_fa->functional_groups, kv.first)){
                curr_fa->functional_groups->insert({kv.first, vector<FunctionalGroup*>()});
            }
            for (auto &func_group : kv.second){
                curr_fa->functional_groups->at(kv.first).push_back(func_group);
            }
        }
        
        
        
        
        for (auto &kv : fa->double_bonds->double_bond_positions) curr_fa->double_bonds->double_bond_positions.insert({kv.first + start_pos - 1, kv.second});
        curr_fa->double_bonds->num_double_bonds = curr_fa->double_bonds->double_bond_positions.size();
        
        if (!tmp.contains_key("tetrahydrofuran") and tmp.contains_key("furan")){
            curr_fa->double_bonds->num_double_bonds += 2;
            if (uncontains_val(curr_fa->double_bonds->double_bond_positions, 1 + curr_fa->num_carbon)) curr_fa->double_bonds->double_bond_positions.insert({1 + curr_fa->num_carbon, "E"});
            if (uncontains_val(curr_fa->double_bonds->double_bond_positions, 3 + curr_fa->num_carbon)) curr_fa->double_bonds->double_bond_positions.insert({3 + curr_fa->num_carbon, "E"});
        }
        curr_fa->num_carbon += fa->num_carbon;
                
        tmp.set_list("fg_pos", new GenericList());
        tmp.get_list("fg_pos")->add_list(new GenericList());
        tmp.get_list("fg_pos")->add_list(new GenericList());
        tmp.get_list("fg_pos")->get_list(0)->add_int(start_pos);
        tmp.get_list("fg_pos")->get_list(0)->add_string("");
        tmp.get_list("fg_pos")->get_list(1)->add_int(end_pos);
        tmp.get_list("fg_pos")->get_list(1)->add_string("");
        
        add_cyclo(node);
        
        if (tmp.contains_key("cyclo_len")) tmp.remove("cyclo_len");
        if (tmp.contains_key("cyclo")) tmp.remove("cyclo");
    }
        
    else if (tmp.contains_key("cyclo")){
        tmp.set_int("cyclo_yl", 1);
        tmp.set_int("cyclo_len", curr_fa->num_carbon);
        tmp.set_list("fg_pos", new GenericList());
        tmp.get_list("fg_pos")->add_list(new GenericList());
        tmp.get_list("fg_pos")->add_list(new GenericList());
        tmp.get_list("fg_pos")->get_list(0)->add_int(1);
        tmp.get_list("fg_pos")->get_list(0)->add_string("");
        tmp.get_list("fg_pos")->get_list(1)->add_int(curr_fa->num_carbon);
        tmp.get_list("fg_pos")->get_list(1)->add_string("");
        
        //add_cyclo(node);
        tmp.remove("cyclo");
    }
    
    tmp.set_string("length_pattern", "");
    tmp.set_list("length_tokens", new GenericList());
    tmp.set_int("add_lengths", 0);
}


const set<string> FattyAcidParserEventHandler::noic_set{"noic acid", "nic acid", "dioic_acid"};
const set<string> FattyAcidParserEventHandler::nal_set{"nal", "dial"};
const set<string> FattyAcidParserEventHandler::acetate_set{"acetate", "noate", "nate"};

void FattyAcidParserEventHandler::set_fatty_acyl_type(TreeNode *node) {
    string t = node->get_text();
    
    if (endswith(t, "ol")) headgroup = "FOH";
    else if (contains_val(noic_set, t)) headgroup = "FA";
    else if (contains_val(nal_set, t)) headgroup = "FAL";
    else if (contains_val(acetate_set, t)) headgroup = "WE";
    else if (t == "ne"){
        headgroup = "HC";
        fatty_acyl_stack.back()->lipid_FA_bond_type = ETHER;
    }
    else {
        headgroup = t;
    }
}


void FattyAcidParserEventHandler::set_double_bond_information(TreeNode *node) {
    tmp.get_dictionary(FA_I)->set_int("db_position", 0);
    tmp.get_dictionary(FA_I)->set_string("db_cistrans", "");
}


void FattyAcidParserEventHandler::add_double_bond_information(TreeNode *node) {    
    int pos = tmp.get_dictionary(FA_I)->get_int("db_position");
    string str_pos = std::to_string(pos);
    string cistrans = tmp.get_dictionary(FA_I)->get_string("db_cistrans");
    if (cistrans == "" && tmp.get_dictionary(FA_I)->contains_key("fg_pos_summary") && tmp.get_dictionary(FA_I)->get_dictionary("fg_pos_summary")->contains_key(str_pos)){
        cistrans = tmp.get_dictionary(FA_I)->get_dictionary("fg_pos_summary")->get_string(str_pos);
    }
    if (pos == 0) return;
    
    cistrans = to_upper(cistrans);
    
    tmp.get_dictionary(FA_I)->remove("db_position");
    tmp.get_dictionary(FA_I)->remove("db_cistrans");
    
    
    if (cistrans != "E" && cistrans != "Z") cistrans = "";
    if (uncontains_val(fatty_acyl_stack.back()->double_bonds->double_bond_positions, pos) || fatty_acyl_stack.back()->double_bonds->double_bond_positions.at(pos).length() == 0){
        if (uncontains_val(fatty_acyl_stack.back()->double_bonds->double_bond_positions, pos)){
            fatty_acyl_stack.back()->double_bonds->double_bond_positions.insert({pos, cistrans});
        }
        else {
            fatty_acyl_stack.back()->double_bonds->double_bond_positions.at(pos) = cistrans;
        }
        fatty_acyl_stack.back()->double_bonds->num_double_bonds = fatty_acyl_stack.back()->double_bonds->double_bond_positions.size();
    }
}


void FattyAcidParserEventHandler::set_double_bond_position(TreeNode *node) {
    int pos = atoi(node->get_text().c_str());
    int num_db = 0;
    if (tmp.contains_key("reduction")){
        GenericList *gl = tmp.get_list("reduction");
        int l = gl->list.size();
        for (int i = 0; i < l; ++i){
            num_db += gl->get_int(i) < pos;
        }
    }
    
    tmp.get_dictionary(FA_I)->set_int("db_position", pos - num_db);
}


void FattyAcidParserEventHandler::set_cistrans(TreeNode *node) {
    tmp.get_dictionary(FA_I)->set_string("db_cistrans", node->get_text());
}


void FattyAcidParserEventHandler::check_db(TreeNode *node) {
    FattyAcid* curr_fa = fatty_acyl_stack.back();
    if (tmp.get_dictionary(FA_I)->contains_key("fg_pos_summary")){
        for (auto &kv : tmp.get_dictionary(FA_I)->get_dictionary("fg_pos_summary")->dictionary){
            int k = atoi(kv.first.c_str());
            string v = tmp.get_dictionary(FA_I)->get_dictionary("fg_pos_summary")->get_string(kv.first);
            if (k > 0 && uncontains_val(curr_fa->double_bonds->double_bond_positions, k) && (v == "E" || v == "Z" || v == "")){
                curr_fa->double_bonds->double_bond_positions.insert({k, v});
                curr_fa->double_bonds->num_double_bonds = curr_fa->double_bonds->double_bond_positions.size();
            }
        }
    }
}


void FattyAcidParserEventHandler::reset_length(TreeNode *node) {
    tmp.set_int("length", 0);
    tmp.set_string("length_pattern", "");
    tmp.set_list("length_tokens", new GenericList());
    tmp.set_int("add_lengths", 1);
}


void FattyAcidParserEventHandler::set_functional_length(TreeNode *node) {
    if (tmp.get_int("length") != (int)tmp.get_list("fg_pos")->list.size()){
        throw LipidException("Length of functional group '" + std::to_string(tmp.get_int("length")) + "' does not match with number of its positions '" + std::to_string(tmp.get_list("fg_pos")->list.size()) + "'");
    }
}


void FattyAcidParserEventHandler::set_fatty_length(TreeNode *node) {
    tmp.set_int("add_lengths", 0);
}


void FattyAcidParserEventHandler::special_number(TreeNode *node) {
    if (tmp.get_int("add_lengths")){
        tmp.set_int("length", tmp.get_int("length") + special_numbers.at(node->get_text()));
        tmp.set_string("length_pattern", tmp.get_string("length_pattern") + "X");
        tmp.get_list("length_tokens")->add_int(special_numbers.at(node->get_text()));
    }
}


void FattyAcidParserEventHandler::last_number(TreeNode *node) {
    if (tmp.get_int("add_lengths")){
        tmp.set_int("length", tmp.get_int("length") + last_numbers.at(node->get_text()));
        tmp.set_string("length_pattern", tmp.get_string("length_pattern") + "L");
        tmp.get_list("length_tokens")->add_int(last_numbers.at(node->get_text()));
    }
}


void FattyAcidParserEventHandler::second_number(TreeNode *node) {
    if (tmp.get_int("add_lengths")){
        tmp.set_int("length", tmp.get_int("length") + second_numbers.at(node->get_text()));
        tmp.set_string("length_pattern", tmp.get_string("length_pattern") + "S");
        tmp.get_list("length_tokens")->add_int(second_numbers.at(node->get_text()));
    }
}


void FattyAcidParserEventHandler::set_functional_group(TreeNode *node) {
    tmp.set_list("fg_pos", new GenericList());
    tmp.set_string("fg_type", "");
}


void FattyAcidParserEventHandler::add_functional_group(TreeNode *node) {
    if (tmp.contains_key("added_func_group")){
        tmp.remove("added_func_group");
        return;
    }
    
    else if (tmp.contains_key("add_methylene")){ 
        tmp.remove("add_methylene");
        add_cyclo(node);
        return;
    }
    
    string t = tmp.get_string("fg_type");
    
    FunctionalGroup *fg = 0;
    if (t != "acetoxy"){
        if (uncontains_val(func_groups, t)){
            throw LipidException("Unknown functional group: '" + t + "'");
        }
        t = func_groups.at(t);
        if (t.length() == 0) return;
        fg = KnownFunctionalGroups::get_functional_group(t);
        if (fg == 0) {
            throw LipidException("Functional group not registered: '" + t + "'");
        }
    }
    else {
        fg = new AcylAlkylGroup(new FattyAcid("O", 2));
    }
    
    FattyAcid* fa = fatty_acyl_stack.back();
    if (uncontains_val_p(fa->functional_groups, t)) fa->functional_groups->insert({t, vector<FunctionalGroup*>()});
    int l = tmp.get_list("fg_pos")->list.size();
    for (int i = 0; i < l; ++i){
        int pos = tmp.get_list("fg_pos")->get_list(i)->get_int(0);
        
        int num_pos = 0;
        if (tmp.contains_key("reduction")){
            GenericList *gl = tmp.get_list("reduction");
            int l = gl->list.size();
            for (int i = 0; i < l; ++i){
                num_pos += gl->get_int(i) < pos;
            }
        }
        
        FunctionalGroup* fg_insert = fg->copy();
        fg_insert->position = pos - num_pos;
        fa->functional_groups->at(t).push_back(fg_insert);
    }
    delete fg;
}


void FattyAcidParserEventHandler::set_functional_pos(TreeNode *node) {
    GenericList* gl = tmp.get_list("fg_pos");
    int s = gl->list.size();
    gl->get_list(s - 1)->set_int(0, atoi(node->get_text().c_str()));
}


void FattyAcidParserEventHandler::set_functional_position(TreeNode *node) {
    GenericList* gl = new GenericList();
    gl->add_int(0);
    gl->add_string("");
    tmp.get_list("fg_pos")->add_list(gl);
}


void FattyAcidParserEventHandler::set_functional_type(TreeNode *node) {
    tmp.set_string("fg_type", node->get_text());
}


void FattyAcidParserEventHandler::rearrange_cycle(TreeNode *node) {
    if (tmp.contains_key("post_adding")){
        fatty_acyl_stack.back()->num_carbon += tmp.get_list("post_adding")->list.size();
        tmp.remove("post_adding");
    }
        
    FattyAcid* curr_fa = fatty_acyl_stack.back();
    int start = tmp.get_list("fg_pos")->get_list(0)->get_int(0);
    if (contains_val_p(curr_fa->functional_groups, "cy")){
        for (auto &cy : curr_fa->functional_groups->at("cy")){
            int shift_val = start - cy->position;
            if (shift_val == 0) continue;
            ((Cycle*)cy)->rearrange_functional_groups(curr_fa, shift_val);
        }
    }
}


void FattyAcidParserEventHandler::add_epoxy(TreeNode *node) {
    GenericList *gl = tmp.get_list("fg_pos");
    while(gl->list.size() > 1){
        gl->del(gl->list.back());
        gl->list.pop_back();
    }
    tmp.set_string("fg_type", "Epoxy");
}


void FattyAcidParserEventHandler::set_cycle(TreeNode *node) {
    tmp.set_int("cyclo", 1);
}


void FattyAcidParserEventHandler::set_methylene(TreeNode *node) {
    tmp.set_string("fg_type", "methylene");
    GenericList *gl = tmp.get_list("fg_pos");
    if (gl->list.size() > 1){
        if (gl->get_list(0)->get_int(0) < gl->get_list(1)->get_int(0)) {
            gl->get_list(1)->set_int(0, gl->get_list(1)->get_int(0) + 1);
        }
        else if (gl->get_list(0)->get_int(0) > gl->get_list(1)->get_int(0)){
            gl->get_list(0)->set_int(0, gl->get_list(0)->get_int(0) + 1);
        }
        fatty_acyl_stack.back()->num_carbon += 1;
        tmp.set_int("add_methylene", 1);
    }
}


void FattyAcidParserEventHandler::set_dioic(TreeNode *node) {
    headgroup = "FA";
    
    int pos = (tmp.get_list("fg_pos")->list.size() == 2) ? tmp.get_list("fg_pos")->get_list(1)->get_int(0) : fatty_acyl_stack.back()->num_carbon;
    fatty_acyl_stack.back()->num_carbon -= 1;
    if (tmp.contains_key("reduction")) {
        GenericList *gl = tmp.get_list("reduction");
        pos -= gl->list.size();
    }
    FunctionalGroup* func_group = KnownFunctionalGroups::get_functional_group("COOH");
    func_group->position = pos - 1;
    if (uncontains_val_p(fatty_acyl_stack.back()->functional_groups, "COOH")) fatty_acyl_stack.back()->functional_groups->insert({"COOH", vector<FunctionalGroup*>()});
    fatty_acyl_stack.back()->functional_groups->at("COOH").push_back(func_group);
}


void FattyAcidParserEventHandler::set_dial(TreeNode *node) {
    FattyAcid* curr_fa = fatty_acyl_stack.back();
    int pos = curr_fa->num_carbon;
    FunctionalGroup *fg = KnownFunctionalGroups::get_functional_group("oxo");
    fg->position = pos;
    if (uncontains_val_p(curr_fa->functional_groups, "oxo")) curr_fa->functional_groups->insert({"oxo", vector<FunctionalGroup*>()});
    curr_fa->functional_groups->at("oxo").push_back(fg);
}


void FattyAcidParserEventHandler::set_prosta(TreeNode *node) {
    int minus_pos = 0;
    if (tmp.contains_key("reduction")){
        GenericList *gl = tmp.get_list("reduction");
        for (int i = 0; i < (int)gl->list.size(); ++i){
            minus_pos += gl->get_int(i) < 8;
        }
    }
    tmp.set_list("fg_pos", new GenericList());
    tmp.get_list("fg_pos")->add_list(new GenericList());
    tmp.get_list("fg_pos")->add_list(new GenericList());
    tmp.get_list("fg_pos")->get_list(0)->add_int(8 - minus_pos);
    tmp.get_list("fg_pos")->get_list(0)->add_string("");
    tmp.get_list("fg_pos")->get_list(1)->add_int(12 - minus_pos);
    tmp.get_list("fg_pos")->get_list(1)->add_string("");
    tmp.set_string("fg_type", "cy");
}


void FattyAcidParserEventHandler::add_cyclo(TreeNode *node) {
    
    int start = tmp.get_list("fg_pos")->get_list(0)->get_int(0);
    int end = tmp.get_list("fg_pos")->get_list(1)->get_int(0);
    
    
    DoubleBonds *cyclo_db = new DoubleBonds();
    // check double bonds
    if (!fatty_acyl_stack.back()->double_bonds->double_bond_positions.empty()){
        for (auto &kv : fatty_acyl_stack.back()->double_bonds->double_bond_positions){
            if (start <= kv.first && kv.first <= end){
                cyclo_db->double_bond_positions.insert({kv.first, kv.second});
            }
        }
        cyclo_db->num_double_bonds = cyclo_db->double_bond_positions.size();
    
        for (auto &kv : cyclo_db->double_bond_positions){
            fatty_acyl_stack.back()->double_bonds->double_bond_positions.erase(kv.first);
        }
        fatty_acyl_stack.back()->double_bonds->num_double_bonds = fatty_acyl_stack.back()->double_bonds->double_bond_positions.size();
        
    }        
    // check functional_groups
    map<string, vector<FunctionalGroup*> > *cyclo_fg = new map<string, vector<FunctionalGroup*>>();
    set<string> remove_list;
    FattyAcid *curr_fa = fatty_acyl_stack.back();
    
    if (contains_val_p(curr_fa->functional_groups, "noyloxy")){
        vector<int> remove_item;
        int i = 0;
        for (auto &func_group : curr_fa->functional_groups->at("noyloxy")){
            if (start <= func_group->position && func_group->position <= end){
                CarbonChain *cc = new CarbonChain((FattyAcid*)func_group, func_group->position);
                
                if (uncontains_val_p(curr_fa->functional_groups, "cc")) curr_fa->functional_groups->insert({"cc", vector<FunctionalGroup*>()});
                curr_fa->functional_groups->at("cc").push_back(cc);
                remove_item.push_back(i);
            }
            ++i;
        }
        for (int i = remove_item.size() - 1; i >= 0; --i) curr_fa->functional_groups->at("noyloxy").erase(curr_fa->functional_groups->at("noyloxy").begin() + remove_item.at(i));
        if (curr_fa->functional_groups->at("noyloxy").empty()) remove_list.insert("noyloxy");
    }
    
    for (auto &kv : *(curr_fa->functional_groups)){
        vector<int> remove_item;
        int i = 0;
        for (auto &func_group : kv.second){
            if (start <= func_group->position && func_group->position <= end){
                if (uncontains_val_p(cyclo_fg, kv.first)) cyclo_fg->insert({kv.first, vector<FunctionalGroup*>()});
                cyclo_fg->at(kv.first).push_back(func_group);
                remove_item.push_back(i);
            }
            ++i;    
        }
        for (int i = remove_item.size() - 1; i >= 0; --i) kv.second.erase(kv.second.begin() + remove_item.at(i));
        if (kv.second.empty()) remove_list.insert(kv.first);
    }
    for (auto &fg : remove_list) curr_fa->functional_groups->erase(fg);
    vector<Element>* bridge_chain = new vector<Element>{};
    if (tmp.contains_key("furan")){
        tmp.remove("furan");
        bridge_chain->push_back(ELEMENT_O);
    }
    
    Cycle *cycle = new Cycle(end - start + 1 + bridge_chain->size(), start, end, cyclo_db, cyclo_fg, bridge_chain);
    if (uncontains_val_p(fatty_acyl_stack.back()->functional_groups, "cy")) fatty_acyl_stack.back()->functional_groups->insert({"cy", vector<FunctionalGroup*>()});
    fatty_acyl_stack.back()->functional_groups->at("cy").push_back(cycle);
}


void FattyAcidParserEventHandler::reduction(TreeNode *node) {
    int reduction_num = tmp.get_list("fg_pos")->list.size();
    fatty_acyl_stack.back()->num_carbon -= reduction_num;
    for (auto &kv : *(fatty_acyl_stack.back()->functional_groups)){
        for (auto &func_group : kv.second){
            func_group->shift_positions(-reduction_num);
        }
    }
    tmp.set_list("reduction", new GenericList());
    for (int i = 0; i < (int)tmp.get_list("fg_pos")->list.size(); ++i){
        tmp.get_list("reduction")->add_int(tmp.get_list("fg_pos")->get_list(i)->get_int(0));
    }
}


void FattyAcidParserEventHandler::homo(TreeNode *node) {
    tmp.set_list("post_adding", new GenericList());
    for (int i = 0; i < (int)tmp.get_list("fg_pos")->list.size(); ++i){
        tmp.get_list("post_adding")->add_int(tmp.get_list("fg_pos")->get_list(i)->get_int(0));
    }
}


void FattyAcidParserEventHandler::set_recursion(TreeNode *node) {
    tmp.set_list("fg_pos", new GenericList());
    tmp.set_string("fg_type", "");
    fatty_acyl_stack.push_back(new FattyAcid("FA"));
    tmp.set_dictionary(FA_I, new GenericDictionary());
    tmp.get_dictionary(FA_I)->set_int("recursion_pos", 0);
}



void FattyAcidParserEventHandler::set_tetrahydrofuran(TreeNode *node){
    tmp.set_int("furan", 1);
    tmp.set_int("tetrahydrofuran", 1);
    set_cycle(node);
}



void FattyAcidParserEventHandler::set_furan(TreeNode *node){
    tmp.set_int("furan", 1);
    set_cycle(node);
}



void FattyAcidParserEventHandler::add_recursion(TreeNode *node) {
        int pos = tmp.get_dictionary(FA_I)->get_int("recursion_pos");
        
        FattyAcid *fa = fatty_acyl_stack.back();
        fatty_acyl_stack.pop_back();
        fa->position = pos;
        FattyAcid *curr_fa = fatty_acyl_stack.back();
        
        string fname = "";
        if (tmp.contains_key("cyclo_yl")){
            fname = "cyclo";
            tmp.remove("cyclo_yl");
        }
        else {
            fname = headgroup;
        }
        if (uncontains_val_p(curr_fa->functional_groups, fname)) curr_fa->functional_groups->insert({fname, vector<FunctionalGroup*>()});
        curr_fa->functional_groups->at(fname).push_back(fa);
        tmp.set_int("added_func_group", 1);
}


void FattyAcidParserEventHandler::set_recursion_pos(TreeNode *node) {
    tmp.get_dictionary(FA_I)->set_int("recursion_pos", atoi(node->get_text().c_str()));
}


void FattyAcidParserEventHandler::set_yl_ending(TreeNode *node) {
    int l = atoi(node->get_text().c_str()) - 1;
    if (l == 0) return;

    FattyAcid *curr_fa = fatty_acyl_stack.back();
    string fname = "";
    
    if (tmp.contains_key("furan")){
        curr_fa->num_carbon -= l;
        return;
    }
    
    
    FunctionalGroup *fg = 0;
    if (l == 1){
        fname = "Me";
        fg = KnownFunctionalGroups::get_functional_group(fname);
    }
    else if (l == 2){
        fname = "Et";
        fg = KnownFunctionalGroups::get_functional_group(fname);
    }
    else {
        FattyAcid *fa = new FattyAcid("FA", l);
        // shift functional groups
        for (auto &kv : *(curr_fa->functional_groups)){
            vector<int> remove_item;
            int i = 0;
            for (auto &func_group : kv.second){
                if (func_group->position <= l){
                    remove_item.push_back(i);
                    if (uncontains_val_p(fa->functional_groups, kv.first)) fa->functional_groups->insert({kv.first, vector<FunctionalGroup*>()});
                    func_group->position = l + 1 - func_group->position;
                    fa->functional_groups->at(kv.first).push_back(func_group);
                }
            }
            for (int i = remove_item.size() - 1; i >= 0; --i) curr_fa->functional_groups->at(kv.first).erase(curr_fa->functional_groups->at(kv.first).begin() + remove_item.at(i));
        }
        map<string, vector<FunctionalGroup*> > *tmp = curr_fa->functional_groups;
        curr_fa->functional_groups = new map<string, vector<FunctionalGroup*> >();
        for (auto &kv : *tmp){
            if (!kv.second.empty()) curr_fa->functional_groups->insert({kv.first, kv.second});
        }
        delete tmp;
        
        // shift double bonds
        if (!curr_fa->double_bonds->double_bond_positions.empty()){
            delete fa->double_bonds;
            fa->double_bonds = new DoubleBonds();
            for (auto &kv : curr_fa->double_bonds->double_bond_positions){
                if (kv.first <= l) fa->double_bonds->double_bond_positions.insert({l + 1 - kv.first, kv.second});
            }
            fa->double_bonds->num_double_bonds = fa->double_bonds->double_bond_positions.size();
            for (auto &kv : fa->double_bonds->double_bond_positions) curr_fa->double_bonds->double_bond_positions.erase(kv.first);
        }
        fname = "cc";
        fg = new CarbonChain(fa);
    }
    curr_fa->num_carbon -= l;
    fg->position = l;
    curr_fa->shift_positions(-l);
    if (uncontains_val_p(curr_fa->functional_groups, fname)) curr_fa->functional_groups->insert({fname, vector<FunctionalGroup*>()});
    curr_fa->functional_groups->at(fname).push_back(fg);
}


void FattyAcidParserEventHandler::set_acetic_acid(TreeNode *node) {
    fatty_acyl_stack.back()->num_carbon += 2;
    headgroup = "FA";
}


void FattyAcidParserEventHandler::add_hydroxyl(TreeNode *node) {
    int h = atoi(node->get_text().c_str());
    tmp.get_list("hydroxyl_pos")->add_int(h);
}


void FattyAcidParserEventHandler::setup_hydroxyl(TreeNode *node) {
    tmp.set_list("hydroxyl_pos", new GenericList());
}


void FattyAcidParserEventHandler::add_hydroxyls(TreeNode *node) {       
    
    if (tmp.get_list("hydroxyl_pos")->list.size() > 1){
        FunctionalGroup *fg_oh = KnownFunctionalGroups::get_functional_group("OH");
        vector<int> sorted_pos;
        for (int i = 0; i < (int)tmp.get_list("hydroxyl_pos")->list.size(); ++i) sorted_pos.push_back(tmp.get_list("hydroxyl_pos")->get_int(i));
        std::sort(sorted_pos.rbegin(), sorted_pos.rend());
        for (int i = 0; i < (int)sorted_pos.size() - 1; ++i){
            int pos = sorted_pos.at(i);
            FunctionalGroup *fg_insert = fg_oh->copy();
            fg_insert->position = pos;
            if (uncontains_val_p(fatty_acyl_stack.back()->functional_groups, "OH")) fatty_acyl_stack.back()->functional_groups->insert({"OH", vector<FunctionalGroup*>()});
            fatty_acyl_stack.back()->functional_groups->at("OH").push_back(fg_insert);
        }
        delete fg_oh;
    }
}


void FattyAcidParserEventHandler::add_wax_ester(TreeNode *node) {
    FattyAcid *fa = fatty_acyl_stack.back();
    fatty_acyl_stack.pop_back();
    
    fa->lipid_FA_bond_type = ETHER;
    fatty_acyl_stack.insert(fatty_acyl_stack.begin(), fa);
}


void FattyAcidParserEventHandler::set_ate(TreeNode *node) {
    fatty_acyl_stack.back()->num_carbon += ate.at(node->get_text());
    headgroup = "WE";
}


void FattyAcidParserEventHandler::set_iso(TreeNode *node) {
        FattyAcid *curr_fa = fatty_acyl_stack.back();
        curr_fa->num_carbon -= 1;
        FunctionalGroup *fg = KnownFunctionalGroups::get_functional_group("Me");
        fg->position = 2;
        if (uncontains_val_p(curr_fa->functional_groups, "Me")) curr_fa->functional_groups->insert({"Me", vector<FunctionalGroup*>()});
        curr_fa->functional_groups->at("Me").push_back(fg);
}


void FattyAcidParserEventHandler::set_coa(TreeNode *node) {
    headgroup = "CoA";
}


void FattyAcidParserEventHandler::set_methyl(TreeNode *node) {
    fatty_acyl_stack.back()->num_carbon += 1;
}


void FattyAcidParserEventHandler::set_car(TreeNode *node) {
    tmp.set_list("fg_pos", new GenericList());
    tmp.set_string("fg_type", "");
}


void FattyAcidParserEventHandler::add_car(TreeNode *node) {
    headgroup = "CAR";
}


void FattyAcidParserEventHandler::add_ethanolamine(TreeNode *node) {
    headgroup = "NAE";
}


void FattyAcidParserEventHandler::add_amine(TreeNode *node) {
    FattyAcid *fa = fatty_acyl_stack.back();
    fatty_acyl_stack.pop_back();
    
    fa->lipid_FA_bond_type = AMIDE;
    fatty_acyl_stack[fatty_acyl_stack.size() - 1]->lipid_FA_bond_type = AMIDE;
    fatty_acyl_stack.insert(fatty_acyl_stack.begin(), fa);
}


void FattyAcidParserEventHandler::add_amine_name(TreeNode *node) {
    headgroup = "NA";
}


void FattyAcidParserEventHandler::add_summary(TreeNode *node) {    
    tmp.get_dictionary(FA_I)->set_dictionary("fg_pos_summary", new GenericDictionary());
    for (int i = 0; i < (int)tmp.get_list("fg_pos")->list.size(); ++i){
        string k = std::to_string(tmp.get_list("fg_pos")->get_list(i)->get_int(0));
        string v = to_upper(tmp.get_list("fg_pos")->get_list(i)->get_string(1));
        tmp.get_dictionary(FA_I)->get_dictionary("fg_pos_summary")->set_string(k, v);
    }
}


void FattyAcidParserEventHandler::add_func_stereo(TreeNode *node) {
    int l = tmp.get_list("fg_pos")->list.size();
    tmp.get_list("fg_pos")->get_list(l - 1)->set_string(1, node->get_text());
}


void FattyAcidParserEventHandler::open_db_length(TreeNode *node) {
    tmp.set_int("add_lengths", 1);
}


void FattyAcidParserEventHandler::close_db_length(TreeNode *node) {
    tmp.set_int("add_lengths", 0);
}
