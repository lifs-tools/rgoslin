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


#include "cppgoslin/parser/KnownParsers.h"

FattyAcidParser::FattyAcidParser() : Parser<LipidAdduct*>(new FattyAcidParserEventHandler(), GrammarString(fatty_acid_grammar), DEFAULT_QUOTE){
    Headgroup::init();
}


FattyAcidParser::~FattyAcidParser(){
    delete parser_event_handler;
}

LipidAdduct* FattyAcidParser::parse_parallel(string lipid_name, bool b, BaseParserEventHandler<LipidAdduct*>* bpeh){
    FattyAcidParserEventHandler* h = new FattyAcidParserEventHandler();
    LipidAdduct* l = Parser<LipidAdduct*>::parse_parallel(to_lower(lipid_name), b, h);
    delete h;
    return l;
}

LipidAdduct* FattyAcidParser::parse(string lipid_name, bool throw_error){
    return Parser<LipidAdduct*>::parse(to_lower(lipid_name), throw_error);
}







ShorthandParser::ShorthandParser() : Parser<LipidAdduct*>(new ShorthandParserEventHandler(), GrammarString(shorthand_grammar), DEFAULT_QUOTE){
    Headgroup::init();
}

ShorthandParser::~ShorthandParser(){
    delete parser_event_handler;
}

LipidAdduct* ShorthandParser::parse_parallel(string lipid_name, bool b, BaseParserEventHandler<LipidAdduct*>* bpeh){
    ShorthandParserEventHandler* h = new ShorthandParserEventHandler();
    LipidAdduct* l = Parser<LipidAdduct*>::parse_parallel(lipid_name, b, h);
    delete h;
    return l;
}





GoslinParser::GoslinParser() : Parser<LipidAdduct*>(new GoslinParserEventHandler(), GrammarString(goslin_grammar), DEFAULT_QUOTE){
    Headgroup::init();
}

GoslinParser::~GoslinParser(){
    delete parser_event_handler;
}

LipidAdduct* GoslinParser::parse_parallel(string lipid_name, bool b, BaseParserEventHandler<LipidAdduct*>* bpeh){
    GoslinParserEventHandler* h = new GoslinParserEventHandler();
    LipidAdduct* l = Parser<LipidAdduct*>::parse_parallel(lipid_name, b, h);
    delete h;
    return l;
}






LipidMapsParser::LipidMapsParser() : Parser<LipidAdduct*>(new LipidMapsParserEventHandler(), GrammarString(lipid_maps_grammar), DEFAULT_QUOTE){
    Headgroup::init();
}

LipidMapsParser::~LipidMapsParser(){
    delete parser_event_handler;
}

LipidAdduct* LipidMapsParser::parse_parallel(string lipid_name, bool b, BaseParserEventHandler<LipidAdduct*>* bpeh){
    LipidMapsParserEventHandler* h = new LipidMapsParserEventHandler();
    LipidAdduct* l = Parser<LipidAdduct*>::parse_parallel(lipid_name, b, h);
    delete h;
    return l;
}






SwissLipidsParser::SwissLipidsParser() : Parser<LipidAdduct*>(new SwissLipidsParserEventHandler(), GrammarString(swiss_lipids_grammar), DEFAULT_QUOTE){
    Headgroup::init();
}

SwissLipidsParser::~SwissLipidsParser(){
    delete parser_event_handler;
}

LipidAdduct* SwissLipidsParser::parse_parallel(string lipid_name, bool b, BaseParserEventHandler<LipidAdduct*>* bpeh){
    SwissLipidsParserEventHandler* h = new SwissLipidsParserEventHandler();
    LipidAdduct* l = Parser<LipidAdduct*>::parse_parallel(lipid_name, b, h);
    delete h;
    return l;
}






HmdbParser::HmdbParser() : Parser<LipidAdduct*>(new HmdbParserEventHandler(), GrammarString(hmdb_grammar), DEFAULT_QUOTE){
    Headgroup::init();
}

HmdbParser::~HmdbParser(){
    delete parser_event_handler;
}

LipidAdduct* HmdbParser::parse_parallel(string lipid_name, bool b, BaseParserEventHandler<LipidAdduct*>* bpeh){
    HmdbParserEventHandler* h = new HmdbParserEventHandler();
    LipidAdduct* l = Parser<LipidAdduct*>::parse_parallel(lipid_name, b, h);
    delete h;
    return l;
}






LipidParser::LipidParser(){
    parser_list.push_back(new ShorthandParser());
    parser_list.push_back(new GoslinParser());
    parser_list.push_back(new FattyAcidParser());
    parser_list.push_back(new LipidMapsParser());
    parser_list.push_back(new SwissLipidsParser());
    parser_list.push_back(new HmdbParser());
    lastSuccessfulParser = 0;
}

LipidParser::~LipidParser(){
    for (auto parser : parser_list) delete parser;
}
    
    
LipidAdduct* LipidParser::parse(string lipid_name){
    lastSuccessfulParser = 0;
    
    for (auto parser : parser_list) {
        LipidAdduct *lipid = parser->parse(lipid_name, false);
        if (lipid){
            lastSuccessfulParser = parser;
            return lipid;
        }
    }
    throw LipidException("Lipid not found");
}
    
    
LipidAdduct* LipidParser::parse_parallel(string lipid_name){
    lastSuccessfulParser = 0;
    
    for (auto parser : parser_list) {
        LipidAdduct *lipid = parser->parse_parallel(lipid_name, false);
        if (lipid){
            lastSuccessfulParser = parser;
            return lipid;
        }
    }
    throw LipidException("Lipid not found");
}
