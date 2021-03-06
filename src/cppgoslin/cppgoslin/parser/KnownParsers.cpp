/*
MIT License

Copyright (c) 2020 Dominik Kopczynski   -   dominik.kopczynski {at} isas.de
                   Nils Hoffmann  -  nils.hoffmann {at} isas.de

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

GoslinParser::GoslinParser() : Parser<LipidAdduct*>(new GoslinParserEventHandler(), GrammarString(goslin_grammar), DEFAULT_QUOTE){
        
}


GoslinParser::~GoslinParser(){
    delete parser_event_handler;
}


GoslinFragmentParser::GoslinFragmentParser() : Parser<LipidAdduct*>(new GoslinFragmentParserEventHandler(), GrammarString(goslin_fragment_grammar), DEFAULT_QUOTE){
        
}


GoslinFragmentParser::~GoslinFragmentParser(){
    delete parser_event_handler;
}


LipidMapsParser::LipidMapsParser() : Parser<LipidAdduct*>(new LipidMapsParserEventHandler(), GrammarString(lipid_maps_grammar), DEFAULT_QUOTE){
        
}

LipidMapsParser::~LipidMapsParser(){
    delete parser_event_handler;
}


SwissLipidsParser::SwissLipidsParser() : Parser<LipidAdduct*>(new SwissLipidsParserEventHandler(), GrammarString(swiss_lipids_grammar), DEFAULT_QUOTE){
        
}

SwissLipidsParser::~SwissLipidsParser(){
    delete parser_event_handler;
}


HmdbParser::HmdbParser() : Parser<LipidAdduct*>(new HmdbParserEventHandler(), GrammarString(hmdb_grammar), DEFAULT_QUOTE){
        
}

HmdbParser::~HmdbParser(){
    delete parser_event_handler;
}



LipidParser::LipidParser(){
    parser_list.push_back(new GoslinParser());
    parser_list.push_back(new LipidMapsParser());
    parser_list.push_back(new SwissLipidsParser());
    parser_list.push_back(new HmdbParser());
    parser_list.push_back(new GoslinFragmentParser());
}

LipidParser::~LipidParser(){
    for (auto parser : parser_list) delete parser;
}
    
    
LipidAdduct* LipidParser::parse(string lipid_name){
    
    for (auto parser : parser_list) {
        LipidAdduct *lipid = parser->parse(lipid_name, false);
        if (lipid){
            return lipid;
        }
    }
    throw LipidException("Lipid not found");
}
