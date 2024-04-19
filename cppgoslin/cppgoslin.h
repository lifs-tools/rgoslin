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


#ifndef CPPGOSLIN_H
#define CPPGOSLIN_H

#include "cppgoslin/domain/Adduct.h"
#include "cppgoslin/domain/Cycle.h"
#include "cppgoslin/domain/DoubleBonds.h"
#include "cppgoslin/domain/Element.h"
#include "cppgoslin/domain/FattyAcid.h"
#include "cppgoslin/domain/FunctionalGroup.h"
#include "cppgoslin/domain/GenericDatastructures.h"
#include "cppgoslin/domain/Headgroup.h"
#include "cppgoslin/domain/LipidAdduct.h"
#include "cppgoslin/domain/LipidEnums.h"
#include "cppgoslin/domain/LipidExceptions.h"
#include "cppgoslin/domain/LipidCompleteStructure.h"
#include "cppgoslin/domain/LipidFullStructure.h"
#include "cppgoslin/domain/LipidStructureDefined.h"
#include "cppgoslin/domain/LipidSnPosition.h"
#include "cppgoslin/domain/LipidMolecularSpecies.h"
#include "cppgoslin/domain/LipidSpecies.h"
#include "cppgoslin/domain/LipidSpeciesInfo.h"
#include "cppgoslin/domain/StringFunctions.h"
#include "cppgoslin/parser/LipidBaseParserEventHandler.h"
#include "cppgoslin/parser/BaseParserEventHandler.h"
#include "cppgoslin/parser/GoslinParserEventHandler.h"
#include "cppgoslin/parser/ShorthandParserEventHandler.h"
#include "cppgoslin/parser/SwissLipidsParserEventHandler.h"
#include "cppgoslin/parser/SumFormulaParserEventHandler.h"
#include "cppgoslin/parser/HmdbParserEventHandler.h"
#include "cppgoslin/parser/KnownGrammars.h"
#include "cppgoslin/parser/KnownParsers.h"
#include "cppgoslin/parser/LipidMapsParserEventHandler.h"
#include "cppgoslin/parser/FattyAcidParserEventHandler.h"
#include "cppgoslin/parser/Parser.h"
#include "cppgoslin/parser/SumFormulaParser.h"

#endif /* CPPGOSLIN_H */
