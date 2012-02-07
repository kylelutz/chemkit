/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
** All rights reserved.
**
** This file is a part of the chemkit project. For more information
** see <http://www.chemkit.org>.
**
** Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions
** are met:
**
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in the
**     documentation and/or other materials provided with the distribution.
**   * Neither the name of the chemkit project nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
******************************************************************************/

#include <string>
#include <iostream>

#include <boost/scoped_ptr.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <chemkit/chemkit.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>
#include <chemkit/lineformat.h>
#include <chemkit/moleculefile.h>
#include <chemkit/substructurequery.h>

void printHelp(char *argv[], const boost::program_options::options_description &options)
{
    std::cout << "Usage: " << argv[0] << " [OPTIONS] PATTERN FILE\n";
    std::cout << "\n";
    std::cout << "Search for molecules matching PATTERN in FILE. PATTERN is a line\n";
    std::cout << "representation (e.g. InChI or SMILES) of a molecule to search\n";
    std::cout << "for. A matching molecule is either an exact match or a\n";
    std::cout << "superstructure of PATTERN.\n";
    std::cout << "\n";
    std::cout << "Options:\n";
    std::cout << options << "\n";
}

int main(int argc, char *argv[])
{
    std::string formula;
    std::string fileName;

    boost::program_options::options_description options;
    options.add_options()
        ("formula",
            boost::program_options::value<std::string>(&formula),
            "Input formula to match against.")
        ("file",
            boost::program_options::value<std::string>(&fileName),
            "Input file to search.")
        ("composition,c",
            "Match composition rather than structure.")
        ("exact-match,e",
            "Return molecules that exactly match PATTERN.")
        ("invert-match,v",
            "Return only non-matching molecules.")
        ("names-only,n",
            "Output only the names of matching molecules.")
        ("help,h",
            "Shows this help message");

    boost::program_options::positional_options_description positionalOptions;
    positionalOptions.add("formula", 1).add("file", 1);

    boost::program_options::variables_map variables;
    boost::program_options::store(
        boost::program_options::command_line_parser(argc, argv)
            .options(options)
            .positional(positionalOptions).run(),
        variables);
    boost::program_options::notify(variables);

    if(variables.count("help")){
        printHelp(argv, options);
        return 0;
    }
    else if(formula.empty()){
        printHelp(argv, options);
        std::cerr << "Error: no input formula given." << std::endl;
        return -1;
    }
    else if(fileName.empty()){
        printHelp(argv, options);
        std::cerr << "Error: no input file given." << std::endl;
        return -1;
    }

    // select input line format based on the input pattern given
    std::string inputFormat;
    if(boost::algorithm::starts_with(formula, "InChI=") || isdigit(formula[0])){
        inputFormat = "inchi";
    }
    else{
        inputFormat = "smiles";
    }

    // create input line format
    boost::scoped_ptr<chemkit::LineFormat> patternFormat(chemkit::LineFormat::create(inputFormat));
    if(!patternFormat){
        std::cerr << "Error: failed to create line format." << std::endl;
        return -1;
    }

    // read input molecule
    boost::scoped_ptr<chemkit::Molecule> patternMolecule(patternFormat->read(formula));
    if(!patternMolecule){
        std::cerr << "Error: failed to read pattern molecule: " << patternFormat->errorString() << std::endl;
        return -1;
    }

    // read input file
    chemkit::MoleculeFile inputFile(fileName);
    if(!inputFile.read()){
        std::cerr << "Error: failed to read input file: " << inputFile.errorString() << std::endl;
        return -1;
    }

    // get options
    bool compositionOnly = variables.find("composition") != variables.end();
    bool exactMatch = variables.find("exact-match") != variables.end();
    bool invertMatch = variables.find("invert-match") != variables.end();
    bool namesOnly = variables.find("names-only") != variables.end();

    int flags = 0;
    if(compositionOnly){
        flags |= chemkit::SubstructureQuery::CompareAtomsOnly;
    }
    if(exactMatch){
        flags |= chemkit::SubstructureQuery::CompareExact;
    }

    chemkit::SubstructureQuery query;
    query.setMolecule(patternMolecule.get());
    query.setFlags(flags);

    chemkit::MoleculeFile outputFile;

    std::vector<chemkit::Molecule *> matchingMolecules;

    foreach(const boost::shared_ptr<chemkit::Molecule> &molecule, inputFile.molecules()){
        bool match = query.matches(molecule.get());

        if((match && !invertMatch) || (!match && invertMatch)){
            outputFile.addMolecule(molecule);
        }
    }

    if(namesOnly){
        foreach(const boost::shared_ptr<chemkit::Molecule> &molecule, outputFile.molecules()){
            std::cout << molecule->name() << "\n";
        }
    }
    else{
        bool ok = outputFile.write(std::cout, inputFile.formatName());
        if(!ok){
            std::cerr << "Error: failed to write output file: " << outputFile.errorString() << std::endl;
            return -1;
        }
    }

    return 0;
}
