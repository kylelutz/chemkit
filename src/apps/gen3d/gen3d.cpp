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

#include <chemkit/atom.h>
#include <chemkit/point3.h>
#include <chemkit/chemkit.h>
#include <chemkit/foreach.h>
#include <chemkit/forcefield.h>
#include <chemkit/lineformat.h>
#include <chemkit/moleculefile.h>
#include <chemkit/coordinatepredictor.h>
#include <chemkit/moleculegeometryoptimizer.h>

void printHelp(char *argv[], const boost::program_options::options_description &options)
{
    std::cout << "Usage: " << argv[0] << " [OPTIONS] formula file\n";
    std::cout << "\n";
    std::cout << "Generates 3D coordinates for a molecule from its SMILES formula.\n";
    std::cout << "\n";
    std::cout << "Options:\n";
    std::cout << options << "\n";
}

int main(int argc, char *argv[])
{
    std::string inputFormula;
    std::string inputFormatName;
    std::string outputFileName;
    std::string outputFormatName;

    boost::program_options::options_description options;
    options.add_options()
        ("formula",
            boost::program_options::value<std::string>(&inputFormula),
            "The input formula.")
        ("output-file",
            boost::program_options::value<std::string>(&outputFileName),
            "The output file.")
        ("input-format,i",
            boost::program_options::value<std::string>(&inputFormatName),
            "Sets the input format.")
        ("output-format,o",
            boost::program_options::value<std::string>(&outputFormatName),
            "Sets the output format.")
        ("no-optimization",
            "Do not perform geometry optimization.")
        ("help,h",
            "Shows this help message");

    boost::program_options::positional_options_description positionalOptions;
    positionalOptions.add("formula", 1).add("output-file", 1);

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
    else if(inputFormula.empty()){
        printHelp(argv, options);
        std::cerr << "Error: No input formula specified." << std::endl;
        return -1;
    }
    else if(outputFileName.empty()){
        printHelp(argv, options);
        std::cerr << "Error: No output file specified." << std::endl;
        return -1;
    }

    // check input format
    if(inputFormatName.empty()){
        // default to smiles if not format specified
        inputFormatName = "smiles";
    }

    // create input line format
    boost::scoped_ptr<chemkit::LineFormat> inputFormat(chemkit::LineFormat::create(inputFormatName));
    if(!inputFormat){
        if(inputFormatName.empty()){
            std::cerr << "Failed to guess input format." << std::endl;
        }
        else{
            std::cerr << "Input format: " << inputFormatName << " is not supported." << std::endl;
        }

        return -1;
    }

    // read input formula
    boost::shared_ptr<chemkit::Molecule> molecule(inputFormat->read(inputFormula));
    if(!molecule){
        std::cerr << "Failed to parse formula: " << inputFormat->errorString() << std::endl;
        return -1;
    }

    // generate 3d coordinates
    chemkit::CoordinatePredictor::predictCoordinates(molecule.get());

    if(!variables.count("no-optimization")){
        // optimize 3d coordinates
        chemkit::MoleculeGeometryOptimizer::optimizeCoordinates(molecule.get());
    }

    // set center to origin
    molecule->setCenter(0, 0, 0);

    // write output file
    chemkit::MoleculeFile outputFile(outputFileName);
    outputFile.addMolecule(molecule);
    if(!outputFormatName.empty()){
        if(!outputFile.setFormat(outputFormatName)){
            std::cerr << "File format '" << outputFormatName << "' is not supported." << std::endl;
            return -1;
        }
    }

    if(!outputFile.write()){
        std::cerr << "Error: failed to write output file: " << outputFile.errorString() << std::endl;
        return -1;
    }

    return 0;
}
