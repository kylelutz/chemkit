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

#include <iostream>

#include <boost/scoped_ptr.hpp>
#include <boost/program_options.hpp>

#include <chemkit/molecule.h>
#include <chemkit/lineformat.h>

int main(int argc, char *argv[])
{
    std::string inputFormula;
    std::string inputFormatName;
    std::string outputFormatName;

    boost::program_options::options_description options;
    options.add_options()
        ("formula",
            boost::program_options::value<std::string>(&inputFormula))
        ("input,i",
            boost::program_options::value<std::string>(&inputFormatName),
            "Sets the input format.")
        ("output,o",
            boost::program_options::value<std::string>(&outputFormatName),
            "Sets the output format.")
        ("help",
            "Shows this help message.");

    boost::program_options::positional_options_description positionalOptions;
    positionalOptions.add("formula", -1);

    boost::program_options::variables_map variables;
    boost::program_options::store(
        boost::program_options::command_line_parser(argc, argv)
            .options(options)
            .positional(positionalOptions).run(),
        variables);
    boost::program_options::notify(variables);

    if(variables.count("help") || inputFormula.empty()){
        // output help message and exit
        std::cout << "Usage: " << argv[0] << " [OPTION]... FORMULA\n";
        std::cout << "\n";
        std::cout << "Translates a chemical formula to a different format.\n";
        std::cout << "\n";
        std::cout << "Options:\n";
        std::cout << options << "\n";
        return 0;
    }

    if(inputFormatName.empty()){
        std::cerr << "Error: No input format specified." << std::endl;
        return -1;
    }
    else if(outputFormatName.empty()){
        std::cerr << "Error: No output format specified." << std::endl;
        return -1;
    }

    // create input format
    boost::scoped_ptr<chemkit::LineFormat> inputFormat(chemkit::LineFormat::create(inputFormatName));
    if(!inputFormat){
        std::cerr << "Error: Input format '" << inputFormatName << "' is not supported." << std::endl;
        return -1;
    }

    // read input formula
    boost::scoped_ptr<chemkit::Molecule> molecule(inputFormat->read(inputFormula));
    if(!molecule){
        std::cerr << "Error: Failed to read formula: " << inputFormat->errorString() << std::endl;
        return -1;
    }

    // create output format
    boost::scoped_ptr<chemkit::LineFormat> outputFormat(chemkit::LineFormat::create(outputFormatName));
    if(!outputFormat){
        std::cerr << "Error: Output format '" << outputFormatName << "' is not supported." << std::endl;
        return -1;
    }

    // write output formula
    std::string outputFormula = outputFormat->write(molecule.get());
    if(outputFormula.empty()){
        std::cerr << "Error: Failed to write output: " << outputFormat->errorString() << std::endl;
        return -1;
    }

    // print output formula
    std::cout << outputFormula << std::endl;

    return 0;
}
