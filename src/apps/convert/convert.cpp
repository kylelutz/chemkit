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
#include <chemkit/moleculefile.h>

void printHelp(char *argv[], const boost::program_options::options_description &options)
{
    std::cout << "Usage: " << argv[0] << " [OPTIONS] inputFile outputFile\n";
    std::cout << "\n";
    std::cout << "Converts a chemical input file to a new file with\n";
    std::cout << "a different file format.\n";
    std::cout << "\n";
    std::cout << "Options:\n";
    std::cout << options << "\n";
}

int main(int argc, char *argv[])
{
    std::string inputFileName;
    std::string inputFormatName;
    std::string outputFileName;
    std::string outputFormatName;

    boost::program_options::options_description options;
    options.add_options()
        ("input-file",
            boost::program_options::value<std::string>(&inputFileName),
            "The input file.")
        ("output-file",
            boost::program_options::value<std::string>(&outputFileName),
            "The output file.")
        ("input-format,i",
            boost::program_options::value<std::string>(&inputFormatName),
            "Sets the input format.")
        ("output-format,o",
            boost::program_options::value<std::string>(&outputFormatName),
            "Sets the output format.")
        ("help,h",
            "Shows this help message");

    boost::program_options::positional_options_description positionalOptions;
    positionalOptions.add("input-file", 1).add("output-file", 1);

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
    else if(inputFileName.empty()){
        printHelp(argv, options);
        std::cerr << "Error: No input file specified." << std::endl;
        return -1;
    }
    else if(outputFileName.empty()){
        printHelp(argv, options);
        std::cerr << "Error: No output file specified." << std::endl;
        return -1;
    }

    // read input
    chemkit::MoleculeFile inputFile;
    if(!inputFormatName.empty()){
        inputFile.setFormat(inputFormatName);
    }

    bool ok = false;
    if(inputFileName == "-"){
        ok = inputFile.read(std::cin);
    }
    else{
        ok = inputFile.read(inputFileName);
    }

    if(!ok){
        std::cerr << "Error: Failed to read input file: " << inputFile.errorString() << std::endl;
        return -1;
    }

    // write output
    if(outputFormatName.empty()){
        ok = inputFile.write(outputFileName);
    }
    else{
        if(outputFileName == "-"){
            ok = inputFile.write(std::cout, outputFormatName);
        }
        else{
            ok = inputFile.write(outputFileName, outputFormatName);
        }
    }

    if(!ok){
        std::cerr << "Error: failed to write output file: " << inputFile.errorString() << std::endl;
        return -1;
    }

    return 0;
}
