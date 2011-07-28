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

#include <QtCore>

#include <getopt.h>
#include <iostream>

#include <chemkit/chemkit.h>
#include <chemkit/molecule.h>
#include <chemkit/lineformat.h>
#include <chemkit/moleculefile.h>

// option flags
int COMPOSITION_FLAG = 0;
int EXACT_MATCH_FLAG = 0;
int NAMES_ONLY_FLAG = 0;
int INVERT_MATCH_FLAG = 0;

struct option options[] = {
    {"composition", no_argument, &COMPOSITION_FLAG, 'c'},
    {"exact-match", no_argument, &EXACT_MATCH_FLAG, 'e'},
    {"help", no_argument, 0, 'h'},
    {"names-only", no_argument, &NAMES_ONLY_FLAG, 'n'},
    {"invert-match", no_argument, &INVERT_MATCH_FLAG, 'v'},
    {0, 0, 0, 0},
};

void printHelp()
{
    QTextStream out(stdout);

    out << QString("Usage: %1 [OPTIONS] PATTERN FILE\n").arg(qAppName());
    out << "\n";
    out << "Search for molecules matching PATTERN in FILE. PATTERN is a line\n";
    out << "representation (e.g. InChI or SMILES) of a molecule to search\n";
    out << "for. A matching molecule is either an exact match or a\n";
    out << "superstructure of PATTERN.\n";
    out << "\n";
    out << "Options:\n";
    out << "  -c, --composition    Match composition rather than structure.\n";
    out << "  -e, --exact-match    Return molecules that exactly match PATTERN.\n";
    out << "  -h, --help           Display this help and exit.\n";
    out << "  -n, --names-only     Output only the names of matching molecules.\n";
    out << "  -v, --invert-match   Return only non-matching molecules.\n";
    out << "\n";
}

int main(int argc, char *argv[])
{
    QCoreApplication app(argc, argv);
    app.setApplicationName(argv[0]);

    QTextStream out(stdout);
    QTextStream err(stderr);

    QString inputPattern;
    QString inputFileName;

    // parse options
    for(;;){
        int optionIndex = 0;

        int c = getopt_long_only(argc, argv, "cehnv", options, &optionIndex);

        // no more options
        if(c == -1){
            if(argc > optind)
                inputPattern = argv[optind];
            if(argc > optind + 1)
                inputFileName = argv[optind+1];
            break;
        }

        switch(c){
            case '?':
                return -1; // invalid argument
            case 'h':
                printHelp();
                return 0;
            default:
                break;
        }
    }

    if(inputPattern.isEmpty()){
        err << "Error: no input pattern given\n";
        printHelp();
        return -1;
    }

    if(inputFileName.isEmpty()){
        err << "Error: no input file given\n";
        printHelp();
        return -1;
    }

    // select input line format based on the input pattern given
    QString inputFormat;
    if(inputPattern.startsWith("InChI=") || inputPattern[0].isDigit()){
        inputFormat = "inchi";
    }
    else{
        inputFormat = "smiles";
    }

    // create input line format
    chemkit::LineFormat *patternFormat = chemkit::LineFormat::create(inputFormat.toStdString());
    if(!patternFormat){
        err << "Error: failed to create line format\n";
        return -1;
    }

    // read input molecule
    chemkit::Molecule *patternMolecule = patternFormat->read(inputPattern.toStdString());
    if(!patternMolecule){
        err << "Error: failed to read pattern molecule: " << patternFormat->errorString().c_str() << "\n";
        delete patternFormat;
        return -1;
    }

    // read input file
    chemkit::MoleculeFile inputFile(inputFileName.toStdString());
    if(!inputFile.read()){
        err << "Error: failed to read input file: " << inputFile.errorString().c_str() << "\n";
        return -1;
    }

    chemkit::MoleculeFile outputFile;

    int flags = 0;

    if(COMPOSITION_FLAG){
        flags |= chemkit::Molecule::CompareAtomsOnly;
    }

    foreach(chemkit::Molecule *molecule, inputFile.molecules()){
        bool match = false;

        if(EXACT_MATCH_FLAG){
            match = patternMolecule->equals(molecule, flags);
        }
        else{
            match = patternMolecule->isSubstructureOf(molecule, flags);
        }

        if((match && !INVERT_MATCH_FLAG) || (!match && INVERT_MATCH_FLAG)){
            inputFile.removeMolecule(molecule);
            outputFile.addMolecule(molecule);
        }
    }

    if(NAMES_ONLY_FLAG){
        foreach(const chemkit::Molecule *molecule, outputFile.molecules()){
            out << molecule->name().c_str() << "\n";
        }
    }
    else{
        bool ok = outputFile.write(std::cout, inputFile.formatName());
        if(!ok){
            err << "Error: failed to write output file: " << outputFile.errorString().c_str() << "\n";
            return -1;
        }
    }

    delete patternFormat;
    delete patternMolecule;

    return 0;
}
