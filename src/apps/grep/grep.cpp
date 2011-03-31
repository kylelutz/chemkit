/******************************************************************************
**
** Copyright (C) 2009-2011 Kyle Lutz <kyle.r.lutz@gmail.com>
**
** This file is part of chemkit. For more information see
** <http://www.chemkit.org>.
**
** chemkit is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** chemkit is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with chemkit. If not, see <http://www.gnu.org/licenses/>.
**
******************************************************************************/

#include <QtCore>

#include <getopt.h>

#include <chemkit/chemkit.h>
#include <chemkit/molecule.h>
#include <chemkit/lineformat.h>
#include <chemkit/chemicalfile.h>

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
    chemkit::ChemicalFile inputFile(inputFileName.toStdString());
    if(!inputFile.read()){
        err << "Error: failed to read input file: " << inputFile.errorString().c_str() << "\n";
        return -1;
    }

    chemkit::ChemicalFile outputFile;

    chemkit::Molecule::CompareFlags flags;

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
        QFile stdoutFile;
        stdoutFile.open(stdout, QFile::WriteOnly);
        bool ok = outputFile.write(&stdoutFile, inputFile.formatName());
        stdoutFile.close();

        if(!ok){
            err << "Error: failed to write output file: " << outputFile.errorString().c_str() << "\n";
            return -1;
        }
    }

    delete patternFormat;
    delete patternMolecule;

    return 0;
}
