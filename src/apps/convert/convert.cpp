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
#include <chemkit/chemicalfile.h>

struct option options[] = {
    {"input-format", required_argument, 0, 'i'},
    {"output-format", required_argument, 0, 'o'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0},
};

void printHelp()
{
    QTextStream out(stdout);

    out << QString("Usage: %1 [OPTIONS] inputFile outputFile\n").arg(qAppName());
    out << "\n";
    out << "Converts a chemical input file to a new file with\n";
    out << "a different file format.\n";
    out << "\n";
    out << "Options:\n";
    out << "  -i, --input-format=FORMAT    Sets the input format.\n";
    out << "  -o, --output-format=FORMAT   Sets the output format.\n";
    out << "  -h, --help                   Shows this help message.\n";
    out << "\n";

    out.flush();
}

int main(int argc, char *argv[])
{
    QCoreApplication app(argc, argv);
    app.setApplicationName(argv[0]);

    QTextStream out(stdout);
    QTextStream err(stderr);

    QString inputFileName;
    QString inputFormatName;
    QString outputFileName;
    QString outputFormatName;

    // parse options
    for(;;){
        int optionIndex = 0;

        int c = getopt_long_only(argc, argv, "ioh", options, &optionIndex);

        // no more options
        if(c == -1){
            if(argc > optind)
                inputFileName = argv[optind];
            if(argc > optind + 1)
                outputFileName = argv[optind+1];
            break;
        }

        switch(c){
            case '?':
                return -1; // invalid argument
            case 'i':
                inputFormatName = optarg;
                break;
            case 'o':
                outputFormatName = optarg;
                break;
            case 'h':
                printHelp();
                return 0;
            default:
                break;
        }
    }

    if(inputFileName.isEmpty()){
        err << "Error: no input file given\n";
        printHelp();
        return -1;
    }
    if(outputFileName.isEmpty()){
        err << "Error: no output file given\n";
        printHelp();
        return -1;
    }

    if(inputFormatName.isEmpty()){
        inputFormatName = QFileInfo(inputFileName).suffix();
    }
    if(outputFormatName.isEmpty()){
        outputFormatName = QFileInfo(outputFileName).suffix();
    }

    // read input
    chemkit::ChemicalFile inputFile(inputFileName);
    inputFile.setFormat(inputFormatName);
    if(!inputFile.read()){
        err << "Error: failed to read input file: " << inputFile.errorString() << "\n";
        return -1;
    }

    // write output
    if(!inputFile.write(outputFileName, outputFormatName)){
        err << "Error: failed to write output file: " << inputFile.errorString() << "\n";
        return -1;
    }

    return 0;
}
