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

#include <chemkit/chemkit.h>
#include <chemkit/moleculefile.h>

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
    chemkit::io::MoleculeFile inputFile(inputFileName.toStdString());
    inputFile.setFormat(inputFormatName.toStdString());
    if(!inputFile.read()){
        err << "Error: failed to read input file: " << inputFile.errorString().c_str() << "\n";
        return -1;
    }

    // write output
    if(!inputFile.write(outputFileName.toStdString(), outputFormatName.toStdString())){
        err << "Error: failed to write output file: " << inputFile.errorString().c_str() << "\n";
        return -1;
    }

    return 0;
}
