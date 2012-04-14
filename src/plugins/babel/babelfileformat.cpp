/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include "babelfileformat.h"

#include <sstream>

#include <chemkit/moleculefile.h>

#include <QProcess>

BabelFileFormat::BabelFileFormat()
    : chemkit::MoleculeFileFormat("babel")
{
}

bool BabelFileFormat::read(std::istream &input, chemkit::MoleculeFile *file)
{
    // get input format to use
    std::string format = option("format").toString();
    if(format.empty()){
        setErrorString("No format set for Babel conversion.");
        return false;
    }

    // setup babel arguments
    QStringList arguments;
    arguments.append(QString("-i") + format.c_str());
    arguments.append("-");
    arguments.append("-ocml");
    arguments.append("-");

    // create and start the babel process
    QProcess babel;
    babel.start("babel", arguments);
    if(!babel.waitForStarted()){
        setErrorString("Failed to start Babel process.");
        return false;
    }

    // write input data to babel via stdin
    while(!input.eof()){
        char buffer[1024];
        input.read(buffer, sizeof(buffer));
        babel.write(buffer, input.gcount());
    }

    babel.closeWriteChannel();

    // wait until the babel process is finished
    if(!babel.waitForFinished()){
        setErrorString("Babel process never finished.");
        return false;
    }

    // check babel's exit status
    if(babel.exitCode() != QProcess::NormalExit){
        setErrorString("Babel process crashed.");
        return false;
    }

    // read output data to string buffer
    std::stringstream buffer;
    buffer << babel.readAll().constData();

    // parse cml output file
    chemkit::MoleculeFile outputFile;
    outputFile.setFormat("cml");
    bool ok = outputFile.read(buffer);
    if(!ok){
        setErrorString("Failed to parse Babel's CML output: " + outputFile.errorString());
        return false;
    }

    // add molecules to file
    foreach(const boost::shared_ptr<chemkit::Molecule> &molecule, outputFile.molecules()){
        file->addMolecule(molecule);
    }

    return true;
}

bool BabelFileFormat::write(const chemkit::MoleculeFile *file, std::ostream &output)
{
    // get output format to use
    std::string format = option("format").toString();
    if(format.empty()){
        setErrorString("No format set for Babel conversion.");
        return false;
    }

    // write the file to a buffer in the CML format
    std::ostringstream inputBuffer;
    bool ok = const_cast<chemkit::MoleculeFile *>(file)->write(inputBuffer, "cml");
    if(!ok){
        setErrorString("Failed to write CML data for Babel conversion.");
        return false;
    }

    // setup babel arguments
    QStringList arguments;
    arguments.append("-icml");
    arguments.append("-");
    arguments.append(QString("-o") + format.c_str());
    arguments.append("-");

    // create and start the babel process
    QProcess babel;
    babel.start("babel", arguments);
    if(!babel.waitForStarted()){
        setErrorString("Failed to start Babel process.");
        return false;
    }

    // write the CML file data via stdin
    std::string inputData = inputBuffer.str();
    babel.write(inputData.c_str(), inputData.size());
    babel.closeWriteChannel();

    // wait until the babel process is finished
    if(!babel.waitForFinished()){
        setErrorString("Babel process never finished.");
        return false;
    }

    // set output from babel's output
    output << babel.readAll().constData();

    return true;
}
