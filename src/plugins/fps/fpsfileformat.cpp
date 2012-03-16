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

#include "fpsfileformat.h"

#include <ctime>
#include <iomanip>

#include <chemkit/foreach.h>
#include <chemkit/fingerprint.h>
#include <chemkit/moleculefile.h>

namespace {

// Output iterator which writes bitset blocks to an output stream.
struct FingerprintWriter : public std::iterator<std::output_iterator_tag, void, void, void, void>
{
    typedef chemkit::Bitset::block_type value_type;

    FingerprintWriter(std::ostream &output)
        : m_output(output)
    {
    }

    FingerprintWriter(const FingerprintWriter &other)
        : m_output(other.m_output)
    {
    }

    FingerprintWriter& operator=(const value_type &value)
    {
        const unsigned char *bytes =
            reinterpret_cast<const unsigned char *>(&value);

        for(size_t i = 0; i < sizeof(value_type); i++){
            m_output << std::hex
                     << std::setw(2)
                     << std::setfill('0')
                     << static_cast<int>(bytes[i]);
        }

        return *this;
    }

    FingerprintWriter& operator*()
    {
        return *this;
    }

    FingerprintWriter& operator++()
    {
        return *this;
    }

    FingerprintWriter& operator++(int)
    {
        return *this;
    }

    std::ostream &m_output;
};

} // end anonymous namespace

FpsFileFormat::FpsFileFormat()
    : chemkit::MoleculeFileFormat("fps")
{
}

// Writes the fingerprint values for each molecule in the file
// to the output stream.
//
// Reference:
//   http://code.google.com/p/chem-fingerprints/wiki/FPS
bool FpsFileFormat::write(const chemkit::MoleculeFile *file, std::ostream &output)
{
    std::string fingerprintName = option("fingerprint").toString();

    // create fingerprint format
    boost::scoped_ptr<chemkit::Fingerprint> fingerprint(chemkit::Fingerprint::create(fingerprintName));
    if(!fingerprint){
        setErrorString("Fingerprint format '" + fingerprintName + "' is not supported.");
        return false;
    }

    // write header
    output << "#FPS1" << std::endl;
    output << "#num_bits=" << fingerprint->size() << std::endl;
    output << "#type=" << fingerprintTypeString() << std::endl;
    output << "#software=chemkit/" << CHEMKIT_VERSION_STRING << std::endl;
    output << "#date=" << dateTimeString() << std::endl;

    // create output writer
    FingerprintWriter writer(output);

    // write each molecule's fingerprint and identifier
    foreach(const boost::shared_ptr<chemkit::Molecule> &molecule, file->molecules()){
        // write fingerprint
        boost::to_block_range(fingerprint->value(molecule.get()), writer);

        // write identifier (name if available, else molecular formula)
        std::string identifier = molecule->name();
        if(identifier.empty()){
            identifier = molecule->formula();
        }

        output << "\t" << identifier << "\n";
    }

    return true;
}

// Returns the default value for the option specified by name.
chemkit::Variant FpsFileFormat::defaultOption(const std::string &name) const
{
    if(name == "fingerprint"){
        return "fp2";
    }

    return chemkit::Variant();
}

// Returns the date and time as a string formatted according to
// the FPS file format standard.
std::string FpsFileFormat::dateTimeString() const
{
    time_t rawtime = time(NULL);
    struct tm *timeinfo = gmtime(&rawtime);

    char buffer[80];
    strftime(buffer, 80, "%Y-%m-%dT%H:%M:%S", timeinfo);

    return buffer;
}

// Returns a string containing the type of the fingerprint.
std::string FpsFileFormat::fingerprintTypeString() const
{
    std::string fingerprint = option("fingerprint").toString();

    if(fingerprint == "fp2"){
        return "chemkit-FP2/1";
    }
    else if(fingerprint == "pubchem"){
        return "PubChem/1";
    }
    else{
        return fingerprint + "/1";
    }
}
