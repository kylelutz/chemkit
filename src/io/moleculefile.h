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

#ifndef CHEMKIT_MOLECULEFILE_H
#define CHEMKIT_MOLECULEFILE_H

#include "io.h"

#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/range/iterator_range.hpp>

#include "genericfile.h"
#include "moleculefileformat.h"

namespace chemkit {

class Molecule;
class MoleculeFilePrivate;

class CHEMKIT_IO_EXPORT MoleculeFile : public GenericFile<MoleculeFile, MoleculeFileFormat>
{
public:
    // typedefs
    typedef boost::iterator_range<std::vector<boost::shared_ptr<Molecule> >::const_iterator> MoleculeRange;

    // construction and destruction
    MoleculeFile();
    MoleculeFile(const std::string &fileName);
    virtual ~MoleculeFile();

    // properties
    size_t size() const;
    bool isEmpty() const;

    // file contents
    void addMolecule(const boost::shared_ptr<Molecule> &molecule);
    bool removeMolecule(const boost::shared_ptr<Molecule> &molecule);
    MoleculeRange molecules() const;
    size_t moleculeCount() const;
    boost::shared_ptr<Molecule> molecule(size_t index = 0) const;
    boost::shared_ptr<Molecule> molecule(const std::string &name) const;
    bool contains(const boost::shared_ptr<Molecule> &molecule) const;
    void clear();

    // static methods
    static boost::shared_ptr<Molecule> quickRead(const std::string &fileName);
    static void quickWrite(const Molecule *molecule, const std::string &fileName);

private:
    MoleculeFilePrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_MOLECULEFILE_H
