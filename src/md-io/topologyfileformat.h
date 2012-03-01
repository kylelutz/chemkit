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

#ifndef CHEMKIT_TOPOLOGYFILEFORMAT_H
#define CHEMKIT_TOPOLOGYFILEFORMAT_H

#include "md-io.h"

#include <string>
#include <vector>

#include <chemkit/plugin.h>

namespace chemkit {

class TopologyFile;
class TopologyFileFormatPrivate;

class CHEMKIT_MD_IO_EXPORT TopologyFileFormat
{
public:
    // construction and destruction
    virtual ~TopologyFileFormat();

    // properties
    std::string name() const;

    // input and output
    virtual bool read(std::istream &input, TopologyFile *file);
    virtual bool write(const TopologyFile *file, std::ostream &output);

    // error handling
    std::string errorString() const;

    // static methods
    static TopologyFileFormat* create(const std::string &name);
    static std::vector<std::string> formats();

protected:
    TopologyFileFormat(const std::string &name);
    void setErrorString(const std::string &error);

private:
    TopologyFileFormatPrivate* const d;
};

} // end chemkit namespace

/// Registers a topology file format with \p name.
#define CHEMKIT_REGISTER_TOPOLOGY_FILE_FORMAT(name, className) \
    CHEMKIT_REGISTER_PLUGIN_CLASS(name, chemkit::TopologyFileFormat, className)

#endif // CHEMKIT_TOPOLOGYFILEFORMAT_H
