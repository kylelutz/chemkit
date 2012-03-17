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

#ifndef CHEMKIT_LINEFORMAT_H
#define CHEMKIT_LINEFORMAT_H

#include "chemkit.h"

#include <string>
#include <vector>

#include "plugin.h"
#include "variant.h"

namespace chemkit {

class Molecule;
class LineFormatPrivate;

class CHEMKIT_EXPORT LineFormat
{
public:
    // construction and destruction
    virtual ~LineFormat();

    // properties
    std::string name() const;

    // options
    void setOption(const std::string &name, const Variant &value);
    Variant option(const std::string &name) const;

    // input and output
    virtual Molecule* read(const std::string &formula);
    virtual std::string write(const Molecule *molecule);

    // error handling
    std::string errorString() const;

    // static methods
    static LineFormat* create(const std::string &name);
    static std::vector<std::string> formats();
    static std::string convert(const std::string &formula,
                               const std::string &inputFormatName,
                               const std::string &outputFormatName);

protected:
    LineFormat(const std::string &name);
    void setErrorString(const std::string &error);
    virtual Variant defaultOption(const std::string &name) const;

private:
    LineFormatPrivate* const d;
};

} // end chemkit namespace

/// Registers a line format with \p name.
#define CHEMKIT_REGISTER_LINE_FORMAT(name, className) \
    CHEMKIT_REGISTER_PLUGIN_CLASS(name, chemkit::LineFormat, className)

#endif // CHEMKIT_LINEFORMAT_H
