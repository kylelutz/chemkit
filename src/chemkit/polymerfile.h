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

#ifndef CHEMKIT_POLYMERFILE_H
#define CHEMKIT_POLYMERFILE_H

#include "chemkit.h"

#include <string>
#include <vector>

#include <QIODevice>

namespace chemkit {

class Polymer;
class PolymerFileFormat;
class PolymerFilePrivate;

class CHEMKIT_EXPORT PolymerFile
{
    public:
        // construction and destruction
        PolymerFile();
        PolymerFile(const std::string &fileName);
        ~PolymerFile();

        // properties
        void setFileName(const std::string &fileName);
        std::string fileName() const;
        void setFormat(PolymerFileFormat *format);
        bool setFormat(const std::string &name);
        PolymerFileFormat* format() const;
        std::string formatName() const;
        int size() const;
        bool isEmpty() const;

        // file contents
        void addPolymer(Polymer *polymer);
        bool removePolymer(Polymer *polymer);
        bool deletePolymer(Polymer *polymer);
        Polymer* polymer(int index = 0) const;
        std::vector<Polymer *> polymers() const;
        int polymerCount() const;
        bool contains(const Polymer *polymer) const;
        void clear();

        // input and output
        bool read();
        bool read(const std::string &fileName);
        bool read(const std::string &fileName, const std::string &format);
        bool read(QIODevice *iodev, const std::string &format);
        bool write();
        bool write(const std::string &fileName);
        bool write(const std::string &fileName, const std::string &format);
        bool write(QIODevice *iodev);
        bool write(QIODevice *iodev, const std::string &format);

        // error handling
        std::string errorString() const;

        // static methods
        static std::vector<std::string> formats();

    private:
        void setErrorString(const std::string &errorString);

    private:
        PolymerFilePrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_POLYMERFILE_H
