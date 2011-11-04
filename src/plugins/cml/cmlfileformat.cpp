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

#include "cmlfileformat.h"

#include <QtXml>

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/molecule.h>
#include <chemkit/moleculefile.h>

namespace {

class CmlHandler : public QXmlDefaultHandler
{
public:
    CmlHandler(chemkit::MoleculeFile *file);
    ~CmlHandler();

    bool startElement(const QString &namespaceURI, const QString &localName, const QString &qName, const QXmlAttributes &atts);
    bool endElement(const QString &namespaceURI, const QString &localName, const QString &qName);

private:
    chemkit::MoleculeFile *m_file;
    chemkit::Molecule *m_molecule;
    QHash<QString, chemkit::Atom *> m_atomIds;
};

CmlHandler::CmlHandler(chemkit::MoleculeFile *file)
    : QXmlDefaultHandler(),
      m_file(file),
      m_molecule(0)
{
}

CmlHandler::~CmlHandler()
{
}

bool CmlHandler::startElement(const QString &namespaceURI, const QString &localName, const QString &qName, const QXmlAttributes &atts)
{
    CHEMKIT_UNUSED(namespaceURI);
    CHEMKIT_UNUSED(localName);

    if(qName == "molecule"){
        m_molecule = new chemkit::Molecule();
    }
    else if(qName == "atom" && m_molecule){
        QString symbol = atts.value("elementType");

        if(!symbol.isEmpty()){
            chemkit::Atom *atom = m_molecule->addAtom(symbol.toStdString());
            if(!atom){
                qDebug() << "invalid atom symbol: " << symbol;
                return true;
            }

            QString id = atts.value("id");
            if(!id.isEmpty())
                m_atomIds[id] = atom;

            chemkit::Real x = atts.value("x3").toDouble();
            chemkit::Real y = atts.value("y3").toDouble();
            chemkit::Real z = atts.value("z3").toDouble();
            atom->setPosition(x, y, z);
        }
    }
    else if(qName == "bond" && m_molecule){
        QString atomRefs = atts.value("atomRefs2");
        if(!atomRefs.isEmpty()){
            QStringList atomsIds = atomRefs.split(' ', QString::SkipEmptyParts);
            if(atomsIds.size() != 2){
                qDebug() << "atomRefs size != 2";
                return true;
            }

            chemkit::Atom *atom1 = m_atomIds[atomsIds[0]];
            chemkit::Atom *atom2 = m_atomIds[atomsIds[1]];

            if(!atom1){
                qDebug() << "invalid atom ref: " << atomsIds[0];
                return true;
            }
            if(!atom2){
                qDebug() << "invalid atom ref: " << atomsIds[1];
                return true;
            }

            int bondOrder;
            QString order = atts.value("order");
            if(order.isEmpty()){
                bondOrder = chemkit::Bond::Single;
            }
            else{
                bondOrder = order.toInt();

                if(bondOrder == 0){
                    if(order == "S")
                        bondOrder = chemkit::Bond::Single;
                    else if(order == "D")
                        bondOrder = chemkit::Bond::Double;
                    else if(order == "T")
                        bondOrder = chemkit::Bond::Triple;
                    else if(order == "A")
                        bondOrder = chemkit::Bond::Single;
                }
            }

            m_molecule->addBond(atom1, atom2, bondOrder);
        }
    }

    return true;
}

bool CmlHandler::endElement(const QString &namespaceURI, const QString &localName, const QString &qName)
{
    CHEMKIT_UNUSED(namespaceURI);
    CHEMKIT_UNUSED(localName);

    if(qName == "molecule"){
        m_file->addMolecule(m_molecule);
        m_molecule = 0;
        m_atomIds.clear();
    }

    return true;
}

} // end anonymous namespace

CmlFileFormat::CmlFileFormat()
    : chemkit::MoleculeFileFormat("cml")
{
}

CmlFileFormat::~CmlFileFormat()
{
}

bool CmlFileFormat::read(std::istream &input, chemkit::MoleculeFile *file)
{
    QByteArray data;
    while(!input.eof()){
        data += input.get();
    }
    data.chop(1);

    QBuffer buffer;
    buffer.setData(data);
    buffer.open(QBuffer::ReadOnly);

    QXmlSimpleReader xml;
    QXmlInputSource source(&buffer);
    CmlHandler handler(file);
    xml.setContentHandler(&handler);
    xml.setErrorHandler(&handler);

    bool ok = xml.parse(source);
    if(!ok)
        setErrorString(QString("XML Parsing failed: %1").arg(handler.errorString()).toStdString());

    return ok;
}

bool CmlFileFormat::write(const chemkit::MoleculeFile *file, std::ostream &output)
{
    QBuffer buffer;
    buffer.open(QBuffer::WriteOnly);

    QXmlStreamWriter stream(&buffer);
    stream.setAutoFormatting(true);

    stream.writeStartDocument();

    foreach(const chemkit::Molecule *molecule, file->molecules()){
        stream.writeStartElement("molecule");

        stream.writeTextElement("name", molecule->name().c_str());

        stream.writeStartElement("atomArray");
        foreach(const chemkit::Atom *atom, molecule->atoms()){
            stream.writeStartElement("atom");
            stream.writeAttribute("id", QString("a%1").arg(atom->index()+1));
            stream.writeAttribute("elementType", atom->symbol().c_str());
            stream.writeAttribute("x3", QString::number(atom->x()));
            stream.writeAttribute("y3", QString::number(atom->y()));
            stream.writeAttribute("z3", QString::number(atom->z()));
            stream.writeEndElement();
        }
        stream.writeEndElement();

        stream.writeStartElement("bondArray");
        foreach(const chemkit::Bond *bond, molecule->bonds()){
            stream.writeStartElement("bond");
            stream.writeAttribute("atomRefs2", QString("a%1 a%2").arg(bond->atom1()->index()+1).arg(bond->atom2()->index()+1));
            stream.writeAttribute("order", QString::number(bond->order()));
            stream.writeEndElement();
        }
        stream.writeEndElement();

        stream.writeEndElement();
    }

    stream.writeEndDocument();

    QByteArray data = buffer.data();
    output.write(data.constData(), data.size());

    return true;
}
