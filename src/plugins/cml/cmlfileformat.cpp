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

#include "cmlfileformat.h"

#include <QtXml>

namespace {

class CmlHandler : public QXmlDefaultHandler
{
    public:
        CmlHandler(chemkit::ChemicalFile *file);
        ~CmlHandler();

        bool startElement(const QString &namespaceURI, const QString &localName, const QString &qName, const QXmlAttributes &atts);
        bool endElement(const QString &namespaceURI, const QString &localName, const QString &qName);

    private:
        chemkit::ChemicalFile *m_file;
        chemkit::Molecule *m_molecule;
        QHash<QString, chemkit::Atom *> m_atomIds;
};

CmlHandler::CmlHandler(chemkit::ChemicalFile *file)
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
    Q_UNUSED(namespaceURI);
    Q_UNUSED(localName);

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

            chemkit::Float x = atts.value("x3").toDouble();
            chemkit::Float y = atts.value("y3").toDouble();
            chemkit::Float z = atts.value("z3").toDouble();
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
    Q_UNUSED(namespaceURI);
    Q_UNUSED(localName);

    if(qName == "molecule"){
        m_file->addMolecule(m_molecule);
        m_molecule = 0;
        m_atomIds.clear();
    }

    return true;
}

} // end anonymous namespace

CmlFileFormat::CmlFileFormat()
    : chemkit::ChemicalFileFormat("cml")
{
}

CmlFileFormat::~CmlFileFormat()
{
}

bool CmlFileFormat::read(QIODevice *iodev, chemkit::ChemicalFile *file)
{
    QXmlSimpleReader xml;
    QXmlInputSource source(iodev);
    CmlHandler handler(file);
    xml.setContentHandler(&handler);
    xml.setErrorHandler(&handler);

    bool ok = xml.parse(source);
    if(!ok)
        setErrorString(QString("XML Parsing failed: %1").arg(handler.errorString()).toStdString());

    return ok;
}

bool CmlFileFormat::write(const chemkit::ChemicalFile *file, QIODevice *iodev)
{
    QXmlStreamWriter stream(iodev);
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

    return true;
}
