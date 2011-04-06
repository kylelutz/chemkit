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

#include "pubchemquery.h"

namespace chemkit {

// === PubChemQuery ======================================================== //
/// \class PubChemQuery pubchemquery.h
/// \ingroup chemkit-web
/// \internal
/// \brief The PubChemQuery class creates querys to be sent to the
///        PubChem Power User Gateway (PUG).

// --- Construction and Destruction ---------------------------------------- //
PubChemQuery::PubChemQuery()
{
}

PubChemQuery::~PubChemQuery()
{
}

// --- Properties ---------------------------------------------------------- //
/// Sets the data for the query.
void PubChemQuery::setData(const QByteArray &data)
{
    m_data = data;
}

/// Returns the data for the query.
QByteArray PubChemQuery::data() const
{
    return m_data;
}

// --- Static Methods ------------------------------------------------------ //
PubChemQuery PubChemQuery::downloadQuery(const QStringList &cids, const QString &format)
{
    PubChemQuery query;

    QString xml = "<PCT-Data>"
                    "<PCT-Data_input>"
                      "<PCT-InputData>"
                        "<PCT-InputData_download>"
                          "<PCT-Download>"
                            "<PCT-Download_uids>"
                              "<PCT-QueryUids>"
                                "<PCT-QueryUids_ids>"
                                  "<PCT-ID-List>"
                                    "<PCT-ID-List_db>pccompound</PCT-ID-List_db>"
                                        "<PCT-ID-List_uids>"
                                          "%1"
                                        "</PCT-ID-List_uids>"
                                      "</PCT-ID-List>"
                                    "</PCT-QueryUids_ids>"
                                  "</PCT-QueryUids>"
                                "</PCT-Download_uids>"
                              "<PCT-Download_format value=\"%2\"/>"
                            "<PCT-Download_compression value=\"none\"/>"
                          "</PCT-Download>"
                        "</PCT-InputData_download>"
                      "</PCT-InputData>"
                    "</PCT-Data_input>"
                  "</PCT-Data>";

    QString idsXml;

    foreach(const QString &id, cids){
        idsXml += QString("<PCT-ID-List_uids_E>%1</PCT-ID-List_uids_E>").arg(id);
    }

    query.setData(xml.arg(idsXml).arg(format).toAscii());

    return query;
}

PubChemQuery PubChemQuery::standardizationQuery(const std::string &formula, const std::string &inputFormat, const std::string &outputFormat)
{
    PubChemQuery query;

    QString xml = "<PCT-Data>"
                    "<PCT-Data_input>"
                      "<PCT-InputData>"
                        "<PCT-InputData_standardize>"
                          "<PCT-Standardize>"
                            "<PCT-Standardize_structure>"
                              "<PCT-Structure>"
                                "<PCT-Structure_structure>"
                                  "<PCT-Structure_structure_string>%1"
                                  "</PCT-Structure_structure_string>"
                                "</PCT-Structure_structure>"
                                "<PCT-Structure_format>"
                                  "<PCT-StructureFormat value=\"%2\"/>"
                                "</PCT-Structure_format>"
                              "</PCT-Structure>"
                            "</PCT-Standardize_structure>"
                            "<PCT-Standardize_oformat>"
                              "<PCT-StructureFormat value=\"%3\"/>"
                            "</PCT-Standardize_oformat>"
                          "</PCT-Standardize>"
                        "</PCT-InputData_standardize>"
                      "</PCT-InputData>"
                    "</PCT-Data_input>"
                  "</PCT-Data>";

    query.setData(xml.arg(formula.c_str()).arg(inputFormat.c_str()).arg(outputFormat.c_str()).toAscii());

    return query;
}

} // end chemkit namespace
