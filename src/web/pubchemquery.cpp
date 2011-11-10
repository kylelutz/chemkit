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

#include "pubchemquery.h"

namespace chemkit {
namespace web {

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

} // end web namespace
} // end chemkit namespace
