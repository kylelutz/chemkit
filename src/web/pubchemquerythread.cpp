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

#include "pubchemquerythread.h"

#include <QtXml>

#include "pubchemquery.h"

namespace chemkit {

// === PubChemQueryThread ================================================== //
/// \class PubChemQueryThread PubChemQueryThread.h
/// \ingroup chemkit-web
/// \internal
/// \brief The PubChemQueryThread class downloads data using the
///        %PubChem Power User Gateway (PUG).

// --- Construction and Destruction ---------------------------------------- //
PubChemQueryThread::PubChemQueryThread(const PubChemQuery &query)
    : QThread()
{
    m_request = query.data();
    m_manager = 0;
    m_url = QString("http://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi");
}

PubChemQueryThread::~PubChemQueryThread()
{
    delete m_reply;
    delete m_manager;
}

// --- Properties ---------------------------------------------------------- //
QByteArray PubChemQueryThread::response() const
{
    return m_response;
}

// --- Thread -------------------------------------------------------------- //
void PubChemQueryThread::run()
{
    m_manager = new QNetworkAccessManager;
    connect(m_manager, SIGNAL(finished(QNetworkReply*)),
            SLOT(replyFinished(QNetworkReply*)), Qt::DirectConnection);

    m_reply = m_manager->post(QNetworkRequest(m_url), m_request);

    exec();
}

// --- Slots --------------------------------------------------------------- //
void PubChemQueryThread::replyFinished(QNetworkReply *reply)
{
    QByteArray replyData = reply->readAll();

    // check to see if the reply contains a request id. If it does
    // the PUG must be polled again.
    QDomDocument document;
    document.setContent(replyData, false);

    QDomNodeList nodes = document.elementsByTagName("PCT-Waiting_reqid");
    if(!nodes.isEmpty()){
        QDomNode node = nodes.at(0);
        QDomElement element = node.toElement();

        QString waitingId = element.text();
        poll(waitingId);
    }
    else{
        m_response = replyData;
        exit(0);
    }
}

// --- Internal Methods ---------------------------------------------------- //
/// Poll the PUG interface and get an update on the status of the
/// request with \p id.
void PubChemQueryThread::poll(const QString &id)
{
    QString xml = "<PCT-Data>"
                    "<PCT-Data_input>"
                      "<PCT-InputData>"
                        "<PCT-InputData_request>"
                          "<PCT-Request>"
                            "<PCT-Request_reqid>%1</PCT-Request_reqid>"
                            "<PCT-Request_type value=\"status\"/>"
                          "</PCT-Request>"
                        "</PCT-InputData_request>"
                      "</PCT-InputData>"
                    "</PCT-Data_input>"
                  "</PCT-Data>";

    QByteArray request = xml.arg(id).toAscii();

    m_reply = m_manager->post(QNetworkRequest(m_url), request);
}

} // end chemkit namespace
