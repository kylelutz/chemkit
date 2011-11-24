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
