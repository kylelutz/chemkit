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

#include "downloadthread.h"

namespace chemkit {
namespace web {

// === DownloadThread ====================================================== //
/// \class DownloadThread downloadthread.h
/// \ingroup chemkit-web
/// \internal
/// \brief The DownloadThread class provides an interface for
///        downloading data from the web.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new download thread object.
DownloadThread::DownloadThread(const QUrl &url)
    : QThread(),
      m_url(url)
{
    m_manager = 0;
    m_ftp = 0;
    m_reply = 0;
}

/// Destroys the download thread object.
DownloadThread::~DownloadThread()
{
    delete m_ftp;
    delete m_manager;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the downloaded data.
QByteArray DownloadThread::data() const
{
    return m_data;
}

// --- Thread -------------------------------------------------------------- //
/// Starts the download.
///
/// This method should not be called directly. Use QThread::start().
void DownloadThread::run()
{
    if(m_url.scheme() == "ftp"){
        m_dataBuffer.open(QBuffer::ReadWrite);

        m_ftp = new QFtp;
        connect(m_ftp, SIGNAL(done(bool)), SLOT(ftpDone(bool)), Qt::DirectConnection);
        m_ftp->connectToHost(m_url.host());
        m_ftp->login("anonymous", "chemkit");

        QFileInfo fileInfo(m_url.path());
        m_ftp->cd(fileInfo.path() + "/");
        m_ftp->get(fileInfo.fileName(), &m_dataBuffer);
        m_ftp->close();
    }
    else{
        m_manager = new QNetworkAccessManager;
        connect(m_manager, SIGNAL(finished(QNetworkReply*)),
                SLOT(replyFinished(QNetworkReply*)), Qt::DirectConnection);

        m_reply = m_manager->get(QNetworkRequest(m_url));
        connect(m_reply, SIGNAL(error(QNetworkReply::NetworkError)),
                SLOT(replyError(QNetworkReply::NetworkError)), Qt::DirectConnection);
    }

    exec();
}

// --- Static Methods ------------------------------------------------------ //
/// Downloads and returns the data from \p url.
QByteArray DownloadThread::download(const QUrl &url)
{
    DownloadThread thread(url);
    thread.start();
    thread.wait();
    return thread.data();
}

// --- Slots --------------------------------------------------------------- //
void DownloadThread::replyFinished(QNetworkReply *reply)
{
    m_data = reply->readAll();

    reply->deleteLater();
    m_manager->deleteLater();

    exit(0);
}

void DownloadThread::replyError(QNetworkReply::NetworkError error)
{
    Q_UNUSED(error);

    qDebug() << "chemkit-web: DownloadThread Error: " << m_reply->errorString();
}

void DownloadThread::ftpDone(bool error)
{
    if(error){
        qDebug() << "chemkit-web: DownloadThread Ftp Error:" << m_ftp->errorString();
        exit(-1);
    }

    m_data = m_dataBuffer.data();

    exit(0);
}

} // end web namespace
} // end chemkit namespace
