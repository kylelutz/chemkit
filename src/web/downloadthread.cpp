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

#include "downloadthread.h"

namespace chemkit {

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

} // end chemkit namespace
