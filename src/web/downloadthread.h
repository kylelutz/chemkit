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

#ifndef CHEMKIT_DOWNLOADTHREAD_H
#define CHEMKIT_DOWNLOADTHREAD_H

#include <QtCore>
#include <QtNetwork>

namespace chemkit {

class DownloadThread : public QThread
{
    Q_OBJECT

    public:
        // construction and destruction
        DownloadThread(const QUrl &url);
        ~DownloadThread();

        // properties
        QByteArray data() const;

        // thread
        void run();

        // static methods
        static QByteArray download(const QUrl &url);

    private slots:
        void replyFinished(QNetworkReply *reply);
        void replyError(QNetworkReply::NetworkError error);
        void ftpDone(bool error);

    private:
        QUrl m_url;
        QFtp *m_ftp;
        QBuffer m_dataBuffer;
        QNetworkReply *m_reply;
        QNetworkAccessManager *m_manager;
        QByteArray m_data;
};

} // end chemkit namespace

#endif // CHEMKIT_DOWNLOADTHREAD_H
