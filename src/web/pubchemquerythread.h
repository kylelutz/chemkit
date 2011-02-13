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

#ifndef CHEMKIT_PUBCHEMQUERYTHREAD_H
#define CHEMKIT_PUBCHEMQUERYTHREAD_H

#include <QtCore>
#include <QtNetwork>

namespace chemkit {

class PubChemQuery;

class PubChemQueryThread : public QThread
{
    Q_OBJECT

    public:
        // construction and destruction
        PubChemQueryThread(const PubChemQuery &query);
        ~PubChemQueryThread();

        // properties
        QByteArray response() const;

        // thread
        void run();

    private slots:
        void replyFinished(QNetworkReply *reply);

    private:
        void poll(const QString &id);

    private:
        QUrl m_url;
        QByteArray m_request;
        QByteArray m_response;
        QNetworkReply *m_reply;
        QNetworkAccessManager *m_manager;
};

} // end chemkit namespace

#endif // CHEMKIT_PUBCHEMQUERYTHREAD_H
