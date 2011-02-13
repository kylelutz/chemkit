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

#ifndef CHEMKIT_PERIODICTABLEDIALOG_H
#define CHEMKIT_PERIODICTABLEDIALOG_H

#include "widgets.h"

#include <QtGui>

#include <chemkit/element.h>

namespace chemkit {

class PeriodicTableDialogPrivate;

class CHEMKIT_WIDGETS_EXPORT PeriodicTableDialog : public QDialog
{
    Q_OBJECT

    public:
        // construction and destruction
        PeriodicTableDialog(QWidget *parent = 0);
        ~PeriodicTableDialog();

        // properties
        Element element() const;

        // static methods
        static Element getElement(QWidget *parent = 0, const QString &caption = QString());

    private slots:
        void elementClicked(const chemkit::Element &element);

    private:
        PeriodicTableDialogPrivate* const d;
};

} // end chemkit namespace

#endif // CHEMKIT_PERIODICTABLEDIALOG_H
