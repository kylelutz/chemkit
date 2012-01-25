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

#include "periodictabledialog.h"
#include "ui_periodictabledialog.h"

#include "periodictablewidget.h"

namespace chemkit {

// === PeriodicTableDialogPrivate ========================================== //
class PeriodicTableDialogPrivate
{
public:
    Element element;
    bool closeOnClick;
    Ui::PeriodicTableDialog *ui;
};

// === PeriodicTableDialog ================================================= //
/// \class PeriodicTableDialog periodictabledialog.h chemkit/periodictabledialog.h
/// \ingroup chemkit-gui
/// \brief The PeriodicTableDialog class displays the periodic table
///        in a dialog.
///
/// \see PeriodicTableWidget

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new periodic table dialog.
PeriodicTableDialog::PeriodicTableDialog(QWidget *parent)
    : QDialog(parent),
      d(new PeriodicTableDialogPrivate)
{
    d->ui = new Ui::PeriodicTableDialog;
    d->ui->setupUi(this);

    d->closeOnClick = false;

    PeriodicTableWidget *widget = new PeriodicTableWidget(this);
    d->ui->layout->addWidget(widget);

    connect(widget, SIGNAL(elementClicked(chemkit::Element)), SLOT(elementClicked(chemkit::Element)));
}

/// Destroys the periodic table dialog.
PeriodicTableDialog::~PeriodicTableDialog()
{
    delete d->ui;
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the last element on the periodic table that was selected.
Element PeriodicTableDialog::element() const
{
    return d->element;
}

// --- Static Methods ------------------------------------------------------ //
/// Prompts the user to select an element and returns it.
Element PeriodicTableDialog::getElement(QWidget *parent, const QString &caption)
{
    PeriodicTableDialog dialog(parent);
    dialog.setWindowTitle(caption);
    dialog.d->closeOnClick = true;
    dialog.exec();

    return dialog.element();
}

// --- Slots --------------------------------------------------------------- //
void PeriodicTableDialog::elementClicked(const chemkit::Element &element)
{
    d->element = element;

    if(d->closeOnClick){
        close();
    }
}

} // end chemkit namespace
