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

#include "periodictablewidget.h"
#include "ui_periodictablewidget.h"

#include <QtGui>

#include <chemkit/element.h>

namespace chemkit {

namespace {

const unsigned char PeriodicTable[] = {
    1,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,
    3,    4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   5,   6,   7,   8,   9,  10,
    11,  12,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  13,  14,  15,  16,  17,  18,
    19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,
    37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,
    55,  56,   0,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,  82,  83,  84,  85,  86,
    87,  88,   0, 104, 105, 106, 107, 108, 109,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,
     0,   0,   0,  89,  90,  91,  92,  93,  94,  95,  96,  97,  98,  99, 100, 101, 102, 103
};

const int PeriodicTableRows = 9;
const int PeriodicTableColumns = 18;

} // end anonymous namespace

// === PeriodicTableWidgetPrivate ========================================== //
class PeriodicTableWidgetPrivate
{
public:
    Ui::PeriodicTableWidget *ui;
    QSignalMapper *signalMapper;
};

// === PeriodicTableWidget ================================================= //
/// \class PeriodicTableWidget periodictablewidget.h chemkit/periodictablewidget.h
/// \ingroup chemkit-gui
/// \brief The PeriodicTableWidget class displays the periodic table.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new periodic table widget.
PeriodicTableWidget::PeriodicTableWidget(QWidget *parent)
    : QWidget(parent),
      d(new PeriodicTableWidgetPrivate)
{
    d->ui = new Ui::PeriodicTableWidget;
    d->ui->setupUi(this);

    d->signalMapper = new QSignalMapper(this);
    connect(d->signalMapper, SIGNAL(mapped(int)), SLOT(buttonClicked(int)));

    for(int row = 0; row < PeriodicTableRows; row++){
        for(int column = 0; column < PeriodicTableColumns; column++){
            int atomicNumber = PeriodicTable[row * PeriodicTableColumns + column];
            if(atomicNumber == 0){
                continue;
            }

            Element element(atomicNumber);

            QPushButton *button = new QPushButton(element.symbol().c_str());
            button->setMinimumWidth(30);

            connect(button, SIGNAL(clicked()), d->signalMapper, SLOT(map()));
            d->signalMapper->setMapping(button, atomicNumber);

            if(row >= 7){
                if(row == 7){
                    // add spacing item between main elements and the lanthanindes
                    d->ui->gridLayout->addItem(new QSpacerItem(0, 10), row, column);
                }

                d->ui->gridLayout->addWidget(button, row + 1, column);
            }
            else{
                d->ui->gridLayout->addWidget(button, row, column);
            }
        }
    }
}

/// Destroys the periodic table widget.
PeriodicTableWidget::~PeriodicTableWidget()
{
    delete d->ui;
    delete d;
}

// --- Signals ------------------------------------------------------------- //
/// \fn void PeriodicTableWidget::elementClicked(const chemkit::Element &element)
///
/// This signal is emitted when an element is clicked.

// --- Slots --------------------------------------------------------------- //
void PeriodicTableWidget::buttonClicked(int atomicNumber)
{
    emit elementClicked(Element(atomicNumber));
}

} // end chemkit namespace
