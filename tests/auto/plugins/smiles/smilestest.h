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

#ifndef SMILESTEST_H
#define SMILESTEST_H

#include <QtTest>

#include <chemkit/molecule.h>

class SmilesTest : public QObject
{
    Q_OBJECT

    private slots:
        void initTestCase();

        // molecule tests
        void acenaphthylene();
        void aceticAcid();
        void adenine();
        void alanine();
        void ampicillin();
        void anthracene();
        void anthraquinone();
        void arsabenzene();
        void arsole();
        void aspirin();
        void aziridine();
        void azulene();
        void benzene();
        void benzofuran();
        void benzofurazan();
        void benzyne();
        void binol();
        void biphenyl();
        void biphenylene();
        void biperiden();
        void borinine();
        void borole();
        void buckminsterfullerene();
        void butene();
        void caffeine();
        void camphor();
        void carbazole();
        void cholesterol();
        void chrysene();
        void cinnoline();
        void colchicine();
        void copperSulfate();
        void corannulene();
        void coronene();
        void cubane();
        void cyanide();
        void cytosine();
        void decalin();
        void dibenzofuran();
        void dichloroethene();
        void dihydrogen();
        void dinitrogen();
        void ethane();
        void fluorenone();
        void folate();
        void furan();
        void furazan();
        void glucose();
        void guanine();
        void heavyWater();
        void histidine();
        void hydride();
        void hydronium();
        void ibuprofen();
        void indazole();
        void indene();
        void indole();
        void indolizine();
        void ipratropium();
        void isobutane();
        void isoindene();
        void isoindole();
        void melatonin();
        void naphthalene();
        void nicotine();
        void nitrobenzene();
        void ovalene();
        void oxazole();
        void pentacene();
        void pentalene();
        void perylene();
        void phenanthrene();
        void phenothiazine();
        void phenoxazine();
        void phosphole();
        void phosphorine();
        void phthalimide();
        void porphin();
        void proline();
        void proton();
        void purine();
        void pyranium();
        void pyrazole();
        void pyrene();
        void pyridazine();
        void pyridine();
        void pyrimidine();
        void pyrrole();
        void quinoline();
        void rhodizonicAcid();
        void selenophene();
        void sodiumChloride();
        void stilbene();
        void sulfurHexafluoride();
        void taxol();
        void tetralin();
        void tetraphenylene();
        void thiamin();
        void thiirane();
        void thiophene();
        void thujone();
        void thymine();
        void triazole();
        void triphenylene();
        void tropone();
        void tryptophan();
        void uracil();
        void vanillin();

        // feature tests
        void addHydrogens();
        void isotope();
        void kekulize();
        void quadrupleBond();

        // invalid tests
        void extraParenthesis();
        void invalidAtom();
        void wildcardAtom();

        // file tests
        void herg();
        void cox2();

    private:
        void COMPARE_SMILES(const chemkit::Molecule *molecule, const std::string &smiles);
};

#endif // SMILESTEST_H
