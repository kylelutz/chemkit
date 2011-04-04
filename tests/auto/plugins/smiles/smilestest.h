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
