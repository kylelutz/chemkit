/******************************************************************************
**
** Copyright (C) 2009-2012 Kyle Lutz <kyle.r.lutz@gmail.com>
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

#include <chemkit/plugin.h>
#include <chemkit/foreach.h>

#include <QProcess>

#include "babelfileformat.h"

namespace {

/// The BabelFileFormatInstantiator class handles instantiating a
/// BabelFileFormat instance with the correct format name.
class BabelFileFormatInstantiator
{
public:
    BabelFileFormatInstantiator(const std::string &format)
        : m_format(format)
    {
    }

    BabelFileFormat* operator()()
    {
        BabelFileFormat *format = new BabelFileFormat;
        format->setOption("format", m_format);
        return format;
    }

private:
    std::string m_format;
};

} // end anonymous namespace

class BabelPlugin : public chemkit::Plugin
{
public:
    BabelPlugin()
        : chemkit::Plugin("babel")
    {
        // check for availability of the babel program
        QProcess babel;
        babel.start("babel");
        if(!babel.waitForStarted()){
            // babel is not available so return without
            // registering any of the babel file formats
            return;
        }
        babel.waitForFinished();

        // register the generic babel file format
        CHEMKIT_REGISTER_MOLECULE_FILE_FORMAT("babel", BabelFileFormat);

        // list of formats for which to use babel
        std::vector<std::string> babelFormats;
        babelFormats.push_back("acr"); // ACR Carine ASCII Crystal
        babelFormats.push_back("adf"); // ADF Input
        babelFormats.push_back("adfout"); // ADF Output
        babelFormats.push_back("alc"); // Alchemy
        babelFormats.push_back("arc"); // Accelrys/MSI Biosym/Insight II CAR
        babelFormats.push_back("bgf"); // MSI BGF
        babelFormats.push_back("box"); // Dock 3.5 Box
        babelFormats.push_back("bs"); // Ball and Stick
        babelFormats.push_back("c3d1"); // Chem3D Cartesian 1
        babelFormats.push_back("c3d2"); // Chem3D Cartesian 2
        babelFormats.push_back("cac"); // CAChe MolStruct
        babelFormats.push_back("caccrt"); // Cacao Cartesian
        babelFormats.push_back("cache"); // CAChe MolStruct
        babelFormats.push_back("cacint"); // Cacao Internal
        babelFormats.push_back("car"); // Accelrys/MSI Biosym/Insight II CAR
        babelFormats.push_back("ccc"); // CCC
        babelFormats.push_back("cdx"); // ChemDraw CDX
        babelFormats.push_back("cdxml"); // ChemDraw CDXML
        babelFormats.push_back("cht"); // Chemtool
        babelFormats.push_back("cif"); // Crystallographic Information File
        babelFormats.push_back("ck"); // Chemkin
        babelFormats.push_back("com"); // Gaussian 98/03 Cartesian Input
        babelFormats.push_back("crk2d"); // Chemical Resource Kit 2D
        babelFormats.push_back("crk3d"); // Chemical Resource Kit 3D
        babelFormats.push_back("csr"); // Accelrys/MSI Quanta CSR
        babelFormats.push_back("cssr"); // CSD CSSR
        babelFormats.push_back("ct"); // ChemDraw Connection Table
        babelFormats.push_back("dmol"); // DMol3 coordinates
        babelFormats.push_back("dx"); // OpenDX grid
        babelFormats.push_back("fa"); // FASTA
        babelFormats.push_back("fasta"); // FASTA
        babelFormats.push_back("fch"); // Gaussian checkpoint file
        babelFormats.push_back("fchk"); // Gaussian checkpoint file
        babelFormats.push_back("fck"); // Gaussian checkpoint file
        babelFormats.push_back("fract"); // Free Form Fractional
        babelFormats.push_back("fsa"); // FASTA
        babelFormats.push_back("g03"); // Gaussian98/03 Output
        babelFormats.push_back("g92"); // Gaussian98/03 Output
        babelFormats.push_back("g94"); // Gaussian98/03 Output
        babelFormats.push_back("g98"); // Gaussian98/03 Output
        babelFormats.push_back("gal"); // Gaussian98/03 Output
        babelFormats.push_back("gam"); // GAMESS Output
        babelFormats.push_back("gamin"); // GAMESS Input
        babelFormats.push_back("gamout"); // GAMESS Output
        babelFormats.push_back("gau"); // Gaussian 98/03 Cartesian Input
        babelFormats.push_back("gjc"); // Gaussian 98/03 Input
        babelFormats.push_back("gjf"); // Gaussian 98/03 Input
        babelFormats.push_back("gpr"); // Ghemical
        babelFormats.push_back("gr96"); // GROMOS96
        babelFormats.push_back("gukin"); // GAMESS UK Input
        babelFormats.push_back("gukout"); // GAMESS UK Output
        babelFormats.push_back("gzmat"); // Gaussian Z-matrix Input
        babelFormats.push_back("hin"); // HyperChem HIN
        babelFormats.push_back("ins"); // ShelX
        babelFormats.push_back("jin"); // Jaguar input
        babelFormats.push_back("jout"); // Jaguar output
        babelFormats.push_back("mcif"); // mmCIF
        babelFormats.push_back("mmcif"); // mmCIF
        babelFormats.push_back("mmd"); // MacroModel
        babelFormats.push_back("mmod"); // MacroModel
        babelFormats.push_back("molden"); // Molden
        babelFormats.push_back("moo"); // MOPAC Output
        babelFormats.push_back("mopout"); // MOPAC Output
        babelFormats.push_back("mpd"); // Sybyl descriptor
        babelFormats.push_back("mpqc"); // MPQC output
        babelFormats.push_back("mpqcin"); // MPQC simplified input
        babelFormats.push_back("msi"); // Accelrys MSI text
        babelFormats.push_back("msms"); // MSMS input
        babelFormats.push_back("nw"); // NWChem input
        babelFormats.push_back("nwo"); // NWChem output
        babelFormats.push_back("outmol"); // DMol3 coordinates
        babelFormats.push_back("pcm"); // PCModel
        babelFormats.push_back("png"); // PNG (embedded)
        babelFormats.push_back("pov"); // POV-Ray input
        babelFormats.push_back("pqs"); // Parallel Quantum Solutions
        babelFormats.push_back("prep"); // Amber Prep
        babelFormats.push_back("qcin"); // Q-Chem input
        babelFormats.push_back("qcout"); // Q-Chem output
        babelFormats.push_back("res"); // ShelX
        babelFormats.push_back("t41"); // ADF Tape41
        babelFormats.push_back("tdd"); // Thermo
        babelFormats.push_back("therm"); // Thermo
        babelFormats.push_back("tmol"); // TurboMole Coordinate
        babelFormats.push_back("unixyz"); // UniChem XYZ
        babelFormats.push_back("vmol"); // ViewMol
        babelFormats.push_back("xed"); // XED
        babelFormats.push_back("yob"); // YASARA Yob
        babelFormats.push_back("zin"); // Zindo

        // register each specific file format
        foreach(const std::string &format, babelFormats){
            BabelFileFormatInstantiator instantiator(format);

            registerPluginClass<chemkit::MoleculeFileFormat>(format, instantiator);
        }
    }
};

CHEMKIT_EXPORT_PLUGIN(babel, BabelPlugin)
