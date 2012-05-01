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

#include "atomtyper.h"

#include "atom.h"
#include "bond.h"
#include "foreach.h"
#include "molecule.h"
#include "pluginmanager.h"

namespace chemkit {

// === AtomTyperPrivate ==================================================== //
class AtomTyperPrivate
{
public:
    std::string name;
    const Molecule *molecule;
};

// === AtomTyper =========================================================== //
/// \class AtomTyper atomtyper.h chemkit/atomtyper.h
/// \ingroup chemkit
/// \brief The AtomTyper class assigns symbolic types to atoms.
///
/// To create atom typer objects use the AtomTyper::create() method.
///
/// A list of supported atom typers is available at:
/// http://wiki.chemkit.org/Features#Atom_Typers

// --- Construction and Destruction ---------------------------------------- //
AtomTyper::AtomTyper(const std::string &name)
    : d(new AtomTyperPrivate)
{
    d->name = name;
    d->molecule = 0;
}

/// Destroys the atom typer object.
AtomTyper::~AtomTyper()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the name of the atom typer.
std::string AtomTyper::name() const
{
    return d->name;
}

/// Sets the molecule for the atom typer to \p molecule.
void AtomTyper::setMolecule(const Molecule *molecule)
{
    d->molecule = molecule;
}

/// Returns the molecule for the atom typer.
const Molecule* AtomTyper::molecule() const
{
    return d->molecule;
}

// --- Types --------------------------------------------------------------- //
/// Returns the symbolic type for \p atom.
std::string AtomTyper::type(const Atom *atom) const
{
    CHEMKIT_UNUSED(atom);

    return std::string();
}

// --- Interaction Types --------------------------------------------------- //
int AtomTyper::bondedInteractionType(const Atom *a, const Atom *b) const
{
    CHEMKIT_UNUSED(a);
    CHEMKIT_UNUSED(b);

    return 0;
}

int AtomTyper::angleInteractionType(const Atom *a, const Atom *b, const Atom *c) const
{
    CHEMKIT_UNUSED(a);
    CHEMKIT_UNUSED(b);
    CHEMKIT_UNUSED(c);

    return 0;
}

int AtomTyper::torsionInteractionType(const Atom *a, const Atom *b, const Atom *c, const Atom *d) const
{
    CHEMKIT_UNUSED(a);
    CHEMKIT_UNUSED(b);
    CHEMKIT_UNUSED(c);
    CHEMKIT_UNUSED(d);

    return 0;
}

// --- Predicates ---------------------------------------------------------- //
/// Returns \c true if \p atom is a carbon in a carbonyl group.
bool AtomTyper::isCarbonylCarbon(const Atom *atom)
{
    return atom->is(Atom::Carbon) &&
           atom->isBondedTo(Atom::Oxygen, Bond::Double);
}

/// Returns \c true if \p atom is an oxygen in a carbonyl group.
bool AtomTyper::isCarbonylOxygen(const Atom *atom)
{
    return atom->is(Atom::Oxygen) &&
           atom->isTerminal() &&
           atom->isBondedTo(Atom::Carbon, Bond::Double);
}

/// Returns \c true if \p atom is a halogen.
bool AtomTyper::isHalogen(const Atom *atom)
{
    return atom->is(Atom::Fluorine) ||
           atom->is(Atom::Chlorine) ||
           atom->is(Atom::Bromine) ||
           atom->is(Atom::Iodine);
}

/// Returns \c true if \p atom is a hydrogen donor.
bool AtomTyper::isHydrogenDonor(const Atom *atom)
{
    return (atom->is(Atom::Oxygen) ||
            atom->is(Atom::Nitrogen) ||
            atom->is(Atom::Fluorine)) &&
           atom->isBondedTo(Atom::Hydrogen);
}

/// Returns \c true if \p atom is a hydrogen acceptor.
bool AtomTyper::isHydrogenAcceptor(const Atom *atom)
{
    return atom->is(Atom::Oxygen) ||
           atom->is(Atom::Nitrogen) ||
           atom->is(Atom::Fluorine);
}

/// Returns \c true if \p atom is a terminal hydrogen in a hydroxyl
/// group.
bool AtomTyper::isHydroxylHydrogen(const Atom *atom)
{
    return atom->isTerminalHydrogen() &&
           isHydroxylOxygen(atom->neighbor(0));
}

/// Returns \c true if \p atom is an oxygen in a hydroxyl group.
bool AtomTyper::isHydroxylOxygen(const Atom *atom)
{
    return atom->is(Atom::Oxygen) &&
           atom->neighborCount() == 2 &&
           atom->isBondedTo(Atom::Hydrogen, Bond::Single);
}

/// Returns \c true if \p atom is a carbon in a nitrile group.
bool AtomTyper::isNitrileCarbon(const Atom *atom)
{
    return atom->is(Atom::Carbon) &&
           atom->neighborCount() == 2 &&
           atom->isBondedTo(Atom::Nitrogen, Bond::Triple);
}

/// Returns \c true if \p atom is a terminal nitrogen in a nitrile
/// group.
bool AtomTyper::isNitrileNitrogen(const Atom *atom)
{
    return atom->is(Atom::Nitrogen) &&
           atom->isTerminal() &&
           atom->isBondedTo(Atom::Carbon, Bond::Triple);
}

/// Returns \c true if \p atom is an oxygen in a nitro group.
bool AtomTyper::isNitroOxygen(const Atom *atom)
{
    if(!atom->is(Atom::Oxygen) || !atom->isTerminal()){
        return false;
    }

    const Atom *neighbor = atom->neighbor(0);
    if(!neighbor->is(Atom::Nitrogen)){
        return false;
    }

    const Bond *neighborBond = atom->bond(0);
    return (neighborBond->is(Bond::Single) && atom->formalCharge() == -1) ||
           (neighborBond->is(Bond::Double) && atom->formalCharge() == 0);
}

/// Returns \c true if \p atom is a nitrogen in a nitro group.
bool AtomTyper::isNitroNitrogen(const Atom *atom)
{
    return atom->is(Atom::Nitrogen) &&
           atom->neighborCount() == 3 &&
           atom->isBondedTo(Atom::Oxygen, Bond::Single) &&
           atom->isBondedTo(Atom::Oxygen, Bond::Double);
}

/// Returns \c true if \p atom is a terminal hydrogen attached to a
/// polar atom.
bool AtomTyper::isPolarHydrogen(const Atom *atom)
{
    if(!atom->isTerminalHydrogen()){
        return false;
    }

    const Atom *neighbor = atom->neighbor(0);

    return neighbor->is(Atom::Nitrogen) ||
           neighbor->is(Atom::Oxygen) ||
           neighbor->is(Atom::Fluorine);
}

/// Returns \c true if \p atom is a terminal hydrogen attached to a
/// non-polar atom.
bool AtomTyper::isNonpolarHydrogen(const Atom *atom)
{
    return atom->isTerminalHydrogen() && !isPolarHydrogen(atom);
}

/// Returns \c true if \p atom is a hydrogen in a thiol group.
bool AtomTyper::isThiolHydrogen(const Atom *atom)
{
    return atom->isTerminalHydrogen() &&
           isThiolSulfur(atom->neighbor(0));
}

/// Returns \c true if \p atom is a sulfur in a thiol group.
bool AtomTyper::isThiolSulfur(const Atom *atom)
{
    return atom->is(Atom::Sulfur) &&
           atom->neighborCount() == 2 &&
           atom->isBondedTo(Atom::Hydrogen, Bond::Single);
}

// --- Static Methods ------------------------------------------------------ //
/// Creates and returns a new atom typer with \p name. Returns \c 0
/// if \p name is invalid.
AtomTyper* AtomTyper::create(const std::string &name)
{
    return PluginManager::instance()->createPluginClass<AtomTyper>(name);
}

/// Returns a list of names of all the available atom typers.
std::vector<std::string> AtomTyper::typers()
{
    return PluginManager::instance()->pluginClassNames<AtomTyper>();
}

/// This static convenience function assigns atom types for atoms
/// in molecule using the specified \p typer.
///
/// For example, to assign SYBYL atom types for a molecule:
/// \code
/// AtomTyper::assignAtomTypes(molecule, "sybyl");
/// \endcode
bool AtomTyper::assignAtomTypes(Molecule *molecule, const std::string &typer)
{
    boost::scoped_ptr<AtomTyper> atomTyper(create(typer));
    if(!atomTyper){
        return false;
    }

    atomTyper->setMolecule(molecule);

    foreach(Atom *atom, molecule->atoms()){
        atom->setType(atomTyper->type(atom));
    }

    return true;
}

} // end chemkit namespace
