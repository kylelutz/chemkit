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

#include "molecule.h"

#include <map>
#include <queue>
#include <sstream>
#include <algorithm>

#include <boost/scoped_ptr.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "atom.h"
#include "bond.h"
#include "ring.h"
#include "point3.h"
#include "element.h"
#include "foreach.h"
#include "vector3.h"
#include "fragment.h"
#include "geometry.h"
#include "conformer.h"
#include "constants.h"
#include "lineformat.h"
#include "quaternion.h"
#include "coordinates.h"
#include "moleculargraph.h"
#include "moleculeprivate.h"
#include "moleculeobserver.h"
#include "internalcoordinates.h"
#include "moleculardescriptor.h"

namespace chemkit {

// === MoleculePrivate ===================================================== //
MoleculePrivate::MoleculePrivate()
{
    conformer = 0;
    fragmentsPerceived = false;
    ringsPerceived = false;
}

// === Molecule ============================================================ //
/// \class Molecule molecule.h chemkit/molecule.h
/// \ingroup chemkit
/// \brief The Molecule class represents a chemical molecule.
///
/// The diagram below shows the various components contained in a
/// molecule object. The molecule object shown contains 29 Atom%'s,
/// 27 Bond%'s, 3 Ring%'s, and 2 Fragment%'s. The molecule image was
/// created using the GraphicsMoleculeItem class.
/// \image html molecule-labels.png
///
/// Molecules can be created in two different ways. The examples below
/// show two methods for creating a new water molecule:
///
/// 1. By adding every atom and bond explicitly:
/// \code
/// Molecule *molecule = new Molecule();
/// Atom *O1 = molecule->addAtom("O");
/// Atom *H2 = molecule->addAtom("H");
/// Atom *H3 = molecule->addAtom("H");
/// molecule->addBond(O1, H2);
/// molecule->addBond(O1, H3);
/// \endcode
///
/// 2. From a chemical line format formula such as InChI or SMILES:
/// \code
/// Molecule *molecule = new Molecule("InChI=1/H2O/h1H2", "inchi");
/// \endcode
///
/// Molecules can also be read from existing files using the
/// MoleculeFile class.
///
/// Molecule objects take ownership of all the Atom, Bond, Ring,
/// Fragment, and Conformer objects that they contain. Deleting the
/// molecule will also delete all of the objects that it contains.

/// \enum Molecule::CompareFlag
/// Option flags for molecule comparisons.
///    - \c CompareAtomsOnly
///    - \c CompareHydrogens
///    - \c CompareAromaticity

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new, empty molecule.
Molecule::Molecule()
    : d(new MoleculePrivate)
{
    m_stereochemistry = 0;
}

/// Creates a new molecule from its formula.
///
/// The following code creates a new benzene molecule from its InChI
/// formula:
/// \code
/// Molecule *benzene = new Molecule("InChI=1/C6H6/c1-2-4-6-5-3-1/h1-6H", "inchi");
/// \endcode
///
/// \see LineFormat
Molecule::Molecule(const std::string &formula, const std::string &format)
    : d(new MoleculePrivate)
{
    m_stereochemistry = 0;

    boost::scoped_ptr<LineFormat> lineFormat(LineFormat::create(format));
    if(!lineFormat){
        return;
    }

    lineFormat->read(formula, this);
}

/// Creates a new molecule that is a copy of \p molecule.
Molecule::Molecule(const Molecule &molecule)
    : d(new MoleculePrivate)
{
    m_stereochemistry = 0;

    d->name = molecule.name();

    std::map<const Atom *, Atom *> oldToNew;

    foreach(const Atom *atom, molecule.atoms()){
        Atom *newAtom = addAtomCopy(atom);
        oldToNew[atom] = newAtom;
    }

    foreach(const Bond *bond, molecule.bonds()){
        Bond *newBond = addBond(oldToNew[bond->atom1()], oldToNew[bond->atom2()]);
        newBond->setOrder(bond->order());
    }
}

/// Destroys a molecule. This also destroys all of the atoms and
/// bonds that the molecule contains.
Molecule::~Molecule()
{
    foreach(Atom *atom, m_atoms)
        delete atom;
    foreach(Bond *bond, d->bonds)
        delete bond;
    foreach(Ring *ring, d->rings)
        delete ring;
    foreach(Fragment *fragment, d->fragments)
        delete fragment;
    foreach(Conformer *conformer, d->conformers)
        delete conformer;

    delete m_stereochemistry;

    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the name of the molecule.
void Molecule::setName(const std::string &name)
{
    d->name = name;
    notifyObservers(NameChanged);
}

/// Returns the name of the molecule.
std::string Molecule::name() const
{
    return d->name;
}

/// Returns the chemical formula (e.g. "H2O") for the molecule. The
/// formula is formated according to the Hill system.
std::string Molecule::formula() const
{
    // a map of atomic symbols to their quantity
    std::map<std::string, int> composition;
    foreach(const Atom *atom, m_atoms){
        composition[atom->symbol()]++;
    }

    std::stringstream formula;

    if(composition.count("C") != 0){
        formula << "C";
        if(composition["C"] > 1){
            formula << composition["C"];
        }
        composition.erase("C");

        if(composition.count("H") != 0){
            formula << "H";
            if(composition["H"] > 1){
                formula << composition["H"];
            }
        }
        composition.erase("H");
    }

    std::map<std::string, int>::iterator iter;
    for(iter = composition.begin(); iter != composition.end(); ++iter){
        formula << iter->first;

        if(iter->second > 1){
            formula << iter->second;
        }
    }

    return formula.str();
}

/// Returns the the formula of the molecule using the specified
/// format. Returns an empty string if format is not supported or if
/// an error occurs.
///
/// The following example returns the InChI formula for a molecule:
/// \code
/// molecule->formula("inchi");
/// \endcode
///
/// \see LineFormat
std::string Molecule::formula(const std::string &format) const
{
    boost::scoped_ptr<LineFormat> lineFormat(LineFormat::create(format));
    if(!lineFormat){
        return std::string();
    }

    return lineFormat->write(this);
}

/// Calculates and returns the molecular descriptor \p name. If the
/// descriptor is not available or the calculation fails a null
/// Variant is returned.
///
/// For example, to calculate the Randic index of the molecule use:
/// \code
/// double randicIndex = molecule->descriptor("randic-index").toDouble();
/// \endcode
///
/// \see MolecularDescriptor
Variant Molecule::descriptor(const std::string &name) const
{
    boost::scoped_ptr<MolecularDescriptor> descriptor(MolecularDescriptor::create(name));
    if(!descriptor){
        return Variant();
    }

    return descriptor->value(this);
}

/// Returns the total molar mass of the molecule. Mass is in g/mol.
Real Molecule::mass() const
{
    Real mass = 0;

    foreach(const Atom *atom, m_atoms)
        mass += atom->mass();

    return mass;
}

/// Sets the data for the molecule with \p name to \p value.
void Molecule::setData(const std::string &name, const Variant &value)
{
    d->data[name] = value;
}

/// Returns the data for the molecule with \p name.
Variant Molecule::data(const std::string &name) const
{
    std::map<std::string, Variant>::const_iterator iter = d->data.find(name);
    if(iter != d->data.end()){
        return iter->second;
    }

    return Variant();
}

// --- Structure ----------------------------------------------------------- //
/// Adds a new atom of the given \p element to the molecule.
///
/// The Element class has a number of constructors which take
/// either an atomic number or an element symbol. Any of the
/// following lines of code will add a new Carbon atom to the
/// molecule:
///
/// \code
/// // add atom from its symbol
/// molecule->addAtom("C");
/// \endcode
///
/// \code
/// // add atom from its atomic number
/// molecule->addAtom(6);
/// \endcode
///
/// \code
/// // add atom from the atom enum name
/// molecule->addAtom(Atom::Carbon);
/// \endcode
Atom* Molecule::addAtom(const Element &element)
{
    Atom *atom = new Atom(this, m_atoms.size());
    m_atoms.push_back(atom);

    // add atom properties
    m_elements.push_back(element);
    d->atomBonds.push_back(std::vector<Bond *>());
    d->positions.push_back(Point3(0, 0, 0));
    d->partialCharges.push_back(0);

    setFragmentsPerceived(false);
    notifyObservers(atom, AtomAdded);

    return atom;
}

/// Adds a new atom to the molecule. The new atom will have the same
/// properties as atom (atomic number, mass number, etc).
Atom* Molecule::addAtomCopy(const Atom *atom)
{
    Atom *newAtom = addAtom(atom->element());

    newAtom->setMassNumber(atom->massNumber());
    newAtom->setPartialCharge(atom->partialCharge());
    newAtom->setPosition(atom->position());

    if(atom->chirality() != Stereochemistry::None){
        newAtom->setChirality(atom->chirality());
    }

    return newAtom;
}

/// Removes atom from the molecule. This will also remove any bonds
/// to/from the atom.
void Molecule::removeAtom(Atom *atom)
{
    if(!contains(atom)){
        return;
    }

    // remove all bonds to/from the atom first
    std::vector<Bond *> bonds = atom->bonds();
    foreach(Bond *bond, bonds){
        removeBond(bond);
    }

    m_atoms.erase(std::remove(m_atoms.begin(), m_atoms.end(), atom), m_atoms.end());

    // remove atom properties
    m_elements.erase(m_elements.begin() + atom->index());
    d->isotopes.erase(atom);
    d->positions.erase(d->positions.begin() + atom->index());
    d->atomBonds.erase(d->atomBonds.begin() + atom->index());
    d->partialCharges.erase(d->partialCharges.begin() + atom->index());

    // subtract one from the index of all atoms after this one
    for(unsigned int i = atom->m_index; i < m_atoms.size(); i++){
        m_atoms[i]->m_index--;
    }

    atom->m_molecule = 0;
    notifyObservers(atom, AtomRemoved);

    delete atom;
}

/// Returns the number of atoms in the molecule of the given
/// \p element.
int Molecule::atomCount(const Element &element) const
{
    int count = 0;

    foreach(const Atom *atom, m_atoms){
        if(atom->is(element)){
            count++;
        }
    }

    return count;
}

/// Returns \c true if the molecule contains atom.
bool Molecule::contains(const Atom *atom) const
{
    return atom->molecule() == this;
}

/// Returns \c true if the molecule contains an atom of the given
/// \p element.
bool Molecule::contains(const Element &element) const
{
    foreach(const Atom *atom, m_atoms){
        if(atom->is(element)){
            return true;
        }
    }

    return false;
}

/// Adds a new bond between atoms \p a and \p b and returns it. If
/// they are already bonded the existing bond is returned.
Bond* Molecule::addBond(Atom *a, Atom *b, int order)
{
    // ensure that the atoms are not the same
    if(a == b){
        return 0;
    }

    // ensure that this molecule contains both atoms
    if(!contains(a) || !contains(b)){
        return 0;
    }

    // check to see if they are already bonded
    if(a->isBondedTo(b)){
        return bond(a, b);
    }

    Bond *bond = new Bond(this, d->bonds.size());
    d->atomBonds[a->index()].push_back(bond);
    d->atomBonds[b->index()].push_back(bond);
    d->bonds.push_back(bond);

    // add bond properties
    d->bondAtoms.push_back(std::make_pair(a, b));
    d->bondOrders.push_back(order);

    setRingsPerceived(false);
    setFragmentsPerceived(false);

    notifyObservers(bond, BondAdded);

    return bond;
}

/// Adds a new bond between atoms with indicies \p a and \p b.
Bond* Molecule::addBond(int a, int b, int order)
{
    return addBond(atom(a), atom(b), order);
}

/// Removes \p bond from the molecule.
void Molecule::removeBond(Bond *bond)
{
    assert(bond->molecule() == this);

    d->bonds.erase(d->bonds.begin() + bond->index());

    // remove bond from atom bond vectors
    std::vector<Bond *> &bondsA = d->atomBonds[bond->atom1()->index()];
    std::vector<Bond *> &bondsB = d->atomBonds[bond->atom2()->index()];

    bondsA.erase(std::find(bondsA.begin(), bondsA.end(), bond));
    bondsB.erase(std::find(bondsB.begin(), bondsB.end(), bond));

    // remove bond properties
    d->bondAtoms.erase(d->bondAtoms.begin() + bond->index());
    d->bondOrders.erase(d->bondOrders.begin() + bond->index());

    // subtract one from the index of all bonds after this one
    for(unsigned int i = bond->index(); i < d->bonds.size(); i++){
        d->bonds[i]->m_index--;
    }

    setRingsPerceived(false);
    setFragmentsPerceived(false);

    notifyObservers(bond, BondRemoved);

    delete bond;
}

/// Removes the bond between atoms \p a and \p b. Does nothing if
/// they are not bonded.
void Molecule::removeBond(Atom *a, Atom *b)
{
    Bond *bond = this->bond(a, b);

    if(bond){
        removeBond(bond);
    }
}

/// Removes the bond between atoms with indicies \p a and \p b.
void Molecule::removeBond(int a, int b)
{
    removeBond(bond(a, b));
}

/// Returns a list of all the bonds in the molecule.
std::vector<Bond *> Molecule::bonds() const
{
    return d->bonds;
}

/// Returns the number of bonds in the molecule.
int Molecule::bondCount() const
{
    return bonds().size();
}

/// Returns the bond at index.
Bond* Molecule::bond(int index) const
{
    return d->bonds[index];
}

/// Returns the bond between atom \p a and \p b. Returns \c 0 if they
/// are not bonded.
///
/// To create a new bond between the atoms use Molecule::addBond().
Bond* Molecule::bond(const Atom *a, const Atom *b) const
{
    return const_cast<Atom *>(a)->bondTo(b);
}

/// Returns the bond between the atoms with indicies \p a and \p b.
Bond* Molecule::bond(int a, int b) const
{
    return bond(atom(a), atom(b));
}

/// Returns \c true if the molecule contains bond.
bool Molecule::contains(const Bond *bond) const
{
    return contains(bond->atom1());
}

/// Removes all atoms and bonds from the molecule.
void Molecule::clear()
{
    std::vector<Bond *> bonds = d->bonds;
    foreach(Bond *bond, bonds){
        removeBond(bond);
    }

    std::vector<Atom *> atoms = m_atoms;
    foreach(Atom *atom, atoms){
        removeAtom(atom);
    }
}

// --- Comparison ---------------------------------------------------------- //
/// Returns \c true if the molecule equals \p molecule.
bool Molecule::equals(const Molecule *molecule, int flags) const
{
    return contains(molecule, flags) && molecule->contains(this, flags);
}

/// Returns \c true if the molecule contains \p molecule as a
/// substructure.
///
/// For example, this method could be used to create a function that
/// checks if a molecule contains a carboxyl group (-COO):
/// \code
/// bool containsCarboxylGroup(const Molecule *molecule)
/// {
///      Molecule carboxyl;
///      Atom *C1 = carboxyl.addAtom("C");
///      Atom *O2 = carboxyl.addAtom("O");
///      Atom *O3 = carboxyl.addAtom("O");
///      carboxyl.addBond(C1, O2, Bond::Double);
///      carboxyl.addBond(C1, O3, Bond::Single);
///
///      return molecule->contains(&carboxyl);
/// }
/// \endcode
bool Molecule::contains(const Molecule *molecule, int flags) const
{
    if(molecule == this){
        return true;
    }

    if(isEmpty() && molecule->isEmpty()){
        return true;
    }
    else if((flags & CompareAtomsOnly) || (bondCount() == 0 && molecule->bondCount() == 0)){
        return molecule->isSubsetOf(this, flags);
    }

    return !molecule->mapping(this, flags).empty();
}

/// Returns \c true if the molecule is a substructure of \p molecule.
bool Molecule::isSubstructureOf(const Molecule *molecule, int flags) const
{
    return molecule->contains(this, flags);
}

/// Returns a mapping (also known as an isomorphism) between the atoms
/// in the molecule and the atoms in \p molecule.
std::map<Atom *, Atom *> Molecule::mapping(const Molecule *molecule, int flags) const
{
    MolecularGraph *source;
    MolecularGraph *target;

    if(flags & CompareHydrogens){
        source = new MolecularGraph(this);
        target = new MolecularGraph(molecule);
    }
    else{
        source = MolecularGraph::hydrogenDepletedGraph(this);
        target = MolecularGraph::hydrogenDepletedGraph(molecule);
    }

    // label for aroamtic bonds
    const int aromaticBondLabel = 10;

    if(flags & CompareAromaticity){
        for(unsigned int i = 0; i < source->bondCount(); i++){
            const Bond *bond = source->bond(i);

            if(bond->isAromatic()){
                source->setBondLabel(i, aromaticBondLabel);
            }
        }

        for(unsigned int i = 0; i < target->bondCount(); i++){
            const Bond *bond = target->bond(i);

            if(bond->isAromatic()){
                target->setBondLabel(i, aromaticBondLabel);
            }
        }
    }

    std::map<Atom *, Atom *> mapping = MolecularGraph::isomorphism(source, target);

    delete source;
    delete target;

    return mapping;
}

/// Searches the molecule for an occurrence of \p moiety and returns
/// it if found. If not found an empty moiety is returned.
///
/// For example, to find an amide group (NC=O) in the molecule:
/// \code
/// Molecule amide;
/// Atom *C1 = amide.addAtom("C");
/// Atom *N2 = amide.addAtom("N");
/// Atom *O3 = amide.addAtom("O");
/// amide.addBond(C1, N2, Bond::Single);
/// amide.addBond(C1, O3, Bond::Double);
///
/// Moiety amideGroup = molecule.find(&amide);
/// \endcode
Moiety Molecule::find(const Molecule *moiety, int flags) const
{
    std::map<Atom *, Atom *> mapping = moiety->mapping(this, flags);

    // no mapping found, return empty moiety
    if(mapping.empty()){
        return Moiety();
    }

    std::vector<Atom *> moietyAtoms;
    foreach(Atom *atom, moiety->atoms()){
        moietyAtoms.push_back(const_cast<Atom *>(mapping[atom]));
    }

    return Moiety(moietyAtoms);
}

// --- Ring Perception ----------------------------------------------------- //
/// Returns the ring at \p index.
///
/// Equivalent to calling:
/// \code
/// molecule.rings()[index];
/// \endcode
Ring* Molecule::ring(int index) const
{
    return rings()[index];
}

/// Returns a list of all rings in the molecule.
///
/// \warning The list of rings returned from this method is only
///          valid as long as the molecule's structure remains
///          unchanged. If any atoms or bonds in the molecule are
///          added or removed the old results must be discarded and
///          this method must be called again.
std::vector<Ring *> Molecule::rings() const
{
    // only run ring perception if neccessary
    if(!ringsPerceived()){
        // find rings
        d->rings = MolecularGraph::sssr(this);

        // set perceived to true
        setRingsPerceived(true);
    }

    return d->rings;
}

/// Returns the number of rings in the molecule.
int Molecule::ringCount() const
{
    return rings().size();
}

void Molecule::setRingsPerceived(bool perceived) const
{
    if(perceived == d->ringsPerceived){
        return;
    }

    if(perceived == false){
        foreach(Ring *ring, d->rings){
            delete ring;
        }

        d->rings.clear();
    }

    d->ringsPerceived = perceived;
}

bool Molecule::ringsPerceived() const
{
    return d->ringsPerceived;
}

// --- Fragment Perception-------------------------------------------------- //
/// Returns the fragment at \p index.
///
/// Equivalent to calling:
/// \code
/// molecule.fragments()[index];
/// \endcode
Fragment* Molecule::fragment(int index) const
{
    return fragments()[index];
}

/// Returns a list of fragments in the molecule.
///
/// \warning The list of fragments returned from this method is only
///          valid as long as the molecule's structure remains
///          unchanged. If any atoms or bonds in the molecule are
///          added or removed the old results must be discarded and
///          this method must be called again.
std::vector<Fragment *> Molecule::fragments() const
{
    if(!fragmentsPerceived()){
        perceiveFragments();

        setFragmentsPerceived(true);
    }

    return d->fragments;
}

/// Returns the number of fragments in the molecule.
int Molecule::fragmentCount() const
{
    return fragments().size();
}

/// Returns \c true if the molecule is fragmented. (i.e. contains
/// more than one fragment).
bool Molecule::isFragmented() const
{
    return fragmentCount() > 1;
}

/// Removes all of the atoms and bonds contained in \p fragment from
/// the molecule.
void Molecule::removeFragment(Fragment *fragment)
{
    foreach(Atom *atom, fragment->atoms()){
        removeAtom(atom);
    }
}

Fragment* Molecule::fragment(const Atom *atom) const
{
    foreach(Fragment *fragment, fragments()){
        if(fragment->contains(atom)){
            return fragment;
        }
    }

    return 0;
}

void Molecule::setFragmentsPerceived(bool perceived) const
{
    if(perceived == d->fragmentsPerceived)
        return;

    if(!perceived){
        foreach(const Fragment *fragment, d->fragments){
            delete fragment;
        }

        d->fragments.clear();
    }

    d->fragmentsPerceived = perceived;
}

bool Molecule::fragmentsPerceived() const
{
    return d->fragmentsPerceived;
}

void Molecule::perceiveFragments() const
{
    if(isEmpty()){
        // nothing to do
        return;
    }

    // position of the next atom to root the depth-first search
    size_t position = 0;

    // bitset marking each atom not visited yet
    boost::dynamic_bitset<> unvisited(m_atoms.size());
    unvisited.set();

    for(;;){
        // bitset marking the atoms contained in the fragment
        boost::dynamic_bitset<> bitset(m_atoms.size());

        // perform depth-first search
        std::vector<const Atom *> row;
        row.push_back(m_atoms[position]);

        while(!row.empty()){
            std::vector<const Atom *> nextRow;

            foreach(const Atom *atom, row){
                bitset.set(atom->index());
                unvisited.set(atom->index(), false);

                foreach(const Atom *neighbor, atom->neighbors()){
                    if(unvisited[neighbor->index()]){
                        nextRow.push_back(neighbor);
                    }
                }
            }

            row = nextRow;
        }

        // create and add fragment
        Fragment *fragment = new Fragment(const_cast<Molecule *>(this), bitset);
        d->fragments.push_back(fragment);

        // find next unvisited atom
        position = unvisited.find_next(position);
        if(position == boost::dynamic_bitset<>::npos){
            break;
        }
    }
}

// --- Geometry ------------------------------------------------------------ //
/// Sets the coordinates for the atoms in the molecule to
/// \p coordinates.
void Molecule::setCoordinates(const Coordinates *coordinates)
{
    int size = std::min(this->size(), coordinates->size());

    for(int i = 0; i < size; i++){
        m_atoms[i]->setPosition(coordinates->position(i));
    }
}

/// Sets the coordinates for the atoms in the molecule to
/// \p coordinates.
void Molecule::setCoordinates(const InternalCoordinates *coordinates)
{
    Coordinates *cartesianCoordinates = coordinates->toCartesianCoordinates();
    setCoordinates(cartesianCoordinates);
    delete cartesianCoordinates;
}

/// Returns the distance between atoms \p a and \p b. The returned
/// distance is in Angstroms.
Real Molecule::distance(const Atom *a, const Atom *b) const
{
    return chemkit::geometry::distance(a->position(), b->position());
}

/// Returns the angle between atoms \p a, \p b, and \p c. The
/// returned angle is in degrees.
Real Molecule::bondAngle(const Atom *a, const Atom *b, const Atom *c) const
{
    return chemkit::geometry::angle(a->position(), b->position(), c->position());
}

/// Returns the torsion angle (also known as the dihedral angle)
/// between atoms \p a, \p b, \p c, and \p d. The returned angle is
/// in degrees.
Real Molecule::torsionAngle(const Atom *a, const Atom *b, const Atom *c, const Atom *d) const
{
    return chemkit::geometry::torsionAngle(a->position(), b->position(), c->position(), d->position());
}

/// Returns the wilson angle between the plane made by atoms \p a,
/// \p b, \p c and the vector from \p c to \p d. The returned angle
/// is in degrees.
Real Molecule::wilsonAngle(const Atom *a, const Atom *b, const Atom *c, const Atom *d) const
{
    return chemkit::geometry::wilsonAngle(a->position(), b->position(), c->position(), d->position());
}

/// Moves all of the atoms in the molecule so that the center point
/// is at \p position.
void Molecule::setCenter(const Point3 &position)
{
    moveBy(position - center());
}

/// Moves all of the atoms in the molecule so that the new center
/// point is at (\p x, \p y, \p z). This convenience function is
/// equivalent to calling setCenter(Point(\p x, \p y, \p z)).
void Molecule::setCenter(Real x, Real y, Real z)
{
    setCenter(Point3(x, y, z));
}

/// Returns the center point of the molecule. This is also known as
/// the centriod.
///
/// \see centerOfMass()
Point3 Molecule::center() const
{
    if(isEmpty()){
        return Point3(0, 0, 0);
    }

    // sums for each component
    Real sx = 0;
    Real sy = 0;
    Real sz = 0;

    foreach(const Atom *atom, m_atoms){
        sx += atom->x();
        sy += atom->y();
        sz += atom->z();
    }

    int n = atomCount();

    return Point3(sx/n, sy/n, sz/n);
}

/// Returns the center of mass for the molecule.
Point3 Molecule::centerOfMass() const
{
    if(isEmpty()){
        return Point3();
    }

    // sums for each component
    Real sx = 0;
    Real sy = 0;
    Real sz = 0;

    // sum of weights
    Real sw = 0;

    foreach(const Atom *atom, m_atoms){
        Real w = atom->mass();

        sx += w * atom->x();
        sy += w * atom->y();
        sz += w * atom->z();

        sw += w;
    }

    int n = sw * size();

    return Point3(sx/n, sy/n, sz/n);
}

/// Moves all the atoms in the molecule by \p vector.
void Molecule::moveBy(const Vector3 &vector)
{
    foreach(Atom *atom, m_atoms){
        atom->moveBy(vector);
    }
}

/// Moves all of the atoms in the molecule by (\p dx, \p dy, \p dz).
void Molecule::moveBy(Real dx, Real dy, Real dz)
{
    foreach(Atom *atom, m_atoms){
        atom->moveBy(dx, dy, dz);
    }
}

/// Rotates the positions of all the atoms in the molecule
/// by \p angle degrees around \p axis.
void Molecule::rotate(const Vector3 &axis, Real angle)
{
    // convert angle to radians
    angle *= chemkit::constants::DegreesToRadians;

    // build rotation transform
    Eigen::Matrix<Real, 3, 1> axisVector(axis.x(), axis.y(), axis.z());
    Eigen::Transform<Real, 3, 3> transform(Eigen::AngleAxis<Real>(angle, axisVector));

    // rotate each atom
    foreach(Atom *atom, m_atoms){
        Eigen::Matrix<Real, 3, 1> position(atom->x(), atom->y(), atom->z());

        position = transform * position;

        atom->setPosition(position.x(), position.y(), position.z());
    }
}

/// Returns \c true if the molecule has coordinates for any of the
/// atoms.
bool Molecule::hasCoordinates() const
{
    foreach(const Atom *atom, m_atoms){
        if(!atom->position().isZero()){
            return true;
        }
    }

    return false;
}

/// Removes all of the atomic coordinates in the molecule.
void Molecule::clearCoordinates()
{
    foreach(Atom *atom, m_atoms){
        atom->setPosition(Point3());
    }
}

// --- Conformers ---------------------------------------------------------- //
/// Adds a new conformer to the molecule and returns it.
Conformer* Molecule::addConformer()
{
    if(!d->conformer)
        conformers();

    Conformer *conformer = new Conformer(this);
    d->conformers.push_back(conformer);
    return conformer;
}

/// Removes \p conformer from the molecule. The currently active
/// conformer cannot be removed.
void Molecule::removeConformer(Conformer *conformer)
{
    // forbid removal of the currently active conformer
    if(conformer == d->conformer){
        return;
    }

    std::vector<Conformer *>::iterator location = std::find(d->conformers.begin(), d->conformers.end(), conformer);
    if(location == d->conformers.end()){
        return;
    }

    d->conformers.erase(location);

    delete conformer;
}

/// Sets the active conformer for the molecule.
void Molecule::setConformer(Conformer *conformer)
{
    if(conformer == d->conformer){
        return;
    }
    else if(conformer->molecule() != this){
        return;
    }

    foreach(Atom *atom, m_atoms){
        atom->setPosition(conformer->position(atom));
    }

    d->conformer = conformer;
}

/// Returns the active conformer for the molecule.
Conformer* Molecule::conformer() const
{
    if(!d->conformer){
        d->conformer = conformers()[0];
    }

    return d->conformer;
}

/// Returns the conformer at \p index.
///
/// Equivalent to:
/// \code
/// molelcule.conformers()[index];
/// \endcode
Conformer* Molecule::conformer(int index) const
{
    return conformers()[index];
}

/// Returns a list of all conformers in the molecule.
std::vector<Conformer *> Molecule::conformers() const
{
    if(d->conformers.empty()){
        d->conformer = new Conformer(this);
        d->conformers.push_back(d->conformer);
    }

    return d->conformers;
}

/// Returns the number of conformers in the molecule.
int Molecule::conformerCount() const
{
    return conformers().size();
}

// --- Operators ----------------------------------------------------------- //
Molecule& Molecule::operator=(const Molecule &molecule)
{
    if(this != &molecule){
        // clear current molecule
        clear();

        // set new name
        setName(molecule.name());

        std::map<const Atom *, Atom *> oldToNew;

        // add new atoms
        foreach(const Atom *atom, molecule.atoms()){
            Atom *newAtom = addAtomCopy(atom);
            oldToNew[atom] = newAtom;
        }

        // add new bonds
        foreach(const Bond *bond, molecule.bonds()){
            Bond *newBond = addBond(oldToNew[bond->atom1()], oldToNew[bond->atom2()]);
            newBond->setOrder(bond->order());
        }
    }

    return *this;
}

// --- Internal Methods ---------------------------------------------------- //
void Molecule::notifyObservers(ChangeType type)
{
    foreach(MoleculeObserver *observer, d->observers){
        observer->moleculeChanged(this, type);
    }
}

void Molecule::notifyObservers(const Atom *atom, ChangeType type)
{
    foreach(MoleculeObserver *observer, d->observers){
        observer->atomChanged(atom, type);
    }
}

void Molecule::notifyObservers(const Bond *bond, ChangeType type)
{
    foreach(MoleculeObserver *observer, d->observers){
        observer->bondChanged(bond, type);
    }
}

void Molecule::addObserver(MoleculeObserver *observer) const
{
    d->observers.push_back(observer);
}

void Molecule::removeObserver(MoleculeObserver *observer) const
{
    d->observers.erase(std::remove(d->observers.begin(), d->observers.end(), observer));
}

bool Molecule::isSubsetOf(const Molecule *molecule, int flags) const
{
    CHEMKIT_UNUSED(flags);

    std::vector<Atom *> otherAtoms = molecule->atoms();

    foreach(const Atom *atom, m_atoms){
        bool found = false;

        foreach(Atom *otherAtom, otherAtoms){
            if(atom->atomicNumber() == otherAtom->atomicNumber()){
                otherAtoms.erase(std::remove(otherAtoms.begin(), otherAtoms.end(), otherAtom), otherAtoms.end());
                found = true;
                break;
            }
        }

        if(!found){
            return false;
        }
    }

    return true;
}

Stereochemistry* Molecule::stereochemistry()
{
    if(!m_stereochemistry){
        m_stereochemistry = new Stereochemistry(this);
    }

    return m_stereochemistry;
}

} // end chemkit namespace
