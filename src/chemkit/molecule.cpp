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
#include <boost/make_shared.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "vf2.h"
#include "atom.h"
#include "bond.h"
#include "ring.h"
#include "graph.h"
#include "bitset.h"
#include "point3.h"
#include "rppath.h"
#include "element.h"
#include "foreach.h"
#include "vector3.h"
#include "fragment.h"
#include "geometry.h"
#include "constants.h"
#include "lineformat.h"
#include "quaternion.h"
#include "variantmap.h"
#include "fingerprint.h"
#include "moleculeprivate.h"
#include "moleculewatcher.h"
#include "diagramcoordinates.h"
#include "internalcoordinates.h"
#include "moleculardescriptor.h"
#include "cartesiancoordinates.h"

namespace chemkit {

// === MoleculePrivate ===================================================== //
MoleculePrivate::MoleculePrivate()
{
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
/// Fragment, and CoordinateSet objects that they contain. Deleting
/// the molecule will also delete all of the objects that it contains.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new, empty molecule.
Molecule::Molecule()
    : d(new MoleculePrivate)
{
    m_coordinates = 0;
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
/// A list of supported formats is available at:
/// http://wiki.chemkit.org/Features#Line_Formats
///
/// \see LineFormat
Molecule::Molecule(const std::string &formula, const std::string &format)
    : d(new MoleculePrivate)
{
    m_coordinates = 0;
    m_stereochemistry = 0;

    boost::scoped_ptr<LineFormat> lineFormat(LineFormat::create(format));
    if(!lineFormat){
        return;
    }

    boost::scoped_ptr<Molecule> molecule(lineFormat->read(formula));
    if(!molecule){
        return;
    }

    *this = *molecule;
}

/// Creates a new molecule that is a copy of \p molecule.
Molecule::Molecule(const Molecule &molecule)
    : d(new MoleculePrivate)
{
    m_coordinates = 0;
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

        if(bond->stereochemistry() != Stereochemistry::None){
            newBond->setStereochemistry(bond->stereochemistry());
        }
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

    // delete coordinates and all coordinate sets
    bool deletedCoordinates = false;

    foreach(const boost::shared_ptr<CoordinateSet> &coordinateSet, d->coordinateSets){
        if(coordinateSet->type() == CoordinateSet::Cartesian &&
           coordinateSet->cartesianCoordinates() == m_coordinates){
            deletedCoordinates = true;
        }
    }

    if(!deletedCoordinates){
        delete m_coordinates;
    }

    delete m_stereochemistry;

    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the name of the molecule.
void Molecule::setName(const std::string &name)
{
    d->name = name;
    notifyWatchers(MoleculeWatcher::NameChanged);
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
    std::map<std::string, size_t> composition;
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

    std::map<std::string, size_t>::iterator iter;
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
/// A list of supported formats is available at:
/// http://wiki.chemkit.org/Features#Line_Formats
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
/// A list of supported molecular descriptors is available at:
/// http://wiki.chemkit.org/Features#Molecular_Descriptors
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

/// Returns the binary fingerprint for \p name.
///
/// A list of supported fingerprints is available at:
/// http://wiki.chemkit.org/Features#Fingerprints
///
/// \see Fingerprint
Bitset Molecule::fingerprint(const std::string &name) const
{
    boost::scoped_ptr<Fingerprint> fingerprint(Fingerprint::create(name));
    if(!fingerprint){
        return Bitset();
    }

    return fingerprint->value(this);
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
    VariantMap::const_iterator iter = d->data.find(name);
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
    d->partialCharges.push_back(0);

    // set atom position
    if(m_coordinates){
        m_coordinates->append(0, 0, 0);
    }

    setFragmentsPerceived(false);
    notifyWatchers(atom, MoleculeWatcher::AtomAdded);

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
    removeBonds(atom->bonds());

    m_atoms.erase(std::remove(m_atoms.begin(), m_atoms.end(), atom), m_atoms.end());

    // remove atom properties
    m_elements.erase(m_elements.begin() + atom->index());
    d->isotopes.erase(atom);
    d->atomBonds.erase(d->atomBonds.begin() + atom->index());
    d->partialCharges.erase(d->partialCharges.begin() + atom->index());

    if(atom->index() < d->atomTypes.size()){
        d->atomTypes.erase(d->atomTypes.begin() + atom->index());
    }

    if(m_coordinates){
        m_coordinates->remove(atom->index());
    }

    // subtract one from the index of all atoms after this one
    for(size_t i = atom->m_index; i < m_atoms.size(); i++){
        m_atoms[i]->m_index--;
    }

    atom->m_molecule = 0;
    setFragmentsPerceived(false);
    notifyWatchers(atom, MoleculeWatcher::AtomRemoved);

    delete atom;
}

/// Removes each atom in \p atoms from the molecule.
void Molecule::removeAtoms(const std::vector<Atom *> &atoms)
{
    BOOST_REVERSE_FOREACH(Atom *atom, atoms){
        removeAtom(atom);
    }
}

/// Returns the number of atoms in the molecule of the given
/// \p element.
size_t Molecule::atomCount(const Element &element) const
{
    return std::count(m_elements.begin(), m_elements.end(), element);
}

/// Requests that the atom capacity for the molecule be increased to
/// \p capacity.
///
/// \internal
void Molecule::setAtomCapacity(size_t capacity)
{
    m_atoms.reserve(capacity);
    m_elements.reserve(capacity);
    d->atomBonds.reserve(capacity);
    d->partialCharges.reserve(capacity);
}

/// Returns the atom capacity for the molecule.
///
/// \internal
size_t Molecule::atomCapacity() const
{
    return m_atoms.capacity();
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
    return std::find(m_elements.begin(), m_elements.end(), element) != m_elements.end();
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

    notifyWatchers(bond, MoleculeWatcher::BondAdded);

    return bond;
}

/// Adds a new bond between atoms with indices \p a and \p b.
Bond* Molecule::addBond(size_t a, size_t b, int order)
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
    for(size_t i = bond->index(); i < d->bonds.size(); i++){
        d->bonds[i]->m_index--;
    }

    setRingsPerceived(false);
    setFragmentsPerceived(false);

    notifyWatchers(bond, MoleculeWatcher::BondRemoved);

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

/// Removes the bond between atoms with indices \p a and \p b.
void Molecule::removeBond(size_t a, size_t b)
{
    removeBond(bond(a, b));
}

/// Removes each bond in \p bonds from the molecule.
void Molecule::removeBonds(const std::vector<Bond *> &bonds)
{
    BOOST_REVERSE_FOREACH(Bond *bond, bonds){
        removeBond(bond);
    }
}

/// Returns a range containing all of the bonds in the molecule.
Molecule::BondRange Molecule::bonds() const
{
    return boost::make_iterator_range(d->bonds.begin(), d->bonds.end());
}

/// Returns the number of bonds in the molecule.
size_t Molecule::bondCount() const
{
    return d->bonds.size();
}

/// Returns the bond at index.
Bond* Molecule::bond(size_t index) const
{
    return d->bonds[index];
}

/// Returns the bond between atom \p a and \p b. Returns \c 0 if they
/// are not bonded.
///
/// To create a new bond between the atoms use Molecule::addBond().
Bond* Molecule::bond(const Atom *a, const Atom *b) const
{
    return a->bondTo(b);
}

/// Returns the bond between the atoms with indices \p a and \p b.
Bond* Molecule::bond(size_t a, size_t b) const
{
    return bond(atom(a), atom(b));
}

/// Requests that the bond capacity for the molecule be increased to
/// \p capacity.
///
/// \internal
void Molecule::setBondCapacity(size_t capacity)
{
    d->bonds.reserve(capacity);
    d->bondOrders.reserve(capacity);
    d->bondAtoms.reserve(capacity);
}

/// Returns the bond capacity for the molecule.
///
/// \internal
size_t Molecule::bondCapacity() const
{
    return d->bonds.capacity();
}

/// Returns \c true if the molecule contains bond.
bool Molecule::contains(const Bond *bond) const
{
    return contains(bond->atom1());
}

/// Removes all atoms and bonds from the molecule.
void Molecule::clear()
{
    removeBonds(d->bonds);
    removeAtoms(m_atoms);
}

// --- Ring Perception ----------------------------------------------------- //
/// Returns the ring at \p index.
///
/// Equivalent to calling:
/// \code
/// molecule.rings()[index];
/// \endcode
Ring* Molecule::ring(size_t index) const
{
    return rings()[index];
}

/// Returns a range containing all of the rings in the molecule.
///
/// This method implements the
/// \blueobeliskalgorithm{findSmallestSetOfSmallestRings}.
///
/// \warning The range of rings returned from this method is only
///          valid as long as the molecule's structure remains
///          unchanged. If any atoms or bonds in the molecule are
///          added or removed the old results must be discarded and
///          this method must be called again.
Molecule::RingRange Molecule::rings() const
{
    // only run ring perception if necessary
    if(!ringsPerceived()){
        // find rings
        foreach(const std::vector<Atom *> &ring, chemkit::algorithm::rppath(this)){
            d->rings.push_back(new Ring(ring));
        }

        // set perceived to true
        setRingsPerceived(true);
    }

    return boost::make_iterator_range(d->rings.begin(), d->rings.end());
}

/// Returns the number of rings in the molecule.
size_t Molecule::ringCount() const
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
Fragment* Molecule::fragment(size_t index) const
{
    return fragments()[index];
}

/// Returns a range containing all of the fragments in the molecule.
///
/// \warning The range of fragments returned from this method is only
///          valid as long as the molecule's structure remains
///          unchanged. If any atoms or bonds in the molecule are
///          added or removed the old results must be discarded and
///          this method must be called again.
Molecule::FragmentRange Molecule::fragments() const
{
    if(!fragmentsPerceived()){
        perceiveFragments();

        setFragmentsPerceived(true);
    }

    return boost::make_iterator_range(d->fragments.begin(),
                                      d->fragments.end());
}

/// Returns the number of fragments in the molecule.
size_t Molecule::fragmentCount() const
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
    removeAtoms(fragment->atoms());
}

Fragment* Molecule::fragmentForAtom(const Atom *atom) const
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
    Bitset unvisited(m_atoms.size());
    unvisited.set();

    for(;;){
        // bitset marking the atoms contained in the fragment
        Bitset bitset(m_atoms.size());

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
        if(position == Bitset::npos){
            break;
        }
    }
}

// --- Coordinates --------------------------------------------------------- //
/// Returns the coordinates for the molecule.
CartesianCoordinates* Molecule::coordinates() const
{
    if(!m_coordinates){
        if(d->coordinateSets.empty() ||
           d->coordinateSets.front()->type() == CoordinateSet::None){
            // create a new, empty cartesian coordinate set
            m_coordinates = new CartesianCoordinates(atomCount());
            d->coordinateSets.push_back(boost::make_shared<CoordinateSet>(m_coordinates));
        }
        else{
            const boost::shared_ptr<CoordinateSet> &coordinateSet = d->coordinateSets.front();

            switch(coordinateSet->type()){
                case CoordinateSet::Cartesian:
                    m_coordinates = coordinateSet->cartesianCoordinates();
                    break;
                case CoordinateSet::Internal:
                    m_coordinates = coordinateSet->internalCoordinates()->toCartesianCoordinates();
                    break;
                case CoordinateSet::Diagram:
                    m_coordinates = coordinateSet->diagramCoordinates()->toCartesianCoordinates();
                    break;
                default:
                    break;
            }
        }
    }

    return m_coordinates;
}

/// Add \p coordinates to the molecule.
void Molecule::addCoordinateSet(const boost::shared_ptr<CoordinateSet> &coordinates)
{
    d->coordinateSets.push_back(coordinates);
}

/// Add a new coordinate set containing \p coordinates.
void Molecule::addCoordinateSet(CartesianCoordinates *coordinates)
{
    addCoordinateSet(boost::make_shared<CoordinateSet>(coordinates));
}

/// Add a new coordinate set containing \p coordinates.
void Molecule::addCoordinateSet(InternalCoordinates *coordinates)
{
    addCoordinateSet(boost::make_shared<CoordinateSet>(coordinates));
}

/// Add a new coordinate set containing \p coordinates.
void Molecule::addCoordinateSet(DiagramCoordinates *coordinates)
{
    addCoordinateSet(boost::make_shared<CoordinateSet>(coordinates));
}

/// Removes \p coordinates from the molecule. Returns \c true if
/// successful.
bool Molecule::removeCoordinateSet(const boost::shared_ptr<CoordinateSet> &coordinates)
{
    typedef std::vector<boost::shared_ptr<CoordinateSet> >::iterator CoordinateSetIterator;
    CoordinateSetIterator iter = std::find(d->coordinateSets.begin(),
                                           d->coordinateSets.end(),
                                           coordinates);

    if(iter != d->coordinateSets.end()){
        d->coordinateSets.erase(iter);
        return true;
    }

    return false;
}

/// Returns the coordinate set at \p index in the molecule.
///
/// Equivalent to:
/// \code
/// molecule.coordinateSets()[index];
/// \endcode
boost::shared_ptr<CoordinateSet> Molecule::coordinateSet(size_t index) const
{
    assert(index < d->coordinateSets.size());

    return d->coordinateSets[index];
}

/// Returns the first coordinate set in the molecule of the given
/// \p type. Returns null if the molecule contains no coordinate
/// sets of the given \p type.
boost::shared_ptr<CoordinateSet> Molecule::coordinateSet(CoordinateSet::Type type) const
{
    foreach(const boost::shared_ptr<CoordinateSet> &coordinates, d->coordinateSets){
        if(coordinates->type() == type){
            return coordinates;
        }
    }

    return boost::shared_ptr<CoordinateSet>();
}

/// Returns a range containing all of the coordinate sets that the
/// molecule contains.
Molecule::CoordinateSetRange Molecule::coordinateSets() const
{
    return boost::make_iterator_range(d->coordinateSets.begin(),
                                      d->coordinateSets.end());
}

/// Returns the number of coordinate sets stored in the molecule.
///
/// Equivalent to:
/// \code
/// molecule.coordinateSets().size();
/// \endcode
size_t Molecule::coordinateSetCount() const
{
    return d->coordinateSets.size();
}

// --- Geometry ------------------------------------------------------------ //
/// Returns the distance between atoms \p a and \p b. The returned
/// distance is in Angstroms.
Real Molecule::distance(const Atom *a, const Atom *b) const
{
    return coordinates()->distance(a->index(), b->index());
}

/// Returns the angle between atoms \p a, \p b, and \p c. The
/// returned angle is in degrees.
Real Molecule::bondAngle(const Atom *a, const Atom *b, const Atom *c) const
{
    return coordinates()->angle(a->index(), b->index(), c->index());
}

/// Returns the torsion angle (also known as the dihedral angle)
/// between atoms \p a, \p b, \p c, and \p d. The returned angle is
/// in degrees.
Real Molecule::torsionAngle(const Atom *a, const Atom *b, const Atom *c, const Atom *d) const
{
    return coordinates()->torsionAngle(a->index(), b->index(), c->index(), d->index());
}

/// Returns the wilson angle between the plane made by atoms \p a,
/// \p b, \p c and the vector from \p c to \p d. The returned angle
/// is in degrees.
Real Molecule::wilsonAngle(const Atom *a, const Atom *b, const Atom *c, const Atom *d) const
{
    return coordinates()->wilsonAngle(a->index(), b->index(), c->index(), d->index());
}

/// Moves all of the atoms in the molecule so that the center point
/// is at \p position.
void Molecule::setCenter(const Point3 &position)
{
    const Vector3 &vector = position - center();

    foreach(Atom *atom, m_atoms){
        atom->setPosition(atom->position() + vector);
    }
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
    if(!m_coordinates){
        return Point3(0, 0, 0);
    }

    return m_coordinates->center();
}

/// Returns the center of mass for the molecule.
///
/// This method implements the \blueobeliskalgorithm{calculate3DCenterOfMass}.
Point3 Molecule::centerOfMass() const
{
    if(!m_coordinates){
        return Point3(0, 0, 0);
    }

    std::vector<Real> weights;

    foreach(const Atom *atom, m_atoms){
        weights.push_back(atom->mass());
    }

    return m_coordinates->weightedCenter(weights);
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

            if(bond->stereochemistry() != Stereochemistry::None){
                newBond->setStereochemistry(bond->stereochemistry());
            }
        }
    }

    return *this;
}

/// Returns the atom at \p index in the molecule.
///
/// Equivalent to:
/// \code
/// molecule.atom(index);
/// \endcode
Atom* Molecule::operator[](size_t index) const
{
    return atom(index);
}

// --- Internal Methods ---------------------------------------------------- //
void Molecule::notifyWatchers(MoleculeWatcher::ChangeType type)
{
    foreach(MoleculeWatcher *watcher, d->watchers){
        watcher->moleculeChanged(this, type);
    }
}

void Molecule::notifyWatchers(const Atom *atom, MoleculeWatcher::ChangeType type)
{
    foreach(MoleculeWatcher *watcher, d->watchers){
        watcher->atomChanged(atom, type);
    }
}

void Molecule::notifyWatchers(const Bond *bond, MoleculeWatcher::ChangeType type)
{
    foreach(MoleculeWatcher *watcher, d->watchers){
        watcher->bondChanged(bond, type);
    }
}

void Molecule::addWatcher(MoleculeWatcher *watcher) const
{
    d->watchers.push_back(watcher);
}

void Molecule::removeWatcher(MoleculeWatcher *watcher) const
{
    d->watchers.erase(std::remove(d->watchers.begin(), d->watchers.end(), watcher));
}

Stereochemistry* Molecule::stereochemistry()
{
    if(!m_stereochemistry){
        m_stereochemistry = new Stereochemistry(this);
    }

    return m_stereochemistry;
}

} // end chemkit namespace
