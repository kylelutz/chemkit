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

#include "molecule.h"

#include "atom.h"
#include "bond.h"
#include "ring.h"
#include "point3.h"
#include "element.h"
#include "vector3.h"
#include "constants.h"
#include "lineformat.h"
#include "quaternion.h"
#include "coordinates.h"
#include "moleculargraph.h"
#include "moleculewatcher.h"
#include "internalcoordinates.h"
#include "moleculardescriptor.h"

namespace chemkit {

// === MoleculePrivate ===================================================== //
class MoleculePrivate
{
    public:
        MoleculePrivate();

        QString name;
        QList<Bond *> bonds;
        QList<Residue *> residues;
        QList<Conformer *> conformers;
        Conformer *conformer;
        bool ringsPerceived;
        QList<Ring *> rings;
        bool fragmentsPerceived;
        QList<Fragment *> fragments;
        QList<MoleculeWatcher *> watchers;
};

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
/// ChemicalFile class.
///
/// Molecule objects take ownership of all the Atom, Bond, Residue,
/// Ring, Fragment, and Conformer objects that they contain. Deleting
/// the molecule will also delete all of the objects that it
/// contains.

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
Molecule::Molecule(const QString &formula, const QString &format)
    : d(new MoleculePrivate)
{
    LineFormat *lineFormat = LineFormat::create(format);
    if(!lineFormat){
        return;
    }

    lineFormat->read(formula, this);

    delete lineFormat;
}

/// Creates a new molecule that is a copy of \p molecule.
Molecule::Molecule(const Molecule &molecule)
    : d(new MoleculePrivate)
{
    d->name = molecule.name();

    QHash<const Atom *, Atom *> oldToNew;

    foreach(const Atom *atom, molecule.atoms()){
        Atom *newAtom = addAtomCopy(atom);
        oldToNew[atom] = newAtom;
    }

    foreach(const Bond *bond, molecule.bonds()){
        Bond *newBond = addBond(oldToNew[bond->atom1()], oldToNew[bond->atom2()]);
        newBond->setOrder(bond->order());
    }
}

/// Destroys a molecule. This also destroys all of the atoms, bonds,
/// and residues that the molecule contains.
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

    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Sets the name of the molecule.
void Molecule::setName(const QString &name)
{
    d->name = name;
    notifyObservers(NameChanged);
}

/// Returns the name of the molecule.
QString Molecule::name() const
{
    return d->name;
}

/// Returns the chemical formula (e.g. "H2O") for the molecule. The
/// formula is formated according to the Hill system.
QString Molecule::formula() const
{
    // a map of atomic symbols to their quantity
    QMap<QString, int> composition;
    foreach(const Atom *atom, m_atoms){
        composition[atom->symbol()]++;
    }

    QString formula;

    if(composition.contains("C")){
        formula += "C";
        if(composition["C"] > 1)
            formula += QString::number(composition["C"]);
        composition.remove("C");

        if(composition.contains("H")){
            formula += "H";
            if(composition["H"] > 1)
                formula += QString::number(composition["H"]);
        }
        composition.remove("H");
    }

    QMapIterator<QString, int> iter(composition);

    while(iter.hasNext()){
        iter.next();

        formula += iter.key();
        if(iter.value() > 1)
            formula += QString::number(iter.value());
    }

    return formula;
}

/// Returns the the formula of the molecule using the specified
/// format. Returns a null QString if format is not supported or if
/// an error occurs.
///
/// The following example returns the InChI formula for a molecule:
/// \code
/// molecule->formula("inchi");
/// \endcode
///
/// \see LineFormat
QString Molecule::formula(const QString &format) const
{
    LineFormat *lineFormat = LineFormat::create(format);
    if(!lineFormat)
        return QString();

    QString formula = lineFormat->write(this);

    delete lineFormat;

    return formula;
}

/// Calculates and returns the molecular descriptor \p name. If the
/// descriptor is not available or the calculation fails a null
/// QVariant is returned.
///
/// For example, to calculate the Randic index of the molecule use:
/// \code
/// double randicIndex = molecule->descriptor("randic-index").toDouble();
/// \endcode
///
/// \see MolecularDescriptor
QVariant Molecule::descriptor(const QString &name) const
{
    QScopedPointer<MolecularDescriptor> descriptor(MolecularDescriptor::create(name));
    if(!descriptor){
        return QVariant();
    }

    return descriptor->value(this);
}

/// Returns the total molar mass of the molecule. Mass is in g/mol.
Float Molecule::mass() const
{
    Float mass = 0;

    foreach(const Atom *atom, m_atoms)
        mass += atom->mass();

    return mass;
}

// --- Structure ----------------------------------------------------------- //
/// Adds a new atom of the given \p element to the molecule. If
/// the element is invalid 0 will be returned.
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
    if(!element.isValid()){
        return 0;
    }

    Atom *atom = new Atom(this, element);
    m_atoms.append(atom);

    setFragmentsPerceived(false);
    notifyObservers(atom, AtomAdded);

    return atom;
}

/// Adds a new atom to the molecule. The new atom will have the same
/// properties as atom (atomic number, mass number, etc).
Atom* Molecule::addAtomCopy(const Atom *atom)
{
    Atom *newAtom = addAtom(atom->atomicNumber());

    newAtom->setMassNumber(atom->massNumber());
    newAtom->setPartialCharge(atom->partialCharge());
    newAtom->setPosition(atom->position());
    newAtom->setChirality(atom->chirality());

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
    foreach(Bond *bond, atom->bonds()){
        removeBond(bond);
    }

    Q_ASSERT(m_atoms.contains(atom));
    m_atoms.removeOne(atom);

    atom->m_molecule = 0;
    notifyObservers(atom, AtomRemoved);

    delete atom;
}

/// Returns the number of atoms in the molecule of the given
/// \p element.
int Molecule::atomCount(const Element &element) const
{
    int count = 0;

    foreach(Atom *atom, m_atoms){
        if(atom->is(element)){
            count++;
        }
    }

    return count;
}

/// Returns the index of \p atom in the molecule.
int Molecule::indexOf(const Atom *atom) const
{
    return m_atoms.indexOf(const_cast<Atom *>(atom));
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
    // ensure that this molecule contains both atoms
    if(!contains(a) || !contains(b) || a->isBondedTo(b)){
        return 0;
    }

    // check to see if they are already bonded
    if(a->isBondedTo(b)){
        return bond(a, b);
    }

    Bond *bond = new Bond(a, b, order);

    bond->atom1()->addBond(bond);
    bond->atom2()->addBond(bond);
    d->bonds.append(bond);

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
    bool found = d->bonds.removeOne(bond);
    if(!found){
        return;
    }

    bond->atom1()->removeBond(bond);
    bond->atom2()->removeBond(bond);

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
QList<Bond *> Molecule::bonds() const
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
    return d->bonds.value(index, 0);
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

/// Returns the index of \p bond in the molecule.
int Molecule::indexOf(const Bond *bond) const
{
    return d->bonds.indexOf(const_cast<Bond *>(bond));
}

/// Returns \c true if the molecule contains bond.
bool Molecule::contains(const Bond *bond) const
{
    return contains(bond->atom1());
}

/// Adds \p residue to the molecule.
void Molecule::addResidue(Residue *residue)
{
    if(residue->molecule() != this){
        return;
    }

    d->residues.append(residue);

    foreach(Atom *atom, residue->atoms()){
        atom->setResidue(residue);
    }
}

/// Removes \p residue from the molecule.
void Molecule::removeResidue(Residue *residue)
{
    if(residue->molecule() != this){
        return;
    }

    d->residues.removeAll(residue);

    foreach(Atom *atom, residue->atoms()){
        atom->setResidue(0);
    }
}

/// Returns a list of all residues in the molecule.
QList<Residue *> Molecule::residues() const
{
    return d->residues;
}

/// Returns the number of residues in the molecule.
int Molecule::residueCount() const
{
    return residues().size();
}

/// Removes all atoms, bonds, and residues from the molecule.
void Molecule::clear()
{
    foreach(Bond *bond, d->bonds){
        removeBond(bond);
    }
    foreach(Atom *atom, m_atoms){
        removeAtom(atom);
    }
    foreach(Residue *residue, d->residues){
        removeResidue(residue);
    }
}

// --- Comparison ---------------------------------------------------------- //
/// Returns \c true if the molecule equals \p molecule.
bool Molecule::equals(const Molecule *molecule, CompareFlags flags) const
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
bool Molecule::contains(const Molecule *molecule, CompareFlags flags) const
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

    return !molecule->mapping(this, flags).isEmpty();
}

/// Returns \c true if the molecule is a substructure of \p molecule.
bool Molecule::isSubstructureOf(const Molecule *molecule, CompareFlags flags) const
{
    return molecule->contains(this, flags);
}

/// Returns a mapping (also known as an isomorphism) between the atoms
/// in the molecule and the atoms in \p molecule.
AtomMapping Molecule::mapping(const Molecule *molecule, CompareFlags flags) const
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
        for(int i = 0; i < source->bondCount(); i++){
            const Bond *bond = source->bond(i);

            if(bond->isAromatic()){
                source->setBondLabel(i, aromaticBondLabel);
            }
        }

        for(int i = 0; i < target->bondCount(); i++){
            const Bond *bond = target->bond(i);

            if(bond->isAromatic()){
                target->setBondLabel(i, aromaticBondLabel);
            }
        }
    }

    AtomMapping mapping = MolecularGraph::isomorphism(source, target);

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
Moiety Molecule::find(const Molecule *moiety, CompareFlags flags) const
{
    AtomMapping mapping = moiety->mapping(this, flags);

    // no mapping found, return empty moiety
    if(mapping.isEmpty()){
        return Moiety();
    }

    QList<const Atom *> moietyAtoms;
    foreach(Atom *atom, moiety->atoms()){
        moietyAtoms.append(mapping.map(atom));
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
    return rings().value(index, 0);
}

/// Returns a list of all rings in the molecule.
///
/// \warning The list of rings returned from this method is only
///          valid as long as the molecule's structure remains
///          unchanged. If any atoms or bonds in the molecule are
///          added or removed the old results must be discarded and
///          this method must be called again.
QList<Ring *> Molecule::rings() const
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
    return fragments().value(index, 0);
}

/// Returns a list of fragments in the molecule.
///
/// \warning The list of fragments returned from this method is only
///          valid as long as the molecule's structure remains
///          unchanged. If any atoms or bonds in the molecule are
///          added or removed the old results must be discarded and
///          this method must be called again.
QList<Fragment *> Molecule::fragments() const
{
    if(!fragmentsPerceived()){
        foreach(Atom *atom, m_atoms){
            if(!atom->m_fragment){
                d->fragments.append(new Fragment(atom));
            }
        }

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
    QList<Atom *> atoms = fragment->atoms();

    foreach(Atom *atom, atoms){
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

        foreach(Atom *atom, m_atoms){
            atom->m_fragment = 0;
        }
    }

    d->fragmentsPerceived = perceived;
}

bool Molecule::fragmentsPerceived() const
{
    return d->fragmentsPerceived;
}

// --- Geometry ------------------------------------------------------------ //
/// Sets the coordinates for the atoms in the molecule to
/// \p coordinates.
void Molecule::setCoordinates(const Coordinates *coordinates)
{
    int size = qMin(m_atoms.size(), coordinates->size());

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
Float Molecule::distance(const Atom *a, const Atom *b) const
{
    return Point3::distance(a->position(), b->position());
}

/// Returns the angle between atoms \p a, \p b, and \p c. The
/// returned angle is in degrees.
Float Molecule::bondAngle(const Atom *a, const Atom *b, const Atom *c) const
{
    return Point3::angle(a->position(),
                         b->position(),
                         c->position());
}

/// Returns the torsion angle (also known as the dihedral angle)
/// between atoms \p a, \p b, \p c, and \p d. The returned angle is
/// in degrees.
Float Molecule::torsionAngle(const Atom *a, const Atom *b, const Atom *c, const Atom *d) const
{
    return Point3::torsionAngle(a->position(),
                                b->position(),
                                c->position(),
                                d->position());
}

/// Returns the wilson angle between the plane made by atoms \p a,
/// \p b, \p c and the vector from \p c to \p d. The returned angle
/// is in degrees.
Float Molecule::wilsonAngle(const Atom *a, const Atom *b, const Atom *c, const Atom *d) const
{
    return Point3::wilsonAngle(a->position(),
                               b->position(),
                               c->position(),
                               d->position());
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
void Molecule::setCenter(Float x, Float y, Float z)
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
        return Point3();
    }

    // sums for each component
    Float sx = 0;
    Float sy = 0;
    Float sz = 0;

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
    Float sx = 0;
    Float sy = 0;
    Float sz = 0;

    // sum of weights
    Float sw = 0;

    foreach(const Atom *atom, m_atoms){
        Float w = atom->mass();

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
void Molecule::moveBy(Float dx, Float dy, Float dz)
{
    foreach(Atom *atom, m_atoms){
        atom->moveBy(dx, dy, dz);
    }
}

/// Rotates the positions of all the atoms in the molecule
/// by \p angle degrees around \p axis.
void Molecule::rotate(const Vector3 &axis, Float angle)
{
    foreach(Atom *atom, m_atoms){
        atom->setPosition(Quaternion::rotate(atom->position(), axis, angle));
    }
}

/// Returns \c true if the molecule has coordinates for any of the
/// atoms.
bool Molecule::hasCoordinates() const
{
    foreach(const Atom *atom, m_atoms){
        if(!atom->position().isNull()){
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
    d->conformers.append(conformer);
    return conformer;
}

/// Removes \p conformer from the molecule. The currently active
/// conformer cannot be removed.
void Molecule::removeConformer(Conformer *conformer)
{
    if(!d->conformers.contains(conformer)){
        return;
    }

    // forbid removal of the currently active conformer
    if(conformer == d->conformer){
        return;
    }

    d->conformers.removeOne(conformer);

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
    return conformers().value(index, 0);
}

/// Returns a list of all conformers in the molecule.
QList<Conformer *> Molecule::conformers() const
{
    if(d->conformers.isEmpty()){
        d->conformer = new Conformer(this);
        d->conformers.append(d->conformer);
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

        QHash<const Atom *, Atom *> oldToNew;

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
QList<Atom *> Molecule::atomPathBetween(const Atom *a, const Atom *b) const
{
    if(a == b){
        return QList<Atom *>();
    }
    else if(!a->isConnectedTo(b)){
        return QList<Atom *>();
    }
    else if(a->isBondedTo(b)){
        QList<Atom *> path;
        path.append(const_cast<Atom *>(b));
        return path;
    }

    QBitArray visited(atomCount());
    visited.setBit(a->index());

    QList<QList<Atom *> > paths;

    foreach(Atom *neighbor, a->neighbors()){
        visited.setBit(neighbor->index());
        QList<Atom *> path;
        path.append(neighbor);
        paths.append(path);
    }

    while(paths.size()){
        QList<Atom *> path = paths.takeFirst();

        const Atom *lastAtom = path.last();
        if(lastAtom == b){
            return path;
        }
        else{
            foreach(Atom *neighbor, lastAtom->neighbors()){
                if(visited[neighbor->index()])
                    continue;

                visited.setBit(neighbor->index());
                QList<Atom *> nextPath(path);
                nextPath.append(neighbor);
                paths.append(nextPath);
            }
        }
    }

    return QList<Atom *>();
}

int Molecule::atomCountBetween(const Atom *a, const Atom *b) const
{
    return atomPathBetween(a, b).size();
}

int Molecule::atomCountBetween(const Atom *a, const Atom *b, int maxCount) const
{
    int count = atomCountBetween(a, b);

    if(count > maxCount)
        return 0;

    return count;
}

QList<Bond *> Molecule::bondPathBetween(const Atom *a, const Atom *b) const
{
    QList<Atom *> atomPath = atomPathBetween(a, b);
    if(atomPath.isEmpty()){
        return QList<Bond *>();
    }

    QList<Bond *> bondPath;
    bondPath.append(a->bondTo(atomPath[0]));
    for(int i = 0; i < atomPath.size()-1; i++){
        bondPath.append(atomPath[i]->bondTo(atomPath[i+1]));
    }

    return bondPath;
}

int Molecule::bondCountBetween(const Atom *a, const Atom *b) const
{
    return bondPathBetween(a, b).size();
}

int Molecule::bondCountBetween(const Atom *a, const Atom *b, int maxCount) const
{
    int count = bondCountBetween(a, b);

    if(count > maxCount){
        return 0;
    }

    return count;
}

void Molecule::notifyObservers(ChangeType type)
{
    foreach(MoleculeWatcher *watcher, d->watchers){
        watcher->notifyObservers(this, type);
    }
}

void Molecule::notifyObservers(const Atom *atom, ChangeType type)
{
    foreach(MoleculeWatcher *watcher, d->watchers){
        watcher->notifyObservers(atom, type);
    }
}

void Molecule::notifyObservers(const Bond *bond, ChangeType type)
{
    foreach(MoleculeWatcher *watcher, d->watchers){
        watcher->notifyObservers(bond, type);
    }
}

void Molecule::addWatcher(MoleculeWatcher *watcher) const
{
    d->watchers.append(watcher);
}

void Molecule::removeWatcher(MoleculeWatcher *watcher) const
{
    d->watchers.removeOne(watcher);
}

bool Molecule::isSubsetOf(const Molecule *molecule, CompareFlags flags) const
{
    Q_UNUSED(flags);

    QList<Atom *> otherAtoms = molecule->atoms();

    foreach(const Atom *atom, m_atoms){
        bool found = false;

        foreach(Atom *otherAtom, otherAtoms){
            if(atom->atomicNumber() == otherAtom->atomicNumber()){
                otherAtoms.removeOne(otherAtom);
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

} // end chemkit namespace
