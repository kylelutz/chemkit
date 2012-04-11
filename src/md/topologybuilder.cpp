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

#include "topologybuilder.h"

#include <boost/make_shared.hpp>

#include <chemkit/atom.h>
#include <chemkit/bond.h>
#include <chemkit/foreach.h>
#include <chemkit/molecule.h>
#include <chemkit/atomtyper.h>
#include <chemkit/partialchargemodel.h>

#include "topology.h"

namespace chemkit {

namespace {

bool atomsWithinTwoBonds(const Atom *a, const Atom *b)
{
    foreach(const Atom *neighbor, a->neighbors()){
        if(neighbor == b){
            return true;
        }
        else if(neighbor->isBondedTo(b)){
            return true;
        }
    }

    return false;
}

} // end anonymous namespace

// === TopologyBuilderPrivate ============================================== //
class TopologyBuilderPrivate
{
public:
    std::string atomTyper;
    std::string partialChargeModel;
    boost::shared_ptr<Topology> topology;
};

// === TopologyBuilder ===================================================== //
/// \class TopologyBuilder topologybuilder.h chemkit/topologybuilder.h
/// \ingroup chemkit-md
/// \brief The TopologyBuilder class builds molecular dynamics
///        topologies.
///
/// For example, to create a topology suitable for use with the UFF
/// force field:
/// \code
/// // load molecule from file, formula, etc.
/// const Molecule *molecule = ...
///
/// // create topology builder object
/// TopologyBuilder builder;
///
/// // set options
/// builder.setAtomTyper("uff");
///
/// // add molecule
/// builder.addMolecule(molecule);
///
/// // get topology
/// boost::shared_ptr<Topology> = builder.topology();
/// \endcode

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new topology builder.
TopologyBuilder::TopologyBuilder()
    : d(new TopologyBuilderPrivate)
{
    d->topology = boost::make_shared<Topology>();
}

/// Destroys the topology builder object.
TopologyBuilder::~TopologyBuilder()
{
    delete d;
}

// --- Properties ---------------------------------------------------------- //
/// Returns the current size of the topology.
size_t TopologyBuilder::size() const
{
    return d->topology->size();
}

/// Returns \c true if the current topology is empty.
bool TopologyBuilder::isEmpty() const
{
    return d->topology->isEmpty();
}

/// Sets the atom typer to \p atomTyper.
bool TopologyBuilder::setAtomTyper(const std::string &atomTyper)
{
    d->atomTyper = atomTyper;

    return true;
}

/// Sets the partial charge model to \p model.
bool TopologyBuilder::setPartialChargeModel(const std::string &model)
{
    d->partialChargeModel = model;

    return true;
}

// --- Topology ------------------------------------------------------------ //
/// Adds \p molecule to the topology.
void TopologyBuilder::addMolecule(const Molecule *molecule)
{
    // reference to the current topology
    boost::shared_ptr<Topology> &topology = d->topology;

    // create atom typer
    boost::scoped_ptr<AtomTyper> atomTyper(0);
    if(!d->atomTyper.empty()){
        atomTyper.reset(AtomTyper::create(d->atomTyper));
    }

    // create partial charge model
    boost::scoped_ptr<PartialChargeModel> chargeModel(0);
    if(!d->partialChargeModel.empty()){
        chargeModel.reset(PartialChargeModel::create(d->partialChargeModel));
    }

    // increase topology size to store the new molecule
    size_t initialSize = topology->size();
    topology->resize(initialSize + molecule->size());

    // set atom types
    if(atomTyper){
        atomTyper->setMolecule(molecule);

        foreach(const Atom *atom, molecule->atoms()){
            topology->setType(initialSize + atom->index(), atomTyper->type(atom));
        }
    }

    // set atom masses
    foreach(const Atom *atom, molecule->atoms()){
        topology->setMass(initialSize + atom->index(), atom->mass());
    }

    // set atom charges
    if(chargeModel){
        chargeModel->setMolecule(molecule);

        foreach(const Atom *atom, molecule->atoms()){
            topology->setCharge(initialSize + atom->index(), chargeModel->partialCharge(atom));
        }
    }
    else{
        foreach(const Atom *atom, molecule->atoms()){
            topology->setCharge(initialSize + atom->index(), atom->partialCharge());
        }
    }

    // add bonded interactions
    foreach(const Bond *bond, molecule->bonds()){
        topology->addBondedInteraction(initialSize + bond->atom1()->index(),
                                       initialSize + bond->atom2()->index());

        if(atomTyper){
            int type = atomTyper->bondedInteractionType(bond->atom1(), bond->atom2());
            if(type != 0){
                topology->setBondedInteractionType(bond->atom1()->index(),
                                                   bond->atom2()->index(),
                                                   type);
            }
        }
    }

    // add angle interactions
    foreach(const Atom *atom, molecule->atoms()){
        if(!atom->isTerminal()){
            std::vector<const Atom *> neighbors(atom->neighbors().begin(),
                                                atom->neighbors().end());

            for(size_t i = 0; i < neighbors.size(); i++){
                for(size_t j = i + 1; j < neighbors.size(); j++){
                    topology->addAngleInteraction(initialSize + neighbors[i]->index(),
                                                  initialSize + atom->index(),
                                                  initialSize + neighbors[j]->index());

                    if(atomTyper){
                        int type = atomTyper->angleInteractionType(neighbors[i], atom, neighbors[j]);
                        if(type != 0){
                            topology->setAngleInteractionType(initialSize + neighbors[i]->index(),
                                                              initialSize + atom->index(),
                                                              initialSize + neighbors[j]->index(),
                                                              type);
                        }
                    }
                }
            }
        }
    }

    // add torsion interactions
    std::vector<std::pair<const Atom *, const Atom *> > torsionPairs;
    foreach(const Bond *bond, molecule->bonds()){
        if(!bond->atom1()->isTerminal() && !bond->atom2()->isTerminal()){
            torsionPairs.push_back(std::make_pair(bond->atom1(), bond->atom2()));
        }
    }

    std::pair<const Atom *, const Atom *> torsionPair;
    foreach(torsionPair, torsionPairs){
        const Atom *b = torsionPair.first;
        const Atom *c = torsionPair.second;

        foreach(const Atom *a, b->neighbors()){
            if(a == c){
                continue;
            }

            foreach(const Atom *d, c->neighbors()){
                if(d == b || d == a){
                    continue;
                }

                topology->addTorsionInteraction(initialSize + a->index(),
                                                initialSize + b->index(),
                                                initialSize + c->index(),
                                                initialSize + d->index());

                if(atomTyper){
                    int type = atomTyper->torsionInteractionType(a, b, c, d);
                    if(type != 0){
                        topology->setTorsionInteractionType(initialSize + a->index(),
                                                            initialSize + b->index(),
                                                            initialSize + c->index(),
                                                            initialSize + d->index(),
                                                            type);
                    }
                }
            }
        }
    }

    // add improper torsion interactions
    foreach(const Atom *atom, molecule->atoms()){
        if(atom->neighborCount() == 3){
            topology->addImproperTorsionInteraction(initialSize + atom->neighbor(0)->index(),
                                                    initialSize + atom->index(),
                                                    initialSize + atom->neighbor(1)->index(),
                                                    initialSize + atom->neighbor(2)->index());
        }
    }

    // add nonbonded interactions
    std::vector<const Atom *> atoms(molecule->atoms().begin(), molecule->atoms().end());
    for(size_t i = 0; i < atoms.size(); i++){
        for(size_t j = i + 1; j < atoms.size(); j++){
            if(!atomsWithinTwoBonds(atoms[i], atoms[j])){
                topology->addNonbondedInteraction(initialSize + atoms[i]->index(),
                                                  initialSize + atoms[j]->index());
            }
        }
    }
}

/// Returns the constructed topology.
boost::shared_ptr<Topology> TopologyBuilder::topology() const
{
    return d->topology;
}

} // end chemkit namespace
