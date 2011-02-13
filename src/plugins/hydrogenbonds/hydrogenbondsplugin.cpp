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

#include "hydrogenbondsplugin.h"

#include "hydrogenbonddonorsdescriptor.h"
#include "hydrogenbondacceptorsdescriptor.h"

HydrogenBondsPlugin::HydrogenBondsPlugin()
    : chemkit::Plugin("hydrogenbonds")
{
    chemkit::MolecularDescriptor::registerDescriptor("hydrogen-bond-donors", createDonorsDescriptor);
    chemkit::MolecularDescriptor::registerDescriptor("hydrogen-bond-acceptors", createAcceptorsDescriptor);
}

HydrogenBondsPlugin::~HydrogenBondsPlugin()
{
    chemkit::MolecularDescriptor::unregisterDescriptor("hydrogen-bond-donors", createDonorsDescriptor);
    chemkit::MolecularDescriptor::unregisterDescriptor("hydrogen-bond-acceptors", createAcceptorsDescriptor);
}

chemkit::MolecularDescriptor* HydrogenBondsPlugin::createDonorsDescriptor()
{
    return new HydrogenBondDonorsDescriptor;
}

chemkit::MolecularDescriptor* HydrogenBondsPlugin::createAcceptorsDescriptor()
{
    return new HydrogenBondAcceptorsDescriptor;
}

Q_EXPORT_PLUGIN2(hydrogenbonds, HydrogenBondsPlugin);
