#ifndef CHEMKIT_FORCEFIELDENERGYDESCRIPTOR_INLINE_H
#define CHEMKIT_FORCEFIELDENERGYDESCRIPTOR_INLINE_H

#include "forcefieldenergydescriptor.h"

#include "forcefield.h"

namespace chemkit {

// === ForceFieldEnergyDescriptor ========================================== //
/// \class ForceFieldEnergyDescriptor forcefieldenergydescriptor.h chemkit/forcefieldenergydescriptor.h
/// \ingroup chemkit-md
/// \internal
/// \brief The ForceFieldEnergyDescriptor class provides a molecular
///        descriptor for force field energy.

// --- Construction and Destruction ---------------------------------------- //
/// Creates a new force field energy descriptor with \p name.
template<typename ForceField>
ForceFieldEnergyDescriptor<ForceField>::ForceFieldEnergyDescriptor(const std::string &name)
    : MolecularDescriptor(name)
{
}

/// Destroys the force field energy descriptor object.
template<typename ForceField>
ForceFieldEnergyDescriptor<ForceField>::~ForceFieldEnergyDescriptor()
{
}

// --- Descriptor ---------------------------------------------------------- //
/// Returns the energy value for \p molecule.
template<typename ForceField>
Variant ForceFieldEnergyDescriptor<ForceField>::value(const Molecule *molecule) const
{
    ForceField forceField;
    forceField.addMolecule(molecule);

    bool ok = forceField.setup();
    if(!ok){
        return Variant();
    }

    return forceField.energy();
}

} // end chemkit namespace

#endif // CHEMKIT_FORCEFIELDENERGYDESCRIPTOR_INLINE_H
