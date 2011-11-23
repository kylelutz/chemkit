#ifndef CHEMKIT_FORCEFIELDENERGYDESCRIPTOR_H
#define CHEMKIT_FORCEFIELDENERGYDESCRIPTOR_H

#include "md.h"

#include <chemkit/moleculardescriptor.h>

namespace chemkit {
namespace md {

template<typename ForceField>
class ForceFieldEnergyDescriptor : public MolecularDescriptor
{
public:
    // construction and destruction
    ForceFieldEnergyDescriptor(const std::string &name);
    ~ForceFieldEnergyDescriptor();

    // descriptor
    virtual Variant value(const Molecule *molecule) const;
};

} // end md namespace
} // end chemkit namespace

#include "forcefieldenergydescriptor-inline.h"

#endif // CHEMKIT_FORCEFIELDENERGYDESCRIPTOR_H
