Molecule
========

.. currentmodule:: chemkit

.. autoclass:: chemkit.Molecule

    .. automethod:: setName(name)
    .. automethod:: name()
    .. automethod:: formula()
    .. automethod:: size()
    .. automethod:: isEmpty()
    .. automethod:: mass()
    .. automethod:: addAtom(element)
    .. automethod:: removeAtom(atom)
    .. automethod:: atom(index)
    .. automethod:: atoms()
    .. automethod:: atomCount()
    .. automethod:: addBond(a, b[, order = chemkit.Bond.Single])
    .. automethod:: removeBond(bond)
    .. automethod:: bond(index)
    .. automethod:: bonds()
    .. automethod:: bondCount()
    .. automethod:: ring(index)
    .. automethod:: rings()
    .. automethod:: ringCount()
    .. automethod:: fragment(index)
    .. automethod:: fragments()
    .. automethod:: fragmentCount()
    .. automethod:: isFragmented()
    .. automethod:: removeFragment(fragment)

