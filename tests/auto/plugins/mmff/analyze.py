#!/usr/bin/python

# this script compares the mmff94.expected and mmff94.actual files
# and outputs the differences

import os
import sys
import xml.dom.minidom

COLORS_ENABLED = False 

class AtomResults:
    def __init__(self, type, charge):
        self.type = type
        self.charge = charge

class MoleculeResults:
    def __init__(self, name, energy):
        self.name = name
        self.energy = energy
        self.atoms = []

class ResultsFile:
    def __init__(self, fileName):
        self.fileName = fileName
        self.molecules = []

    def read(self):
        doc = xml.dom.minidom.parse(self.fileName)

        for moleculeElem in doc.getElementsByTagName('molecule'):
            name = moleculeElem.getAttribute('name')
            energy = float(moleculeElem.getAttribute('energy'))
            moleculeResults = MoleculeResults(name, energy)

            for atomElem in moleculeElem.getElementsByTagName('atom'):
                type = atomElem.getAttribute('type')
                charge = float(atomElem.getAttribute('charge'))
                moleculeResults.atoms.append(AtomResults(type, charge))

            self.molecules.append(moleculeResults)

if __name__ == '__main__':
    actualResultsFile = 'mmff94.actual'
    expectedResultsFile = 'mmff94.expected'

    if not os.path.exists(actualResultsFile):
        print 'could not find actual results file (%s)' % actualResultsFile
        sys.exit(-1)
    if not os.path.exists(expectedResultsFile):
        print 'could not find expected results file (%s)' % expectedResultsFile
        sys.exit(-1)

    actualResults = ResultsFile(actualResultsFile)
    actualResults.read()

    expectedResults = ResultsFile(expectedResultsFile)
    expectedResults.read()

    # escape codes to color text
    RED_COLOR = '\033[91m'
    END_COLOR = '\033[0m'
    if not COLORS_ENABLED:
        RED_COLOR = ''
        END_COLOR = ''

    ATOMS_FAILED = 0
    MOLECULES_FAILED = 0

    # compare files
    expectedMoleculeIndex = 0
    for i, actualMolecule in enumerate(actualResults.molecules):
        expectedMolecule = expectedResults.molecules[expectedMoleculeIndex]

        while expectedMolecule.name != actualMolecule.name:
            expectedMoleculeIndex += 1
            expectedMolecule = expectedResults.molecules[expectedMoleculeIndex]

        print '%i. %s' % (expectedMoleculeIndex+1, actualMolecule.name)
        
        for j in range(len(actualMolecule.atoms)):
            actualAtom = actualMolecule.atoms[j]
            expectedAtom = expectedMolecule.atoms[j]

            expectedTypeText = ''
            colorCode = ''
            if(actualAtom.type != expectedAtom.type or 
               (abs(actualAtom.charge - expectedAtom.charge) > 0.01)):
                ATOMS_FAILED += 1
                colorCode = RED_COLOR
                expectedTypeText = '%s[%s, %s] -- FAILED%s' % (colorCode, expectedAtom.type, expectedAtom.charge, END_COLOR)

            print '   %i. %s, %s %s' % (j+1, actualAtom.type, actualAtom.charge, expectedTypeText)

        colorCode = ''
        if(int(actualMolecule.energy) != int(expectedMolecule.energy)):
            MOLECULES_FAILED += 1
            colorCode = RED_COLOR

        print 'energy: %f %s[%f]%s' % (actualMolecule.energy, colorCode, expectedMolecule.energy, END_COLOR)

    # print some statistics
    print >> sys.stderr, ''
    print >> sys.stderr, 'atoms: %i failed' % ATOMS_FAILED
    print >> sys.stderr, 'molecules: %i failed' % MOLECULES_FAILED
