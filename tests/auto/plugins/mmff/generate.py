#!/usr/bin/python

# this script parses the MMFF94_opti.log file from the MMFF94 validation
# suite and creates the mmff94.expected file.
# http://www.ccl.net/cca/data/MMFF94/MMFF94_opti.log

class Molecule:
    def __init__(self, name):
        self.name = name
        self.atoms = []
        self.charge_index = 0
        self.energy = 0

    def addAtom(self, symbol):
        self.atoms.append([symbol, 0])

    def addCharge(self, charge):
        self.atoms[self.charge_index][1] = charge
        self.charge_index += 1

def read(filename):
    molecules = []

    molecule = None
    
    in_atom_symbols = False
    in_atom_charges = False

    for line in file(filename):
        line = line.strip()

        if line.startswith("Structure Name:") or line.startswith("New Structure Name"):
            if(molecule):
                molecules.append(molecule)
            molecule = Molecule(line.split(":")[1].strip())

        #elif line.startswith("ATOM NAME SYMBOL"):
        elif line.startswith("ATOM NAME  TYPE"):
            in_atom_symbols = True

        elif line.startswith("ATOM    CHARGE"):
            in_atom_charges = True

        elif in_atom_symbols and line.startswith("OPTIMOL-LIST"):
            in_atom_symbols = False

        elif in_atom_charges and line.startswith("OPTIMOL-LIST"):
            in_atom_charges = False

        elif in_atom_symbols:
            line = line.split()
            for i in range(len(line) / 3):
                molecule.addAtom(line[i*3+2])

        elif in_atom_charges:
            line = line.split()
            for i in range(len(line) / 3):
                molecule.addCharge(line[i*3+2])

        elif line.startswith("Total ENERGY (Kcal)"):
            molecule.energy = line.split(' ')[-1]

    # get the last one
    if molecule:
        molecules.append(molecule)

    return molecules

def write(molecules, filename):
    f = open(filename, "w")
    f.write("<molecules>\n")

    for mol in molecules:
        f.write("  <molecule name=\"%s\" energy=\"%s\" atomCount=\"%i\">\n" % (mol.name, mol.energy, len(mol.atoms)))
        for id, atom in enumerate(mol.atoms):
            f.write("    <atom type=\"%s\" charge=\"%s\"/>\n" % (atom[0], atom[1]))
        f.write("  </molecule>\n")

    f.write("</molecules>\n")
    f.close()

if __name__ == '__main__':
    write(read("MMFF94_opti.log"), "mmff94.expected")

