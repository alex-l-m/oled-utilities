import os.path
from rdkit import Chem
from rdkit.Chem.rdchem import RWMol
from rdkit.Chem.rdmolfiles import MolFromMol2File, MolFromSmarts
from rdkit.Chem.AllChem import ConstrainedEmbed
from rdkit.Chem.rdChemReactions import ReactionFromSmarts
from rdkit.Chem.rdmolops import RemoveStereochemistry

def make_bonds_dative(mol, target_elem = "Ir"):
    editable_mol = RWMol(mol)

    # If you don't make a list, it loops infinitely over the bonds it's creating
    for bond in list(editable_mol.GetBonds()):
        iridium = None
        nitrogen = None
        carbene = None
        if bond.GetBeginAtom().GetSymbol() == target_elem and \
                bond.GetEndAtom().GetSymbol() in ["N", "P"] and \
                bond.GetEndAtom().GetFormalCharge() == 1:
            iridium = bond.GetBeginAtom()
            nitrogen = bond.GetEndAtom()
            start_idx = bond.GetEndAtomIdx()
            end_idx = bond.GetBeginAtomIdx()
        elif bond.GetEndAtom().GetSymbol() == target_elem and \
                bond.GetBeginAtom().GetSymbol() in ["N", "P"] and \
                bond.GetBeginAtom().GetFormalCharge() == 1:
            iridium = bond.GetEndAtom()
            nitrogen = bond.GetBeginAtom()
            start_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
        if bond.GetBeginAtom().GetSymbol() == target_elem and \
                bond.GetEndAtom().GetSymbol() == "C" and \
                bond.GetEndAtom().GetTotalValence() == 3:
            iridium = bond.GetBeginAtom()
            carbene = bond.GetEndAtom()
            start_idx = bond.GetEndAtomIdx()
            end_idx = bond.GetBeginAtomIdx()
        elif bond.GetEndAtom().GetSymbol() == target_elem and \
                bond.GetBeginAtom().GetSymbol() == "C" and \
                bond.GetBeginAtom().GetTotalValence() == 3:
            iridium = bond.GetEndAtom()
            carbene = bond.GetBeginAtom()
            start_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()

        if nitrogen is not None:
            # Replace N+ - Ir with N -> Ir
            nitrogen.SetFormalCharge(0)

        if iridium is not None and (nitrogen is not None or carbene is not None):
            editable_mol.RemoveBond(start_idx, end_idx)
            editable_mol.AddBond(start_idx, end_idx, Chem.rdchem.BondType.DATIVE)

    outmol = editable_mol.GetMol()
    Chem.SanitizeMol(outmol)

    return outmol

fac = make_bonds_dative(MolFromMol2File(os.path.join(__path__[0], "OHUZEW.mol2")))
mer = make_bonds_dative(MolFromMol2File(os.path.join(__path__[0], "OHUZIA.mol2")))

template = MolFromSmarts("[Ir]1~n:[*]~[*]:c~1")

# Extract skeletons of a molecule based on a template, keeping coordinates
# Multiple skeletons because I don't know how to do wildcards
def skeleton(template, mol):
    template_matches = mol.GetSubstructMatches(template)
    matching_indices = set(sum(template_matches, ()))
    editable_mol = RWMol(mol)
    editable_mol.BeginBatchEdit()
    for atom in editable_mol.GetAtoms():
        # Remove all atoms except the matching ones
        if atom.GetIdx() not in matching_indices:
            editable_mol.RemoveAtom(atom.GetIdx())
        # All the atoms should have formal charge 0, not sure why they don't
        atom.SetFormalCharge(0)
    editable_mol.CommitBatchEdit()
    skeleton_mol = editable_mol.GetMol()
    return skeleton_mol

fac_skeleton = skeleton(template, fac)
mer_skeleton = skeleton(template, mer)

template = MolFromSmarts("[Ir]1~n:[*]~[*]:c~1")

def run_three_times(mol, reaction):
    for i in range(3):
        mol = reaction.RunReactants([mol])[0][0]
    return mol
reactions = [
    ReactionFromSmarts("[Ir:1]1~[n:2]:[n:3]~[c:4]:[c:5]~1>>[Ir:1]1~[n:2]:[c:3]~[n:4]:[c:5]~1"),
    ReactionFromSmarts("[Ir:1]1~[n:2]:[n:3]~[c:4]:[c:5]~1>>[Ir:1]1~[n:2]:[n:3]~[n:4]:[c:5]~1"),
    ReactionFromSmarts("[Ir:1]1~[n:2]:[n:3]~[c:4]:[c:5]~1>>[Ir:1]1~[n:2]:[c:3]~[c:4]:[c:5]~1")
]
fac_skeletons = [fac_skeleton] + [run_three_times(fac_skeleton, reaction) for reaction in reactions]
mer_skeletons = [mer_skeleton] + [run_three_times(mer_skeleton, reaction) for reaction in reactions]

def octahedral_embed(mol, isomer):
    # Needed for some of the mol2 files I got from CSD
    # Will not be able to embed with stereochemistry
    RemoveStereochemistry(mol)
    if isomer == "fac":
        skeletons = fac_skeletons
    elif isomer == "mer":
        skeletons = mer_skeletons
    else:
        raise ValueError(f"Isomer should be \"mer\" or \"fac\", given {isomer}")
    finished = False
    for skeleton in skeletons:
        if len(mol.GetSubstructMatch(skeleton)) > 0:
            ConstrainedEmbed(mol, skeleton)
            finished = True
    if not finished:
        raise ValueError("Doesn't match templates")
