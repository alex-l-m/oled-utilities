from rdkit.Chem.rdchem import KekulizeException
from rdkit.Chem.rdChemReactions import ReactionFromSmarts
from rdkit.Chem.AllChem import ConstrainedEmbed
from rdkit.Chem.rdmolops import SanitizeMol, AddHs, RemoveHs
from tblite.ase import TBLite
from ase.constraints import FixAtoms
from annotate_rdkit_with_ase import optimize_geometry

def substitute(mol, substitution_string, removal_string):

    mol = RemoveHs(mol)

    for atom in mol.GetAtoms():
        atom.SetBoolProp("original", True)

    core_reaction = ReactionFromSmarts(removal_string)
    substitution_reaction = ReactionFromSmarts(substitution_string)
    core_results = ((mol,),)
    subs_results = ((mol,),)
    new_core_results = core_reaction.RunReactants([mol])
    new_subs_results = substitution_reaction.RunReactants([mol])
    while len(new_subs_results) > 0:
        SanitizeMol(new_core_results[0][0])
        SanitizeMol(new_subs_results[0][0])
        core_results = new_core_results
        subs_results = new_subs_results
        new_core_results = core_reaction.RunReactants([core_results[0][0]])
        new_subs_results = substitution_reaction.RunReactants([subs_results[0][0]])

    subs_mol = subs_results[0][0]
    core = core_results[0][0]

    SanitizeMol(subs_mol)
    SanitizeMol(core)

    subs_mol = AddHs(subs_mol)
    # Using settings from here:
    # https://sourceforge.net/p/rdkit/mailman/rdkit-discuss/thread/C761AFBF8DEB604DB8D72CB1B301A1EB21C651DB%40MBX07.ad.oak.ox.ac.uk/#msg32082674
    # Prevents some cases of "ValueError: Could not embed molecule."
    ConstrainedEmbed(subs_mol, core, ignoreSmoothingFailures=True, useRandomCoords=True, useTethers=True)

    # Set the "original" property for all atoms
    # Has to be done after adding hydrogens, or they won't get it
    for atom in subs_mol.GetAtoms():
        if not atom.HasProp("original"):
            atom.SetBoolProp("original", False)
    is_original = [atom.GetBoolProp("original") for atom in subs_mol.GetAtoms()]
    ase_constraint = FixAtoms(mask = is_original)
    optimize_geometry(TBLite(), subs_mol,
                      conformation_index = 0,
                      constraints = [ase_constraint])

    return subs_mol
