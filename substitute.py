from rdkit.Chem.rdChemReactions import ReactionFromSmarts
from rdkit.Chem.AllChem import ConstrainedEmbed
from rdkit.Chem.rdmolops import SanitizeMol

def substitute(mol, substitution_string, removal_string):
    for atom in mol.GetAtoms():
        atom.SetBoolProp("original", True)

    core_reaction = ReactionFromSmarts(removal_string)
    substitution_reaction = ReactionFromSmarts(substitution_string)
    core_results = ((mol,))
    subs_results = ((mol,))
    new_core_results = core_reaction.RunReactants([mol])
    new_subs_results = substitution_reaction.RunReactants([mol])
    while len(new_subs_results) > 0:
        core_results = new_core_results
        subs_results = new_subs_results
        new_core_results = core_reaction.RunReactants([core_results[0][0]])
        new_subs_results = substitution_reaction.RunReactants([subs_results[0][0]])

    subs_mol = subs_results[0][0]
    core = core_results[0][0]
    for atom in subs_mol.GetAtoms():
        if not atom.HasProp("original"):
            atom.SetBoolProp("original", False)

    SanitizeMol(subs_mol)
    ConstrainedEmbed(subs_mol, core)

    return subs_mol
