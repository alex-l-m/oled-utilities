Geometry generation for iridium-based phosphorescent emitter molecules.

RDKit's EmbedMolecule is not parametrized for metals, and will not generate the correct octahedral geometry around the iridium.
This package provides a replacement function, "octahedral\_embed", which generates a geometry constrained to have the octahedral geometry.
To generate the geometry, the input molecule is matched to a template "core", containing the iridium and surrounding atoms, extracted from a crystal structure.

Furthermore, octahedral\_embed provides control over the isomer, through the "isomer" argument, set to "fac" or "mer".
Tridentate carbenes are handled separately by setting it to "tridentate".
