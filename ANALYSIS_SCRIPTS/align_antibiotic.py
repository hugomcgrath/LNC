import MDAnalysis as mda
from MDAnalysis.analysis.align import alignto
import shared_data as sd


alignment_keyword = f"not segid 9 and {sd.HEAVY_ATOMS_SELECTION_KEYWORD} and {sd.NUCLEIC_PRUNED_SELECTION_KEYWORD}"
system_path = sd.SYSTEMS["LNC"]
reference_pdb_path = f"{system_path}/common_atoms.pdb"
reference = mda.Universe(reference_pdb_path)
selection_reference = reference.select_atoms("all")
pdb_in = f"{system_path}/system.pdb"
pdb_out = f"{system_path}/antibiotic_aligned.pdb"
mobile = mda.Universe(pdb_in)
selection_mobile = mobile.select_atoms("all")

print(f"aligning {pdb_in}")
alignto(selection_mobile, selection_reference, select=alignment_keyword)
print(f"writing to {pdb_out}")
selection_mobile.write(pdb_out)
