# Select HEC residues and rename to HEME
set hec [atomselect top "resname HEC"]
$hec set resname HEME
$hec delete

# Get the name of the current top molecule
set name [molinfo top get name]

# Write out the modified PDB with the molecule name included
set allsel [atomselect top all]
$allsel writepdb "${name}_HEME.pdb"
$allsel delete

# Reload the modified structure
mol new "${name}_HEME.pdb"
