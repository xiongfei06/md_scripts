Important scripts for going from raw pdb --> pdb/psf --> ...
Scripts listed in order of usage in the workflow:

0. HEC_to_HEME_pdb.tcl: takes all HEC named groups and renames to HEME. creates a new file called "${name}_HEME.pdb"

1-4 psfgen sequence:
  Once complete, these scripts will take in a pdb file with renamed HEC groups and then return a new pdb/psf pair with defined ligand connectivity
1. 1_psfgen_heme_build.tcl: loads top molecule, uses 3 topology files using pathing, creates pdb/psf pair with defined segments and chains
2. 2_heme_histidine_patch.tcl: calculates closest histidine to each heme and creates {$name}_patch.txt document with patch information, defining patches for the two closest histidines.
3. 3_HIS_HSO.tcl: reads {$name}_hispatch.txt document and converts corresponding residues to HSO residues

need to make another script that finds the closest cysteines to hemes and another script that appends new patch information to {$name}_patch.txt
4. 4_heme_cysteine_patch.tcl
5. 5_CYS_CYO.tcl
