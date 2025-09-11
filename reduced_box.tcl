# ==== Reduce top molecule ====
# Keeps: protein + surface (resname AUM) + waters inside bounding box
# Writes: *_reduced.pdb (first frame) and *_reduced.dcd (all frames)

# ---- Box bounds ----
set xmin -30.0
set xmax  30.0
set ymin -40.0
set ymax  40.0
set zmin -20.0
set zmax  30.0

# ---- Setup ----
set mid [molinfo top]
if {$mid < 0} { puts "ERROR: No top molecule loaded."; return }
set nframes [molinfo $mid get numframes]
if {$nframes < 1} { puts "ERROR: No trajectory loaded."; return }

# ---- Output names ----
set rawname [molinfo $mid get name]
if {$rawname eq ""} { set rawname "molid$mid" }
set rawname [file tail $rawname]
set ext [string tolower [file extension $rawname]]
if {$ext in {".psf" ".pdb" ".dcd"}} { set rawname [file rootname $rawname] }
set safename [string map {" " "_" "/" "_" "\\" "_" ":" "_" ";" "_"} $rawname]
set outbase "${safename}_reduced"
set out_pdb "${outbase}.pdb"
set out_dcd "${outbase}.dcd"

puts "Reducing molecule: [molinfo $mid get name]"
puts "Box: x $xmin..$xmax  y $ymin..$ymax  z $zmin..$zmax"
puts "Frames: $nframes"

# ---- Build selection ----
# Always keep protein and AUM surface atoms
set keepCore [atomselect $mid "protein or resname AUM" frame 0]
set coreIdx [$keepCore get index]
$keepCore delete

# Waters inside the box (only evaluated on frame 0)
set keepWat [atomselect $mid "water and x>$xmin and x<$xmax and y>$ymin and y<$ymax and z>$zmin and z<$zmax" frame 0]
set watIdx [$keepWat get index]
$keepWat delete

# Combine
set keepIdx [lsort -unique [concat $coreIdx $watIdx]]
if {[llength $keepIdx] == 0} {
    puts "ERROR: Empty keep selection."; return
}

puts "Keeping atoms: [llength $keepIdx]"

# ---- Write reduced PDB for frame 0 ----
set refSel [atomselect $mid "index [join $keepIdx " "]" frame 0]
$refSel writepdb $out_pdb
$refSel delete
puts "Wrote reduced PDB: $out_pdb"

# ---- Write reduced DCD for all frames ----
set selAll [atomselect $mid "index [join $keepIdx " "]"]
animate write dcd $out_dcd sel $selAll waitfor all
$selAll delete
puts "Wrote reduced DCD: $out_dcd"

puts "=== Done: $outbase ==="
