# ===== Monomer keep + half SOD + NO water under gold + 3 Å solvate + final PSF/PDB =====
# Outputs:
#   <basename>_monomer.psf/.pdb       (trimmed, SOD-thinned)
#   <basename>_monomer_solv.psf/.pdb  (solvated box around protein; no water below gold)

# --- USER PATHS (Windows paths in braces!) ---
set psfPath {C:\Users\xiong\Desktop\labstuff\OmcS_Au\omcs_au_ox_wb_ions.psf}
set pdbPath {C:\Users\xiong\Desktop\labstuff\OmcS_Au\omcs_au_ox_wb_ions.pdb}

# --- USER OPTIONS ---
# Chains to KEEP (space-separated). Fallback to segname with the same tokens if chain is empty.
set KEEP_CHAINS "B I Q O S U Z D"
# Which chain holds the SOD ions to thin by 50% (fallback to segname).
set SOD_CHAIN "I"
# Water padding (Å) around protein (X/Y/Z). Bottom Z is clamped above gold.
set PAD 3.0
# Gold selection (edit if your naming differs)
set GOLD_SEL "resname AUM or element AU"
# Include these cofactors with protein when measuring min/max for solvation box:
set PROT_WITH_COF "protein or resname HEM HEC HEME HEO"

# --- helpers ---
proc _basename_noext {name} {
    set base [file rootname [file tail $name]]
    regsub -all {\s+} $base "_" base
    return $base
}
proc _rand_pick_half {L} {
    # Fisher–Yates via random insertion; returns floor(n/2) elements
    set A {}
    foreach x $L {
        set r [expr {int(rand()*([llength $A]+1))}]
        set A [linsert $A $r $x]
    }
    set n [llength $A]
    return [lrange $A 0 [expr {($n/2)-1}]]
}
proc _triplets {mol seltext} {
    set s [atomselect $mol $seltext]
    set t {}
    if {[$s num] > 0} { set t [$s get {segname resid name}] }
    $s delete
    return $t
}
proc _indices {mol seltext} {
    set s [atomselect $mol $seltext]
    set ids {}
    if {[$s num] > 0} { set ids [$s get index] }
    $s delete
    return $ids
}
proc _keep_by_chains_or_segnames {mol keeplist} {
    set s [atomselect $mol "chain $keeplist"]
    if {[$s num] > 0} { set idx [$s get index]; $s delete; return $idx }
    $s delete
    set s [atomselect $mol "segname $keeplist"]
    set idx {}
    if {[$s num] > 0} { set idx [$s get index] }
    $s delete
    return $idx
}

# --- sanity ---
if {![file exists $psfPath]} { error "PSF not found: $psfPath" }
if {![file exists $pdbPath]} { error "PDB not found: $pdbPath" }

# --- 0) Load source into a fresh molecule for robust selections ---
set mid [mol new $psfPath type psf waitfor all]
mol addfile $pdbPath type pdb waitfor all top $mid

# --- 1) KEEP only chosen chains (fallback segname) ---
set baseKeepIdx [_keep_by_chains_or_segnames $mid $KEEP_CHAINS]
if {[llength $baseKeepIdx] == 0} {
    mol delete $mid
    error "No atoms found for {$KEEP_CHAINS} in chain or segname."
}

# --- 2) Randomly remove HALF of SOD residues in SOD_CHAIN (whole residues; fallback to segname) ---
set sodSel [atomselect $mid "resname SOD and chain $SOD_CHAIN"]
if {[$sodSel num] == 0} {
    $sodSel delete
    set sodSel [atomselect $mid "resname SOD and segname $SOD_CHAIN"]
}
set rmSODatomIdx {}
if {[$sodSel num] > 0} {
    set resid_list  [$sodSel get resid]
    set chain_list  [$sodSel get chain]
    set seg_list    [$sodSel get segname]
    $sodSel delete

    # unique residue keys (seg|chain|resid)
    set keys {}
    for {set i 0} {$i < [llength $resid_list]} {incr i} {
        lappend keys "[lindex $seg_list $i]|[lindex $chain_list $i]|[lindex $resid_list $i]"
    }
    set keys   [lsort -unique $keys]
    set rmkeys [_rand_pick_half $keys]

    foreach k $rmkeys {
        lassign [split $k "|"] seg ch rs
        set seltext "resname SOD and resid $rs"
        if {$seg ne ""} { append seltext " and segname $seg" }
        if {$ch  ne ""} { append seltext " and chain $ch" }
        set s [atomselect $mid $seltext]
        if {[$s num] > 0} { foreach i [$s get index] { lappend rmSODatomIdx $i } }
        $s delete
    }
}

# --- 3) Final KEEP = base keep minus chosen SOD atoms ---
array set rmset {}
foreach i $rmSODatomIdx { set rmset($i) 1 }
set finalKeepIdx {}
foreach i $baseKeepIdx {
    if {![info exists rmset($i)]} { lappend finalKeepIdx $i }
}
if {[llength $finalKeepIdx] == 0} {
    mol delete $mid
    error "Final keep set is empty after SOD thinning."
}

# --- 4) Build deletion list (segname/resid/name) WITHOUT long 'index' selections ---
set allSel   [atomselect $mid all]
set all_idx  [$allSel get index]
set all_seg  [$allSel get segname]
set all_res  [$allSel get resid]
set all_nam  [$allSel get name]
$allSel delete
array set keepmap {}
foreach i $finalKeepIdx { set keepmap($i) 1 }

set kill_triplets {}
for {set j 0} {$j < [llength $all_idx]} {incr j} {
    set ai [lindex $all_idx $j]
    if {[info exists keepmap($ai)]} { continue }
    set seg [lindex $all_seg $j]
    set rs  [lindex $all_res $j]
    set nm  [lindex $all_nam $j]
    if {$seg eq "" || $nm eq ""} {
        mol delete $mid
        error "Atom index $ai missing segname/name; ensure your PSF provides segment IDs."
    }
    lappend kill_triplets [list $seg $rs $nm]
}

# --- 5) Prune in psfgen to write the trimmed pair (<stem>_monomer) ---
package require psfgen
resetpsf
if {[catch {readpsf $psfPath} msg]}  { mol delete $mid ; error "psfgen readpsf failed: $msg" }
if {[catch {coordpdb $pdbPath} msg]} { mol delete $mid ; error "psfgen coordpdb failed: $msg" }

# also delete **all water** so nothing remains under gold before we solvate
set water_trip [_triplets $mid "water"]

set removed 0
foreach trip $kill_triplets {
    lassign $trip seg rs nm
    if {![catch {delatom $seg $rs $nm}]} { incr removed }
}
foreach trip $water_trip {
    lassign $trip seg rs nm
    if {![catch {delatom $seg $rs $nm}]} { incr removed }
}

set stem "[_basename_noext $psfPath]_monomer"
writepsf ${stem}.psf
writepdb ${stem}.pdb
puts ">> Trimmed & dried: ${stem}.psf / ${stem}.pdb   (removed $removed atoms)"

# --- 6) Measure protein box and gold top from the trimmed/dry system ---
set mid2 [mol new ${stem}.psf type psf waitfor all]
mol addfile ${stem}.pdb type pdb waitfor all top $mid2

# protein + cofactors for bounding box
set prot [atomselect $mid2 $PROT_WITH_COF]
if {[$prot num] == 0} {
    $prot delete
    set prot [atomselect $mid2 "not (water or ions or ($GOLD_SEL))"]
}
lassign [measure minmax $prot] bbmin bbmax
$prot delete
set xmin [lindex $bbmin 0]; set ymin [lindex $bbmin 1]; set zmin [lindex $bbmin 2]
set xmax [lindex $bbmax 0]; set ymax [lindex $bbmax 1]; set zmax [lindex $bbmax 2]

# gold top (max Z)
set gold [atomselect $mid2 $GOLD_SEL]
set zgoldtop ""
if {[$gold num] > 0} {
    set zlist [$gold get z]
    # tcl::mathfunc::max may not be present on some builds; compute manually
    set zgoldtop [lindex $zlist 0]
    foreach z $zlist { if {$z > $zgoldtop} { set zgoldtop $z } }
}
$gold delete
mol delete $mid2

# build custom solvation box: X/Y +/- PAD; Z top = zmax+PAD; Z bottom = max(zmin-PAD, zgoldtop+0.5)
set out_xmin [expr {$xmin - $PAD}]
set out_ymin [expr {$ymin - $PAD}]
set out_xmax [expr {$xmax + $PAD}]
set out_ymax [expr {$ymax + $PAD}]
set out_zmax [expr {$zmax + $PAD}]
if {$zgoldtop eq ""} {
    set out_zmin [expr {$zmin - $PAD}]
} else {
    set out_zmin [expr {max($zmin - $PAD, $zgoldtop + 0.5)}]
}

# --- 7) Solvate using the custom box, then write the final solvated pair ---
package require solvate
solvate -psf ${stem}.psf -pdb ${stem}.pdb \
        -o   ${stem}_solv \
        -minmax "{{$out_xmin $out_ymin $out_zmin} {$out_xmax $out_ymax $out_zmax}}"

puts ">> Solvated (3 Å padding; clamped above gold): ${stem}_solv.psf / ${stem}_solv.pdb"

# --- cleanup ---
mol delete $mid
