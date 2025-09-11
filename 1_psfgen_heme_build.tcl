# ===== Universal PSFGEN builder (order-safe, gap-safe) =====
# Builds all categories, prevents “long diagonal bonds”, writes & loads:
#   <top_molecule_name>_build.psf/.pdb
# No patches are applied here; paste them after it loads if needed.

# --- Topology files (edit paths if needed) ---
set topo_prot {C:/Program Files/VMD/plugins/noarch/tcl/readcharmmtop1.2/top_all36_prot.rtf}
set topo_heme {C:/Users/xiong/Desktop/labstuff/toppar_all36_prot_heme.rtf}
set topo_wion {C:/Program Files/VMD/plugins/noarch/tcl/readcharmmtop1.2/toppar_water_ions_namd.str}

# --- Need a molecule; remember original top mol to delete later ---
set origMol [molinfo top]
if {$origMol < 0} { puts "ERROR: No molecule loaded."; return }

# --- Output name from top molecule name, with _build suffix ---
proc __sanitize_name {s} {
    set s [file tail $s]
    regsub -all {\.(pdb|psf|gro|mol2)$} $s "" s
    regsub -all {[^A-Za-z0-9._-]} $s "_" s
    if {$s eq ""} { set s "structure" }
    return $s
}
set outbase "[__sanitize_name [molinfo $origMol get name]]_build"

# --- Small helpers ---
proc __nonempty {lst} {
    set out {}
    foreach x $lst { if {[string trim $x] ne ""} { lappend out $x } }
    return $out
}
proc __sel_count {mol seltext} {
    set s [atomselect $mol $seltext]
    set n [$s num]
    $s delete
    return $n
}
proc __write_sel_pdb {mol seltext path} {
    set s [atomselect $mol $seltext]
    set n [$s num]
    if {$n > 0} { $s writepdb $path }
    $s delete
    return $n
}

# ---- Order-preserving split into “occurrence groups” (handles duplicate PDB resids) ----
# Returns a list of groups; each group is an ordered list of VMD 'residue' ids.
proc __occurrence_groups {mol seltext} {
    set s [atomselect $mol $seltext]
    # tuples: (residue, pdbresid, atom index) — we sort by atom index to keep original file order
    set tuples [$s get {residue resid index}]
    $s delete
    set tuples [lsort -integer -index 2 $tuples]

    array unset seenOcc
    array unset seenResidInGroup
    set groups {}

    foreach tup $tuples {
        lassign $tup rid pdbres idx
        if {![info exists seenOcc($pdbres)]} { set seenOcc($pdbres) 0 }
        set k $seenOcc($pdbres)
        incr seenOcc($pdbres)

        while {[llength $groups] <= $k} { lappend groups {} }
        if {![info exists seenResidInGroup($k,$rid)]} {
            lset groups $k [concat [lindex $groups $k] [list $rid]]
            set seenResidInGroup($k,$rid) 1
        }
    }
    return $groups
}

# ---- Split at peptide C..N discontinuities (prevents psfgen spanning gaps) ----
# Input: ordered list of VMD 'residue' ids; Output: list of fragments (each is ordered rid list)
proc __split_by_backbone_breaks {mol sel_prefix ordered_rids} {
    set MAX_CN_DIST 1.9   ;# Å; relax to 2.1 if your PDB is coarse
    set frags {}
    set cur {}

    proc __firstIndex {mol sel} {
        set s [atomselect $mol $sel]
        set idx -1
        if {[$s num] > 0} { set idx [lindex [$s get index] 0] }
        $s delete
        return $idx
    }

    foreach rid $ordered_rids {
        if {[llength $cur] == 0} { lappend cur $rid ; continue }
        set prev [lindex $cur end]
        set Cidx [__firstIndex $mol "$sel_prefix and residue $prev and name C"]
        set Nidx [__firstIndex $mol "$sel_prefix and residue $rid  and name N"]
        set break 0
        if {$Cidx < 0 || $Nidx < 0} {
            set break 1
        } else {
            set sc [atomselect $mol "index $Cidx"]
            set sn [atomselect $mol "index $Nidx"]
            set rc [$sc get {x y z}]
            set rn [$sn get {x y z}]
            $sc delete; $sn delete
            set dx [expr {[lindex $rc 0]-[lindex $rn 0]}]
            set dy [expr {[lindex $rc 1]-[lindex $rn 1]}]
            set dz [expr {[lindex $rc 2]-[lindex $rn 2]}]
            set d  [expr {sqrt($dx*$dx + $dy*$dy + $dz*$dz)}]
            if {$d > $MAX_CN_DIST} { set break 1 }
        }
        if {$break} {
            if {[llength $cur] > 0} { lappend frags $cur }
            set cur [list $rid]
        } else {
            lappend cur $rid
        }
    }
    if {[llength $cur] > 0} { lappend frags $cur }
    return $frags
}

# ---- Build a category (segname-first, then chain; auto-split; protein-safe) ----
# code = P (protein), O (heme), W (water), I (ion), L (other)
# firstlast = "protein" -> NTER/CTER; "none" -> no termini patches
proc __build_category {mol basesel code firstlast tmpdir} {
    # decide grouping key: segname if present; else chain; else single group
    set sel_all [atomselect $mol $basesel]
    set all_segnames [__nonempty [lsort -unique [$sel_all get segname]]]
    set all_chains   [__nonempty [lsort -unique [$sel_all get chain]]]
    $sel_all delete

    set groupspecs {}
    if {[llength $all_segnames] > 0} {
        foreach S $all_segnames { lappend groupspecs [list segname $S] }
    } elseif {[llength $all_chains] > 0} {
        foreach C $all_chains { lappend groupspecs [list chain $C] }
    } else {
        lappend groupspecs [list all X]
    }

    foreach spec $groupspecs {
        lassign $spec mode val
        if {$mode eq "segname"} {
            set base [string toupper [string range $val 0 1]]
            if {$base eq ""} { set base "X" }
            set subsel "$basesel and segname $val"
        } elseif {$mode eq "chain"} {
            set base [string toupper [string range $val 0 1]]
            if {$base eq ""} { set base "X" }
            set subsel "$basesel and chain $val"
        } else {
            set base "X"
            set subsel $basesel
        }

        # 1) split by duplicate PDB resids, preserving original order
        set groups [__occurrence_groups $mol $subsel]

        if {$firstlast eq "protein"} {
            # 2) for each group, further split at peptide gaps
            set gi 0
            foreach g $groups {
                incr gi
                set frags [__split_by_backbone_breaks $mol $subsel $g]
                set fj 0
                foreach frag $frags {
                    incr fj
                    set tag [string range "${base}$code$gi$fj" 0 3]   ;# e.g., AP11, AP12, ...
                    set ridlist [join $frag " "]
                    set pdb [file join $tmpdir "${tag}.pdb"]
                    if {[__write_sel_pdb $mol "$subsel and residue $ridlist" $pdb] == 0} { continue }
                    segment $tag { first NTER ; last CTER ; pdb $pdb }
                    coordpdb $pdb $tag
                    puts "  Built segment $tag   ($subsel ; group $gi frag $fj ; residues [llength $frag])"
                }
            }
        } else {
            # non-proteins: build direct (no peptide-gap logic)
            set gi 0
            foreach g $groups {
                incr gi
                set tag [string range "${base}$code$gi" 0 3]          ;# e.g., AO1, AW2, ...
                set ridlist [join $g " "]
                set pdb [file join $tmpdir "${tag}.pdb"]
                if {[__write_sel_pdb $mol "$subsel and residue $ridlist" $pdb] == 0} { continue }
                segment $tag { first none ; last none ; pdb $pdb }
                coordpdb $pdb $tag
                puts "  Built segment $tag   ($subsel ; group $gi ; residues [llength $g])"
            }
        }
    }
}

# --- Hard stop if any HEC present ---
set hecSel [atomselect $origMol "resname HEC"]
set hecN [$hecSel num]
if {$hecN > 0} {
    set info [$hecSel get {resid chain segname}]
    $hecSel delete
    puts "ERROR: Found $hecN residue(s) with resname HEC. Rename/handle HEC, then re-run."
    foreach t $info { lassign $t r c s ; puts "  HEC resid $r chain $c segname $s" }
    return
}
$hecSel delete

# --- Define categories (everything gets built) ---
set heme_names {HEM HEME HEO HEMO PHEM FHEM}
set sel_prot   "protein"
set sel_heme   "resname [join $heme_names " "]"
set sel_water  "water"
# If your VMD lacks 'ions' macro, replace with explicit list, e.g.:
# set sel_ions "resname NA CL K MG CA ZN LI RB CS BA CD"
set sel_ions   "ions"
set sel_other  "not ($sel_prot) and not ($sel_water) and not ($sel_ions) and not ($sel_heme)"

# --- psfgen start + topos ---
package require psfgen
if {[llength [info commands psfcontext]]} { psfcontext reset }
resetpsf
topology $topo_prot
topology $topo_heme
topology $topo_wion

# Useful aliases
pdbalias atom ILE CD1 CD
pdbalias residue HIE HSE
pdbalias residue HID HSD
pdbalias residue HIS HSD
pdbalias residue CYX CYS
pdbalias residue CYM CYS
pdbalias residue CYO CYS

# Temp dir for intermediate PDBs
set tmpdir __tmp_psfgen_build
file mkdir $tmpdir

# ===================== BUILD ALL CATEGORIES =====================
if {[__sel_count $origMol $sel_prot]  > 0} { __build_category $origMol $sel_prot  P protein $tmpdir } else { puts "NOTE: no protein detected." }
if {[__sel_count $origMol $sel_heme]  > 0} { __build_category $origMol $sel_heme  O none    $tmpdir } else { puts "NOTE: no heme-like residues detected (HEM/HEME/…; excluding HEC)." }
if {[__sel_count $origMol $sel_water] > 0} { __build_category $origMol $sel_water W none    $tmpdir } else { puts "NOTE: no water detected." }
if {[__sel_count $origMol $sel_ions]  > 0} { __build_category $origMol $sel_ions  I none    $tmpdir } else { puts "NOTE: no ions detected." }
if {[__sel_count $origMol $sel_other] > 0} { __build_category $origMol $sel_other L none    $tmpdir } else { puts "NOTE: no other hetero-ligands detected." }

# Regenerate angles/dihedrals after proteins
regenerate angles dihedrals

# ===================== Write baseline (no patches yet) =====================
# Tip: you can comment out 'guesscoord' here and run it AFTER you paste heme patches to reduce warnings.
guesscoord
writepsf ${outbase}.psf
writepdb ${outbase}.pdb
puts "Baseline written: ${outbase}.psf  ${outbase}.pdb"

# ===================== Replace original molecule with the new build =====================
mol delete $origMol
mol new ${outbase}.psf type psf waitfor all
mol addfile ${outbase}.pdb type pdb waitfor all
puts "Replaced original molecule with: ${outbase}.psf/.pdb"

# Cleanup temp files
foreach f [glob -nocomplain [file join __tmp_psfgen_build *.pdb]] { file delete -force $f }
file delete -force __tmp_psfgen_build

puts "\nNext: paste patch lines (e.g., 'patch PHEM AP1:501 AO1:38' or 'patch PHEM A1:501 AO1:38'), then run:"
puts "  guesscoord"
puts "  writepsf ${outbase}.psf"
puts "  writepdb ${outbase}.pdb"
