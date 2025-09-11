proc _one_atom_xyz {mol seltext} {
    set sel [atomselect $mol $seltext]
    if {[$sel num] < 1} { $sel delete ; return {} }
    set xyz [lindex [$sel get {x y z}] 0]
    $sel delete
    return $xyz
}

proc _dist {a b} {
    if {$a eq {} || $b eq {}} { return 1e99 }
    lassign $a ax ay az
    lassign $b bx by bz
    set dx [expr {$ax-$bx}]
    set dy [expr {$ay-$by}]
    set dz [expr {$az-$bz}]
    return [expr {sqrt($dx*$dx+$dy*$dy+$dz*$dz)}]
}

proc heme_cys_binding {} {
    set mol [molinfo top]
    if {$mol < 0} { puts "No molecule loaded"; return }

    # Heme Fe atoms
    set fes [atomselect $mol "resname HEM HEME HEC and name FE"]
    if {[$fes num] == 0} { puts "No heme Fe atoms"; $fes delete; return }
    set fe_idx_list [$fes get index]
    set fe_resids   [$fes get resid]
    set fe_segnames [$fes get segname]
    set fe_chains   [$fes get chain]
    $fes delete

    # All cysteine sulfurs (SG or wildcard S*)
    set cands [atomselect $mol "resname CYS CYO and (name SG or name \"S*\")"]
    if {[$cands num] == 0} { puts "No cysteine sulfurs"; $cands delete; return }
    set cys_idx   [$cands get index]
    set cys_resid [$cands get resid]
    set cys_resnm [$cands get resname]
    $cands delete

    foreach fe_idx $fe_idx_list fe_res $fe_resids fe_seg $fe_segnames fe_chain $fe_chains {
        set fe_xyz [_one_atom_xyz $mol "index $fe_idx"]

        # Build safe base selection
        set base "resid $fe_res and resname HEM HEME HEC"
        if {$fe_seg ne ""} { append base " and segname $fe_seg" }
        if {$fe_chain ne ""} { append base " and chain $fe_chain" }

        set cab_xyz [_one_atom_xyz $mol "$base and name CAB C1B C4B CB"]
        set cad_xyz [_one_atom_xyz $mol "$base and name CAD C1D C4D CD"]

        if {$cab_xyz eq {} || $cad_xyz eq {}} {
            puts "Heme $fe_seg:$fe_res chain $fe_chain -> binding sites not found"
            continue
        }

        # Distances to cysteine sulfurs
        set dlist {}
        set n [llength $cys_idx]
        for {set i 0} {$i < $n} {incr i} {
            set idx [lindex $cys_idx $i]
            set s_xyz [_one_atom_xyz $mol "index $idx"]
            if {$s_xyz eq {}} { continue }
            set d [_dist $fe_xyz $s_xyz]
            lappend dlist [list $d $idx $s_xyz [lindex $cys_resid $i] [lindex $cys_resnm $i]]
        }
        set dlist [lsort -real -increasing -index 0 $dlist]
        set top2 [lrange $dlist 0 1]

        puts "Heme $fe_seg:$fe_res chain $fe_chain"
        foreach rec $top2 {
            lassign $rec dFeS s_idx s_xyz resid resnm
            set dSB [_dist $s_xyz $cab_xyz]
            set dSD [_dist $s_xyz $cad_xyz]
            puts "   $resnm:$resid atom $s_idx   S-CAB = $dSB A   S-CAD = $dSD A   (Fe-S = $dFeS A)"
        }
    }
}

# Run immediately
heme_cys_binding
