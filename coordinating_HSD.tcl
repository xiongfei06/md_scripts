# Collect all HSD NE2 atoms once
set hsds [atomselect top "resname HSD and name NE2"]
set h_xyzs   [$hsds get {x y z}]
set h_resids [$hsds get resid]
set n_hsd [llength $h_resids]

if {$n_hsd < 1} {
    puts "No HSD NE2 atoms found."
} else {
    # Loop over all HEME Fe atoms
    set fes [atomselect top "resname HEME and name FE"]
    set fe_indices [$fes get index]
    $fes delete

    foreach fe_idx $fe_indices {
        set fe [atomselect top "index $fe_idx"]
        set heme_resid [lindex [$fe get resid] 0]
        set fe_xyz [lindex [$fe get {x y z}] 0]

        # Build list of {distance hsd_resid}
        set dist_list {}
        for {set i 0} {$i < $n_hsd} {incr i} {
            set h_xyz [lindex $h_xyzs $i]
            set d [veclength [vecsub $fe_xyz $h_xyz]]
            lappend dist_list [list $d [lindex $h_resids $i]]
        }

        # Sort by distance and take the first two unique resids
        set sorted_list [lsort -real -index 0 $dist_list]
        set picked {}
        foreach entry $sorted_list {
            set d [lindex $entry 0]
            set r [lindex $entry 1]
            if {[lsearch -exact [lmap x $picked {lindex $x 0}] $r] < 0} {
                lappend picked [list $r $d]
            }
            if {[llength $picked] == 2} {break}
        }

        # Print line: "HEME <resid>: HSD<resid1>(<dist1> Å) HSD<resid2>(<dist2> Å)"
        set out "HEME $heme_resid:"
        foreach p $picked {
            set r [lindex $p 0]
            set d [format "%.2f" [lindex $p 1]]
            append out " HSD$r($d Å)"
        }
        puts $out

        $fe delete
    }
}
$hsds delete
