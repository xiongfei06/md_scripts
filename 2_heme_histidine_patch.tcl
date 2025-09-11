# === HEME–HIS scan → patch-file writer (self-contained) ===
# Paste into VMD Tk Console; it runs immediately.

# ---- helpers ----
proc __xyz1 {mol seltext} {
    set s [atomselect $mol $seltext]
    if {[$s num] == 0} { $s delete ; return {} }
    set xyz [lindex [$s get {x y z}] 0]
    $s delete
    return $xyz
}
proc __dist {a b} {
    if {$a eq {} || $b eq {}} { return 1e99 }
    lassign $a ax ay az
    lassign $b bx by bz
    set dx [expr {$ax-$bx}] ; set dy [expr {$ay-$by}] ; set dz [expr {$az-$bz}]
    return [expr {sqrt($dx*$dx+$dy*$dy+$dz*$dz)}]
}
proc __seg {segid segname} {
    return [expr {$segid ne "" ? $segid : $segname}]
}
proc __fmt_sr {seg resid} { return "${seg}:$resid" }

# ---- main ----
set mol [molinfo top]
if {$mol < 0} { puts "ERROR: No molecule loaded."; return }

# Derive output filename base from top molecule name
set topname_raw [molinfo $mol get name]
set topname     [file rootname $topname_raw]
if {$topname eq ""} { set topname "top" }
set patch_file "${topname}_patch.txt"   ;# <-- output: exact patch commands

# Collect all HEME/HEC FE atoms
set fe_sel [atomselect $mol "resname HEME HEC and name FE"]
if {[$fe_sel num] == 0} { $fe_sel delete ; puts "No HEME/HEC FE atoms found."; return }
set fe_info [$fe_sel get {resname resid segid segname chain index}]
$fe_sel delete

# Collect all candidate histidine NE2 atoms (common variants)
set his_sel [atomselect $mol "resname HSD HSE HIE HIS HSP and name NE2"]
set nhis [$his_sel num]
if {$nhis == 0} { $his_sel delete ; puts "No HSD/HIS NE2 atoms found."; return }
set his_resid    [$his_sel get resid]
set his_resname  [$his_sel get resname]
set his_segid    [$his_sel get segid]
set his_segname  [$his_sel get segname]
set his_chain    [$his_sel get chain]
set his_xyz      [$his_sel get {x y z}]
$his_sel delete

# Parameters
set cutoff 2.5   ;# Å

# ---- Section 1: list HEME/HEC groups (screen only) ----
puts "Found HEME/HEC groups in the structure:"
puts "======================================"
set g 0
foreach rec $fe_info {
    incr g
    lassign $rec he_rn he_res he_segid he_segname he_chain he_idx
    puts "Group $g: $he_rn:$he_res (segid=$he_segid, chain=$he_chain)"
}
puts "Total HEME/HEC groups found: $g\n"

# ---- Section 2: per-heme nearest histidines ----
puts "HEME/HEC to HSD/HIS Distance Analysis:"
puts "====================================="

set heme_best {}
set used_pairs {}
set excluded_pairs {}

# Track unique patch lines to avoid duplicates
array set seen {}

set hidxs [list]
for {set i 0} {$i < $nhis} {incr i} { lappend hidxs $i }

foreach rec $fe_info {
    lassign $rec he_rn he_res he_segid he_segname he_chain he_idx
    set he_seg [__seg $he_segid $he_segname]
    set fe_xyz [__xyz1 $mol "index $he_idx"]
    if {$fe_xyz eq {}} { continue }

    # compute all distances
    set cands {}
    foreach i $hidxs {
        set d [__dist $fe_xyz [lindex $his_xyz $i]]
        set rresid  [lindex $his_resid   $i]
        set rname   [lindex $his_resname $i]
        set rseg    [__seg [lindex $his_segid $i] [lindex $his_segname $i]]
        set rchain  [lindex $his_chain   $i]
        lappend cands [list $d $rresid $rname $rseg $rchain]
    }
    set cands [lsort -real -increasing -index 0 $cands]

    puts ""
    puts "Analyzing $he_rn:$he_res (segid: $he_segid, chain: $he_chain):"
    set top2 [lrange $cands 0 1]
    set keep_for_patches {}
    set rank 1
    foreach c $top2 {
        lassign $c d rresid rname rseg rchain
        puts [format "Closest %d HSD/HIS NE2: {%s %s} at %.6f Å" $rank $rresid $rname $d]
        if {$d <= $cutoff} {
            lappend keep_for_patches [list $d $rresid $rname $rseg $rchain]
            lappend used_pairs [list $rresid $rseg]
        } else {
            lappend excluded_pairs [list $rresid $rseg]
        }
        incr rank
    }
    dict set heme_best [list $he_seg $he_res] $keep_for_patches
}
puts "\nDistance analysis complete.\n"

# ---- Section 2b: Summary (screen only) ----
puts "Summary of histidine residues for patches (cutoff = ${cutoff} Å):"
puts "----------------------------------------------------------------"
puts "Used for patches:"
if {[llength $used_pairs] == 0} {
    puts "  (none)"
} else {
    foreach pair [lsort -unique $used_pairs] {
        lassign $pair rresid rseg
        puts "  ${rresid}/${rseg}"
    }
}
puts "Excluded (> cutoff):"
if {[llength $excluded_pairs] == 0} {
    puts "  (none)"
} else {
    foreach pair [lsort -unique $excluded_pairs] {
        lassign $pair rresid rseg
        puts "  ${rresid}/${rseg}"
    }
}
puts ""

# ---- Section 3: Emit exact patch commands (screen + file) ----
# Format per line:  patch PHEM <hemeSeg:hemeResid> <hisSeg:hisResid>
set fh [open $patch_file w]
puts $fh "# Auto-generated histidine patches"
puts $fh "# Source molecule: $topname_raw"
puts $fh "# Format: patch PHEM <hemeSeg:hemeResid> <hisSeg:hisResid>"
puts "Step 4: Patches (also written to $patch_file)"
puts [format "%-6s %-12s %-12s" "Patch" "Seg:Resid(HEME)" "Seg:Resid(HIS)"]

# Stable iteration order: by heme, then by distance
foreach rec $fe_info {
    lassign $rec he_rn he_res he_segid he_segname he_chain he_idx
    set he_seg [__seg $he_segid $he_segname]
    set kept {}
    if {[dict exists $heme_best [list $he_seg $he_res]]} {
        set kept [dict get $heme_best [list $he_seg $he_res]]
        set kept [lsort -real -increasing -index 0 $kept]
    }
    foreach k $kept {
        lassign $k d rresid rname rseg rchain
        set lhs [__fmt_sr $he_seg $he_res]
        set rhs [__fmt_sr $rseg   $rresid]
        set line [format "patch PHEM %s %s" $lhs $rhs]
        # De-duplicate identical lines across hemes if any
        if {![info exists seen($line)]} {
            set seen($line) 1
            puts [format "%-6s %-12s %-12s" "PHEM" $lhs $rhs]
            puts $fh $line
        }
    }
}
close $fh

if {[array size seen] == 0} {
    puts ">> No patch lines were generated (no HIS within cutoff)."
    # still leave an empty file with header for transparency
}

puts ">> Patch file written: $patch_file"

