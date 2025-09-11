# === File 3: Rename HIS-family (HIS/HSD/HSE/HSP/HIE) → HSO using patch map ===

proc _warn {m} { puts "!! $m" }
set his_aliases {HIS HSD HSE HSP HIE}  ;# accepted histidine forms

# --- CLI arg fetch ---
proc _arg {flag {def ""}} {
    set i [lsearch $::argv $flag]
    if {$i < 0 || $i == [expr {[llength $::argv]-1}]} { return $def }
    return [lindex $::argv [expr {$i+1}]]
}

# --- Parse map file ---
proc _parse_map {mapfile} {
    set fh [open $mapfile r]
    set out {}
    while {[gets $fh line] >= 0} {
        set line [string trim [string map {"\r" ""} $line]]
        if {$line eq ""} continue
        if {[regexp {^\s*#} $line]} continue
        set hash [string first "#" $line]
        if {$hash >= 0} {
            set line [string trim [string range $line 0 [expr {$hash-1}]]]
            if {$line eq ""} continue
        }
        set parts [split $line]
        if {[llength $parts] < 4} continue
        if {[string tolower [lindex $parts 0]] ne "patch"} continue
        if {[string toupper [lindex $parts 1]] ne "PHEM"} continue
        set rhs [lindex $parts 3]
        if {![regexp {^([A-Za-z0-9_]+):([A-Za-z0-9\-]+)$} $rhs -> seg resid]} continue
        lappend out [list $seg $resid]
    }
    close $fh
    return [lsort -unique $out]
}

# --- Rewrite PSF ---
proc _rewrite_psf_resnames {inpsf outpsf pairs aliases} {
    array set targets {}
    foreach p $pairs { lassign $p s r ; set targets("${s}|${r}") 1 }

    set fin  [open $inpsf r]
    set fout [open $outpsf w]
    set in_atoms 0
    set atoms_left 0
    set saw_natom 0

    while {[gets $fin line] >= 0} {
        if {!$in_atoms && [regexp {^\s*(\d+)\s+!NATOM} $line -> nA]} {
            set atoms_left $nA
            set in_atoms 1
            set saw_natom 1
            puts $fout $line
            continue
        }
        if {$in_atoms && $atoms_left > 0} {
            set toks [regexp -all -inline {\S+} $line]
            if {[llength $toks] >= 7} {
                set segid [lindex $toks 1]
                set resid [lindex $toks 2]
                set rname [lindex $toks 3]
                if {[info exists targets("${segid}|${resid}")] && [lsearch -exact $aliases $rname] >= 0} {
                    set toks [lreplace $toks 3 3 "HSO"]
                    set line [join $toks " "]
                }
            }
            puts $fout $line
            incr atoms_left -1
            if {$atoms_left == 0} { set in_atoms 0 }
            continue
        }
        puts $fout $line
    }
    close $fin
    close $fout
    if {!$saw_natom} { error "PSF parse error: !NATOM block not found in $inpsf" }
}

# ---------- MAIN ----------
set mol [molinfo top]
if {$mol < 0} { error "No molecule loaded" }

set topname_raw [molinfo $mol get name]
set topname     [file rootname $topname_raw]
if {$topname eq ""} { set topname "top" }

set map_file "${topname}_hispatch.txt"
if {![file exists $map_file]} { error "Map file not found: $map_file" }

# Input PSF/PDB (args or guesses)
set in_psf [_arg -psf ""]
set in_pdb [_arg -pdb ""]
if {$in_psf eq ""} {
    foreach g [list "${topname}.psf" "${topname}_build.psf"] {
        if {[file exists $g]} { set in_psf $g ; break }
    }
}
if {$in_pdb eq ""} {
    foreach g [list "${topname}.pdb" "${topname}_build.pdb"] {
        if {[file exists $g]} { set in_pdb $g ; break }
    }
}
if {$in_psf eq "" || $in_pdb eq ""} {
    error "Could not infer PSF/PDB. Provide with -psf and -pdb."
}
puts ">> Using input:\n   PSF: $in_psf\n   PDB: $in_pdb"

# Parse HIS targets
set his_targets [_parse_map $map_file]
if {[llength $his_targets] == 0} { error "No HIS targets parsed" }
puts ">> Parsed [llength $his_targets] HIS targets"

# Apply renames in-memory
set psf_pairs {}
set n_changed 0
set n_notfound 0
set n_notAlias 0

foreach pair $his_targets {
    lassign $pair seg resid
    set sel [atomselect $mol "segname $seg and resid $resid"]
    if {[$sel num] == 0} {
        $sel delete
        set sel [atomselect $mol "segid $seg and resid $resid"]
    }
    if {[$sel num] == 0} {
        _warn "Not found: $seg:$resid"
        $sel delete
        incr n_notfound
        continue
    }
    set curr_names [lsort -unique [$sel get resname]]
    if {[llength $curr_names] == 1 && [lsearch -exact $his_aliases [lindex $curr_names 0]] >= 0} {
        $sel set resname HSO
        set segids [lsort -unique [$sel get segid]]
        lappend psf_pairs [list [lindex $segids 0] $resid]
        incr n_changed
    } else {
        _warn "$seg:$resid not HIS alias (found $curr_names)"
        incr n_notAlias
    }
    $sel delete
}

set psf_pairs [lsort -unique $psf_pairs]

# Write outputs
set out_psf "${topname}_HSO.psf"
set out_pdb "${topname}_HSO.pdb"

set allsel [atomselect $mol all]
$allsel writepdb $out_pdb
$allsel delete
puts ">> Wrote $out_pdb"

_rewrite_psf_resnames $in_psf $out_psf $psf_pairs $his_aliases
puts ">> Wrote $out_psf"

# Load new molecule (keep original visible)
set newmol [mol new $out_psf type psf waitfor all]
mol addfile $out_pdb type pdb -molid $newmol waitfor all
mol rename $newmol "${topname}_HSO"

puts "== DONE =="
puts "Renamed: $n_changed, Skipped not found: $n_notfound, Skipped not HIS-alias: $n_notAlias"
puts "Loaded new molecule: [molinfo $newmol get name]"
