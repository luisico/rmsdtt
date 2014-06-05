# Test measure_rmsd patch

set mol top
set sel_text "residue 1 to 5 and name C N CA"

set sel1 [atomselect $mol $sel_text frame 1]
set sel2 [atomselect $mol $sel_text frame 2]
lassign [measure rmsd $sel1 $sel2 byatom] r_hack_global r_hack_byatom
lassign [measure rmsd $sel1 $sel2 byres] r_hack_global2 r_hack_byres

set r_stock_global [measure rmsd $sel1 $sel2]

set r_stock_byres {}
set residues [lsort -unique -integer [$sel1 get residue]]
set nres [llength $residues]
for {set i 0} {$i < $nres} {incr i} {
  set sel1_res [atomselect $mol "$sel_text and residue [lindex $residues $i]" frame 1]
  set sel2_res [atomselect $mol "$sel_text and residue [lindex $residues $i]" frame 2]
  lappend r_stock_byres [measure rmsd $sel1_res $sel2_res]
}

set natoms [$sel1 num]
for {set i 0} {$i < $natoms} {incr i} {
  set sel1_at [atomselect $mol "index [lindex [$sel1 get index] $i]" frame 1]
  set sel2_at [atomselect $mol "index [lindex [$sel2 get index] $i]" frame 2]
  lappend r_stock_byatom [measure rmsd $sel1_at $sel2_at]
}

puts ""
for {set i 0} {$i < $nres} {incr i} {
  if {[::tcl::mathop::!= [lindex $r_hack_byres $i] [lindex $r_stock_byres $i]]} {
    puts "Error with res $i -> [lindex $r_hack_byres $i] --- [lindex $r_stock_byres $i]"
  }
}

puts ""
for {set i 0} {$i < $natoms} {incr i} {
  if {[::tcl::mathop::!= [lindex $r_hack_byatom $i] [lindex $r_stock_byatom $i]]} {
    puts "atom $i -> [lindex $r_hack_byatom $i] --- [lindex $r_stock_byatom $i]"
  }
}
