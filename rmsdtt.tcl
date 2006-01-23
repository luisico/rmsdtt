##
## RMSD Trajectory Tool
##
## A GUI interface for RMSD alignment
##
## Authors: Alexander Balaeff, Pascal Mercier, John Stone
##
## Modifed by Luis Gracia to include trajectories
##
##

##
## Example code to add this plugin to the VMD extensions menu:
##
#  if { [catch {package require rmsdtt} msg] } {
#    puts "VMD RMSDTT package could not be loaded:\n$msg"
#  } elseif { [catch {menu tk register "rmsdtt" rmsdtt} msg] } {
#    puts "VMD RMSDTT could not be started:\n$msg"
#  }

# Plugin-ized version of a modified version of the VMD RMSD script

# Authors:
#   Alexander Balaeff, PhD, TB Group, UIUC
#   Pascal Mercier, PhD, Biochemistry Dept, Univ. of Alberta, Edmonton, Canada
#   John Stone, TB Group, UIUC
#   Luis Gracia, PhD, Weill Medical College, Cornel University, NY

# Tell Tcl that we're a package and any dependencies we may have
package provide rmsdtt 1.0

namespace eval ::RMSDTT:: {
  namespace export rmsdtt
  variable w                         ;# handle to main window
  variable rms_ave                   ;# array of rms ave
  variable rms_list                  ;# list of rms values
  variable pdb_list                  ;# list of pdb files
  variable bb_only                   ;# backbon-only settings (C CA N)
  variable trace_only                ;# Trace settings
  variable noh                       ;# No hydrogens
  variable skip_sw
  variable skip_ini
  variable skip_steps
  variable rms_sel                   ;# rms selection text
  variable rmsd_base                 ;# which molecule is the reference for rmsd calc.
  variable tot_rms                   ;# total rms value 
  variable RMSDhistory		     ;# history of selections
  variable rms_values
  
  option add *Font {Helvetica -12}
  # Original color scheme by Alexander/Pascal
  #  variable ftr_bgcol   \#2a4
  #  variable ftr_fgcol   \#ff0
  #  variable calc_bgcol  \#24a
  #  variable calc_fgcol  \#ff0
  #  variable sel_bgcol   \#333
  #  variable sel_fgcol   \#ccc
  #  variable act_bgcol   \#0b4
  #  variable act_fgcol   \#f00
  #  variable entry_bgcol \#a64da6
  #  variable but_abgcol  \#26c
  #  variable but_bgcol   \#338
  #  variable scr_trough  \#126

  # Calmer color scheme that looks like other plugins
  variable ftr_bgcol   \#d9d9d9
  variable ftr_fgcol   \#000
  variable calc_bgcol  \#d9d9d9
  variable calc_fgcol  \#000
  variable sel_bgcol   \#ff0
  variable sel_fgcol   \#000
  variable act_bgcol   \#d9d9d9
  variable act_fgcol   \#000
  variable entry_bgcol \#fff
  variable but_abgcol  \#d9d9d9
  variable but_bgcol   \#d9d9d9
  variable scr_trough  \#c3c3c3
}

proc ::RMSDTT::get_molid_from_pdb_list_to_operate_on {} { 
  variable pdb_list
  set all_mols {}
  set listedpdbs [$pdb_list get 0 end]

  foreach pdb $listedpdbs { lappend all_mols [lindex $pdb 0] }
    return $all_mols
}


proc ::RMSDTT::molactive {} {
  set mol_active_list {}
  set all_mols [molinfo list]
    
  foreach i $all_mols  {
    if { [molinfo $i get active] } {
      set mol_active_list [concat $mol_active_list $i]
    }
  }

  return $mol_active_list
}

proc ::RMSDTT::align { sel1 thisframe sel2} {
  # tries to align sel1 with sel2, using the atoms of each selection.
  # If sel1 and sel2 contain a different number of atoms, it fails.

  $sel1 frame $thisframe
  set tmatrix [measure fit $sel1 $sel2]
  set molid [$sel1 molid]
  set move_sel [atomselect $molid "all" frame $thisframe]
  $move_sel move $tmatrix
  return $tmatrix
}


proc ::RMSDTT::align_all {sel_text rmsd_base frames_ref} {
  variable pdb_list
  variable frames_sw
  set all_mols [get_molid_from_pdb_list_to_operate_on]
  
  
  foreach i $all_mols  {set sel($i) [atomselect $i $sel_text]}
  # do the top molecule as well, just in case it's not in the list
  set mol_on_top [molinfo top] 
  set sel($mol_on_top) [atomselect $mol_on_top $sel_text]

  switch $rmsd_base {
    top {set reference_for_alignment $mol_on_top}
    selected {set reference_for_alignment [lindex [$pdb_list get active] 0]}
  }

  set sel_ref [atomselect $reference_for_alignment $sel_text frame $frames_ref]

  foreach i $all_mols  {
    if {$frames_sw == 0} {
      if { $i != $reference_for_alignment } {
	set j [molinfo $i get frame]
	align $sel($i) $j $sel($reference_for_alignment)
      }
    } else {
      set jmax [molinfo $i get numframes]
      for {set j 0} {$j < $jmax} {incr j} {
	if { $i == $reference_for_alignment && $j == $frames_ref } {
	} else {
  	  align $sel($i) $j $sel_ref
	}
      }
    }
  }
}


proc ::RMSDTT::ini_zero {length} {
  if {$length<=0} return {}

  set half_l [expr $length>>1]
  set s {{0 0 0}}
  for {set i 1} {$i<$half_l} {set i [expr $i<<1]} {
    set s [concat $s $s]
  }

  set more [expr $length-$i-1]
  set s [concat $s [lrange $s 0 $more]]

  return $s
}


proc ::RMSDTT::ave_struc {sel_text} {
  set all_mols [get_molid_from_pdb_list_to_operate_on]
  set n_mols [llength $all_mols]
  set factor [expr 1./$n_mols]
  set mol_on_top [molinfo top]

  foreach i $all_mols {
    set all_coor($i) [[atomselect $i $sel_text] get {x y z}]
  }
  
  set ave_coor {}
  set m0 [lindex $all_mols 0]
  set len [llength $all_coor($m0)]

  for {set j 0} {$j<$len} {incr j} {
    set b [veczero]
    foreach i $all_mols {
      set b [vecadd $b [lindex $all_coor($i) $j]]
    }
    lappend ave_coor [vecscale $b $factor]
  }

  return $ave_coor
}


proc ::RMSDTT::compute_rms {sel_text {frames_ref 0}} {
  variable rms_ave 
  variable rms_disp 
  variable pdb_list
  variable rms_values
  variable rmsd_base
  
  variable skip_sw
  variable skip_ini
  variable skip_steps
  variable frames_sw
  variable frames_all
  variable rms_aveframe
  # delete this variable mol_ref?
  variable mol_ref
  
  array unset rms_values

  set all_mols [get_molid_from_pdb_list_to_operate_on]
  set mol_on_top [molinfo top]
  set n_mols [llength $all_mols]
  # 1 will be subtracted from n_mols if the RMSD is calculated against the top molecule 
  # and if the top molecule is part of the listed pdbs
  # or the selected item in the pdb list, see below
    
  if {$sel_text == ""} {
    showMessage "Selection is empty selection!"
    return -code return
  }

  # Check if the top molecule is listed in the pdb list.
  set topisthere [lsearch $all_mols $mol_on_top]
  foreach i $all_mols {
    set jmax [molinfo $i get numframes]
    set natoms($i) [[atomselect $i $sel_text frame 0] num]
  }
  # Check same number of atoms
  foreach i $all_mols {
    foreach j $all_mols {
      if {$i < $j} {
	if {$natoms($i) != $natoms($j)} {
	  showMessage "Selections differ for molecules $i ($natoms($i)) and $j ($natoms($j))"
	  return -code return
	}
      }
    }
  }
  if {$topisthere < 0} {
    set jmax [molinfo $mol_on_top get numframes]
    set natomstop [[atomselect $mol_on_top $sel_text frame 0] num]
    foreach i $all_mols {
      if {$natoms($i) != $natomstop } {
	showMessage "Selections differ for molecules $i ($natoms($i)) and top ($natomstop)"
	return -code return
      }
    }
  }
    
  
  switch $rmsd_base {
    top {
      set mol_ref $mol_on_top
      if {$frames_sw == 0} {
	set frames_ref [molinfo $mol_ref get frame]
      } else {
	set j $frames_ref
	set n0 [molinfo $mol_ref get numframes]
	if {$frames_ref >= $n0} {
	  showMessage "Frame ref out of range (max is [expr $n0-1])"
	  return -code return
	}
      }
      if {$topisthere >= 0 && $frames_sw == 0} {incr n_mols -1}
    }
  
    ave {
      set mol_ref "ave"
      set ref_coor {}

      if {$frames_sw == 0} {
	set m0 [lindex $all_mols 0]
	set len [llength [[atomselect $m0 $sel_text frame 0] get {x y z}]]
	set factor [expr 1./$n_mols]
	for {set j 0} {$j<$len} {incr j} {
	  set b [veczero]
	  foreach i $all_mols {
	    set k [molinfo $i get frame]
	    set b [vecadd $b [lindex [[atomselect $i $sel_text frame $k] get {x y z}] $j]]
	  }
	  lappend ref_coor [vecscale $b $factor]
	}
      } else {
	set len [llength [[atomselect $mol_on_top $sel_text frame 0] get {x y z}]]
	set nframes [molinfo $mol_on_top get numframes]
	set factor [expr 1./$nframes]
	for {set j 0} {$j<$len} {incr j} {
	  set b [veczero]
	  for {set i 0} {$i < $nframes} {incr i} {
	    set b [vecadd $b [lindex [[atomselect $mol_on_top $sel_text frame $i] get {x y z}] $j]]
	  }
	  lappend ref_coor [vecscale $b $factor]
	}

      }

      puts [lindex $ref_coor 0]
    }
  
    selected {
      set index [lindex [$pdb_list get active] 0]
      set mol_ref $index
      if {$frames_sw == 0} {
	incr n_mols -1
	set j [molinfo $index get frame]
      } else {
	set j $frames_ref
	set n0 [molinfo $mol_ref get numframes]
	if {$frames_ref >= $n0} {
	  showMessage "Frame ref out of range (max is [expr $n0-1])"
	  return -code return
	}
      }
    }
  }

  set tot_rms 0
  
  # for each molecule
  foreach i $all_mols {
    set jmax [molinfo $i get numframes]
    if {$frames_sw == 0} { set jmax 1 }

    set rms_ave($i) 0

    # for each frame
    if {$frames_sw == 1 && $skip_sw == 1} {
      set skipped_tot $skip_ini
      set skipped $skip_steps
    }
    for {set j 0} {$j < $jmax} {incr j} {
      if {$frames_sw == 1 && $skip_sw == 1} {
	if {$j < $skip_ini} {continue}
	if {$skipped < $skip_steps} {
	  incr skipped
	  incr skipped_tot
	  continue
	} else {
	  set skipped 0
	}
      }
      if {$frames_sw == 0} {set j [molinfo $i get frame]}
      if {$rmsd_base == "ave"} {
	set rmsd [get_rmsd_ave $ref_coor $i $j $sel_text]
      } else {
	set rmsd [get_rmsd $mol_ref $frames_ref $i $j $sel_text]
      }
      set rms_values($i,$j) $rmsd
      set rms_ave($i) [expr $rms_ave($i) + $rmsd]
    }
    
    if {$frames_sw == 1 && $skip_sw == 1} {set jmax [expr $jmax - $skipped_tot]}
    if {$frames_sw == 1 && $mol_ref == $i && $jmax > 1} {incr jmax -1}
    set rms_ave($i) [expr $rms_ave($i)/$jmax]
    set tot_rms [expr $tot_rms + $rms_ave($i)]
  }

  set tot_rms [expr $tot_rms/$n_mols]

  return $tot_rms
}

proc ::RMSDTT::get_rmsd { mol1 frame1 mol2 frame2 sel_text } {
  set sel1 [atomselect $mol1 $sel_text frame $frame1] 
  set sel2 [atomselect $mol2 $sel_text frame $frame2]
  set rmsd [measure rmsd $sel1 $sel2]
  return $rmsd
}

proc ::RMSDTT::get_rmsd_ave { ref_coor mol2 frame2 sel_text } {
  set sel2 [atomselect $mol2 $sel_text frame $frame2]
  set coor2 [$sel2 get {x y z}]
  set numatoms [llength $ref_coor]

  set rmsd 0
  # for each atom
  for {set k 0} {$k < $numatoms} {incr k} {
    set v1 [lindex $ref_coor $k]
    set v2 [lindex $coor2 $k]
    set rmsd [expr $rmsd + [veclength2 [vecsub $v1 $v2]]]
  }
  set rmsd [expr $rmsd/$numatoms]
  set rmsd [expr sqrt($rmsd)]
  
  return $rmsd
}

proc ::RMSDTT::reveal_rms {} {
  variable rms_ave 
  variable rms_list 
 
  set all_mols [get_molid_from_pdb_list_to_operate_on]
  $rms_list delete 0 end
  foreach i $all_mols {
    $rms_list insert end $rms_ave($i)
  }
}



proc ::RMSDTT::two_scroll args {
  variable pdb_list 
  variable rms_list
  eval "$pdb_list yview $args"
  eval "$rms_list yview $args"
}


proc ::RMSDTT::onlyactive args {
  variable pdb_list 
  variable rms_list
  $pdb_list delete 0 end
  $rms_list delete 0 end
  for {set i 0} {$i < [molinfo num]} {incr i} {
    set molid [molinfo index $i]
    if {[molinfo $molid get active]} {$pdb_list insert end [format "%-4s%-10s" $molid [molinfo $molid get name]]}      
  }
}


# This gets called by VMD the first time the menu is opened.
proc rmsdtt_tk_cb {} {
  variable foobar
  # Don't destroy the main window, because we want to register the window
  # with VMD and keep reusing it.  The window gets iconified instead of
  # destroyed when closed for any reason.
  #set foobar [catch {destroy $::RMSDTT::w  }]  ;# destroy any old windows

  ::RMSDTT::rmsdtt   ;# start the RMSD Tool
  return $RMSDTT::w
}

proc RMSDTT::showMessage {mess} {
  bell
  toplevel .messpop 
  grab .messpop
  option add *Font {Helvetica -12 bold}
  wm title .messpop "Warning"
    message .messpop.msg -relief groove -bd 2 -text $mess -aspect 400 -justify center -padx 20 -pady 20
  
  button .messpop.okb -text OK -command {destroy .messpop ; return 0}
  pack .messpop.msg .messpop.okb -side top 
}


proc RMSDTT::set_sel {} {
  variable w
  variable bb_only
  variable trace_only
  variable noh

#  set a [$w.top.left.inner.selfr.sel get 1.0 end]
#  puts "a <$a>"
  regsub -all "\#.*?\n" [$w.top.left.inner.selfr.sel get 1.0 end] "" temp1
  regsub -all "\n" $temp1 " " temp2
  regsub -all " $" $temp2 "" temp3
#  puts "c <$temp3>"
  if { $trace_only } {
    append rms_sel "($temp3) and name CA"
  } elseif { $bb_only } {
    append rms_sel "($temp3) and name C CA N"
  } elseif { $noh } {
    append rms_sel "($temp3) and noh"
  } else {
    append rms_sel $temp3
  }
  return $rms_sel
}

proc RMSDTT::ListHisotryPullDownMenu {} {
  variable RMSDhistory
  set rc_win .rmsdtt.top.left.inner
  $rc_win.selectionhistory.m delete 0 end
  foreach sel $RMSDhistory {
  	$rc_win.selectionhistory.m add command -label $sel \
  	 -command [list RMSDTT::chooseHistoryItem $sel]
  }
}


proc ::RMSDTT::chooseHistoryItem {sel} {
  variable rms_sel
  set rms_sel $sel
}


proc ::RMSDTT::saveDataAll {frame_ref file_out_id time_ref} {
  variable rms_values
  variable time_sw
  variable time_ini
  variable time_step

  set all_mols [get_molid_from_pdb_list_to_operate_on]
  set n_mols [llength $all_mols]
  
  if {$n_mols > 1} {
    foreach i $all_mols {
      set nframes [molinfo $i get numframes]
      for {set j 0} {$j < $nframes} {incr j} {
	if {$time_sw} {
	  set time [expr $time_ini + $time_step * $j]
	  puts $file_out_id [format "%8.2f\t%3d\t%8.2f\t%10.7f" $frame_ref $i $time $rms_values($i,$j)]
	} else {
	  puts $file_out_id [format "%5d\t%3d\t%5d\t%10.7f" $frame_ref $i $j $rms_values($i,$j)]
	}
      }
    }
  } else {
    set nframes [molinfo $all_mols get numframes]
    for {set j 0} {$j < $nframes} {incr j} {
      puts -nonewline $file_out_id [format "%6.2f " $rms_values($all_mols,$j)]
    }
    puts $file_out_id ""
  }
  
}


proc ::RMSDTT::saveData {} {
  variable rms_values
  variable file_out
  variable time_sw
  variable time_ini
  variable time_step
  variable skip_sw
  variable skip_ini
  variable skip_steps
  
  set file_out_id [open $file_out w]
  fconfigure $file_out_id -buffering line

  set all_mols [get_molid_from_pdb_list_to_operate_on]
  set n_mols [llength $all_mols]

  if {$time_sw} {
    puts -nonewline $file_out_id "time"
  } else {
    puts -nonewline $file_out_id "frame"
  }
  set jmaxmax 1
  foreach i $all_mols {
    set jmax($i) [molinfo $i get numframes]
    if {$jmax($i) > $jmaxmax} {set jmaxmax $jmax($i)}
    puts -nonewline $file_out_id [format " %10s" "mol$i"]
  }
  puts $file_out_id ""
  
  for {set j 0} {$j < $jmaxmax} {incr j} {
    if {![info exists rms_values([lindex $all_mols 0],$j)]} {continue}
    if {$time_sw} {
      set time [expr $time_ini + $time_step * $j]
      puts -nonewline $file_out_id [format "%8.2f" $time]
    } else {
      puts -nonewline $file_out_id [format "%5d" $j]
    }
    foreach i $all_mols {
      set jmax($i) [molinfo $i get numframes]
#      puts -nonewline "$j $i $jmax($i) $jmaxmax"
      if {$j < $jmax($i)} {
#	puts [format " %10.7f" $rms_values($i,$j)]
	puts -nonewline $file_out_id [format " %10.7f" $rms_values($i,$j)]
      } else {
#	puts [format " %10s" {NA}]
	puts -nonewline $file_out_id [format " %10s" {NA}]
      }
     }
     puts $file_out_id ""
  }
  puts $file_out_id ""
    
  close $file_out_id
  
}

proc ::RMSDTT::tempfile {prefix suffix} {
  # From wiki.tcl.tk/772
  set chars "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
  set nrand_chars 10
  set maxtries 10
  set access [list RDWR CREAT EXCL TRUNC]
  set permission 0600
  set channel ""
  set checked_dir_writable 0
  set mypid [pid]
  for {set i 0} {$i < $maxtries} {incr i} {
    set newname $prefix
    for {set j 0} {$j < $nrand_chars} {incr j} {
      append newname [string index $chars \
			[expr ([clock clicks] ^ $mypid) % 62]]
    }
    append newname $suffix
    if {[file exists $newname]} {
      after 1
    } else {
      if {[catch {open $newname $access $permission} channel]} {
	if {!$checked_dir_writable} {
	  set dirname [file dirname $newname]
	  if {![file writable $dirname]} {
	    error "Directory $dirname is not writable"
	  }
	  set checked_dir_writable 1
	}
      } else {
	# Success
	return [list $newname $channel]
      }
    }
  }
  if {[string compare $channel ""]} {
    error "Failed to open a temporary file: $chanel"
  } else {
    error "Failed to find an unused temporary file name"
  }
}

proc ::RMSDTT::doRmsd {} {
  variable w
  variable rmsd_base
  variable rms_sel
  variable frames_sw
  variable frames_all
  variable frames_ref
  variable file_out_sw
  variable file_out
  variable plot_sw
  variable RMSDhistory
  variable tot_rms
  variable rms_values
  variable time_sw
  variable time_ini
  variable time_step

  if {$frames_sw && $file_out_sw && $file_out == ""} {
    showMessage "Filename is missing!"
    return -code return
  }

  set rms_sel [set_sel]

  if {$frames_all} {
    switch $rmsd_base {
      top {
	set mol_ref [molinfo top]
	set nframes [molinfo $mol_ref get numframes]
      }
      ave {
	set mol_ref "ave"
	set nframes [molinfo [molinfo top] get numframes]
      }
      selected {
	set index [lindex [$pdb_list get active] 0]
	set mol_ref $index
	set nframes [molinfo $mol_ref get numframes]
      }
    }
    set all_mols [get_molid_from_pdb_list_to_operate_on]
    set n_mols [llength $all_mols]
    if {$frames_sw && $file_out_sw} {
      set file_out_id [open $file_out w]
      fconfigure $file_out_id -buffering line
      if {$n_mols > 1} {
	if {$time_sw} {
	  puts $file_out_id "time_ref mol\ttime\trmsd"
	} else {
	  puts $file_out_id "frame_ref mol\tframe\trmsd"
	}
      }
    }
    for {set k 0} {$k < $nframes} {incr k} {
      set tot_rms [compute_rms $rms_sel $k]
      if {$frames_sw} {
	if {$file_out_sw} {
	  set time_ref [expr $time_ini + $time_step * $k]
	  saveDataAll $k $file_out_id $time_ref
	}
	#if {$n_mols == 1} {
	#  if {$plot_sw} {doPlotAll}
	#}
      }
      
    }
    if {$frames_sw && $file_out_sw} {
      close $file_out_id
    }
    
  } else {
    set tot_rms [compute_rms $rms_sel $frames_ref]
    reveal_rms
    if {$frames_sw} {
      if {$file_out_sw} {saveData}
      if {$plot_sw} {doPlot}
    }
  }
  lappend RMSDhistory $rms_sel
  ListHisotryPullDownMenu
}

proc ::RMSDTT::doAlign {} {
  variable w
  variable rmsd_base
  variable rms_sel
  variable frames_ref
  
  if {$rmsd_base=="ave"} {
    showMessage "Average option not available for Alignment in this version"
    return -code return
  }
  set rms_sel [set_sel]
  align_all $rms_sel $rmsd_base $frames_ref
}

proc ::RMSDTT::doPlot {} {
  variable rms_values
  variable rms_sel
  variable time_sw
  variable time_ini
  variable time_step
  global tcl_platform
  
  set all_mols [get_molid_from_pdb_list_to_operate_on]
  set n_mols [llength $all_mols]

  switch $tcl_platform(platform) {
    unix {
      # parray rms_values
      # set pipe_id [open "| xmgrace -pipe &" w]
      
      set f [tempfile rmsdtt .tmp]
      set filename [lindex $f 0]
      set pipe_id [lindex $f 1]
      fconfigure $pipe_id -buffering line
      
      puts $pipe_id "@ page size 576, 432"
      puts $pipe_id "@ g0 on"
      puts $pipe_id "@ with g0"
      puts $pipe_id "@ subtitle \"$rms_sel\""
      if {$time_sw} {
	puts $pipe_id "@ title \"Rmsd vs Time\""
	puts $pipe_id "@ xaxis  label \"Time (ps)\""
      } else {
	puts $pipe_id "@ title \"Rmsd vs Frame\""
	puts $pipe_id "@ xaxis  label \"Frame\""
      }
      puts $pipe_id "@ yaxis  label \"Rmsd (A)\""
      puts $pipe_id "@ TYPE xy"
      puts $pipe_id "@ view 0.15, 0.15, 0.75, 0.85"
      puts $pipe_id "@ legend on"
      puts $pipe_id "@ legend box on"
      set k 0
      foreach i $all_mols {
	set iname [molinfo $i get name]
	puts $pipe_id "@ s$k legend \"$iname ($i)\""
	set jmax [molinfo $i get numframes]
	if {$jmax == 1} {
	  puts $pipe_id "@ s$k symbol 1"
	}
	for {set j 0} {$j < $jmax} {incr j} {
	  if {![info exists rms_values($i,$j)]} {continue}
	  if {$time_sw} {
	    set time [expr $time_ini + $time_step * $j]
	    puts $pipe_id "$time $rms_values($i,$j)"
	  } else {
	    puts $pipe_id "$j $rms_values($i,$j)"
	  }
	}
	puts $pipe_id ""
	set k [expr $k+1]
	
      }
      #  foreach i [lsort -dictionary [array names rms_values]] {
      #    #puts $i:$rms_values($i)
      #  }
      
      close $pipe_id
      set status [catch {exec xmgrace $filename &} msg]
      if { $status } {
	showMessage "Could not open xmgrace. Error returned:\n $msg"
	file delete -force $filename
	return -code return
      } 
    }
    windows {
      package require tcom
      
      set excel [::tcom::ref createobject "Excel.Application"]
      #set excel [::tcom::ref getactiveobject "Excel.Application"]
      $excel Visible 1
      
      set workbooks [$excel Workbooks]
      set workbook [$workbooks Add]
      set worksheets [$workbook Worksheets]
      [$worksheets Item [expr 2]] Delete
      [$worksheets Item [expr 2]] Delete
      set worksheet [$worksheets Item [expr 1]]
      $worksheet Name "RMSDTT data"
      
      set cells [$worksheet Cells]
      set jmaxmax 1
      foreach i $all_mols {
	set jmax($i) [molinfo $i get numframes]
	if {$jmax($i) > $jmaxmax} {set jmaxmax $jmax($i)}
      }
      if {$time_sw} {
	$cells Item 1 1 "Time (ps)"
	set time $time_ini
	for {set j 0} {$j < $jmaxmax} {incr j} {
	  $cells Item [expr $j+2] 1 $time
	  set time [expr $time + $time_step]
	}
	set exceltitle "Rmsd vs Time"
      } else {
	$cells Item 1 1 "Frame"
	for {set j 0} {$j < $jmaxmax} {incr j} {
	  $cells Item [expr $j+2] 1 $j
	}
	set exceltitle "Rmsd vs Frame"
      }
      set k 1
      foreach i $all_mols {
	incr k
	set iname [molinfo $i get name]
	$cells Item 1 $k "$iname ($i)"

	for {set j 0} {$j < $jmax($i)} {incr j} {
	  if {![info exists rms_values($i,$j)]} {continue}
	  $cells Item [expr $j+2] $k $rms_values($i,$j)
	}
      }

      set charts [$workbook Charts]
      set chart [$charts Add]
      $chart Name "RMSDTT graph"
      set endrange [int2word $k]
      append endrange [expr $jmaxmax+2]
      $chart ChartWizard [$worksheet Range "A1" $endrange] -4169 [::tcom::na] 2 1 [expr $k-1] 1 "$exceltitle\n($rms_sel)"
      $chart ChartType 75
      [[$chart PlotArea] Interior] ColorIndex 0
      
      set axes [$chart Axes]
      set xaxis [$axes Item 1]
      $xaxis HasMajorGridlines 0
      $xaxis HasTitle 1
      if {$time_sw} {
	[$xaxis AxisTitle] Text "Time (ps)"
      } else {
	[$xaxis AxisTitle] Text "Frame"
      }
      set yaxis [$axes Item 2]
      $yaxis HasMajorGridlines 0
      $yaxis HasTitle 1
      [$yaxis AxisTitle] Text "Rmsd (A)"

#      set interface [::tcom::info interface $xaxis]
#      set properties [$interface properties]
#      foreach property $properties {
#	puts "property $property"
#      }
#      set methods [$interface methods]
#      foreach method $methods {
#			       puts "method [lrange $method 0 2] \{"
#			       set parameters [lindex $method 3]
#			       foreach parameter $parameters {
#				 puts "    \{$parameter\}"
#			       }
#			       puts "\}"
#			     }
      #$excel Quit
    }
  }


}

proc RMSDTT::int2word {int} {
  # http://wiki.tcl.tk/10915
  set alphabet [a-z]
  set word ""
  set la [llength $alphabet]
  while {$int > 0} {
    incr int -1
    set word  [lindex $alphabet [expr {$int % $la}]]$word
    set int   [expr {$int/$la}]
  }
  set word
}
proc a-z {} {list a b c d e f g h i j k l m n o p q r s t u v w x y z}

proc ::RMSDTT::ctrlgui {} {
  variable w
  variable frames_sw
  variable frames_all
  variable file_out_sw
  variable plot_sw
  variable rmsd_base
  variable time_sw
  variable skip_sw
  variable skip_ini
  variable skip_steps

  if {$frames_sw} {
    if {$frames_all} {
      $w.top.right.traj.file.plot config -state disable
    } else {
      $w.top.right.traj.file.plot config -state normal
    }
    $w.top.right.traj.file.0 config -state normal
    if {$file_out_sw} {
      $w.top.right.traj.file.name config -state normal
    } else {
      $w.top.right.traj.file.name config -state disable
    }
    if {$rmsd_base == "ave"} {
      $w.top.right.traj.frames.reflabel config -state disable
      $w.top.right.traj.frames.all config -state disable
      $w.top.right.traj.frames.ref config -state disable
    } else {
      $w.top.right.traj.frames.reflabel config -state normal
      $w.top.right.traj.frames.all config -state normal
      if {$frames_all} {
	$w.top.right.traj.frames.ref config -state disable
      } else {
	$w.top.right.traj.frames.ref config -state normal
      }
    }
    $w.top.right.traj.time.0 config -state normal
    $w.top.right.traj.skip.0 config -state normal
  } else {
    $w.top.right.traj.file.plot config -state disable
    $w.top.right.traj.file.0 config -state disable
    $w.top.right.traj.file.name config -state disable
    $w.top.right.traj.frames.reflabel config -state disable
    $w.top.right.traj.frames.all config -state disable
    $w.top.right.traj.frames.ref config -state disable
    $w.top.right.traj.time.0 config -state disable
    $w.top.right.traj.skip.0 config -state disable
  }

  if {$time_sw && $frames_sw} {
    $w.top.right.traj.time.inilabel config -state normal
    $w.top.right.traj.time.inival config -state normal
    $w.top.right.traj.time.steplabel config -state normal
    $w.top.right.traj.time.stepval config -state normal
  } else {
    $w.top.right.traj.time.inilabel config -state disable
    $w.top.right.traj.time.inival config -state disable
    $w.top.right.traj.time.steplabel config -state disable
    $w.top.right.traj.time.stepval config -state disable
  }

  if {$skip_sw && $frames_sw} {
    $w.top.right.traj.skip.inilabel config -state normal
    $w.top.right.traj.skip.ini config -state normal
    $w.top.right.traj.skip.stepslabel config -state normal
    $w.top.right.traj.skip.steps config -state normal
  } else {
    $w.top.right.traj.skip.inilabel config -state disable
    $w.top.right.traj.skip.ini config -state disable
    $w.top.right.traj.skip.stepslabel config -state disable
    $w.top.right.traj.skip.steps config -state disable
  }

}


proc ::RMSDTT::ctrlbb { obj } {
  variable w
  variable bb_only
  variable trace_only
  variable noh

  if {$obj == "bb"} {
    set trace_only 0
    set noh 0
  } elseif {$obj == "trace"} {
    set bb_only 0
    set noh 0
  } elseif {$obj == "noh"} {
    set trace_only 0
    set bb_only 0
  }
}


proc ::RMSDTT::rmsdtt {} {
  variable w ;# Tk window
  variable bb_only 
  variable trace_only
  variable noh
  variable skip_sw
  variable skip_ini
  variable skip_steps
  variable pdb_list
  variable rms_ave 
  variable rms_list 
  variable rms_sel 
  variable rmsd_base
  variable tot_rms 
  variable rms_disp
  variable RMSDhistory

  variable frames_sw
  variable frames_ref
  variable frames_all
  variable time_sw
  variable time_ini
  variable time_step
  variable file_out_sw
  variable file_out
  variable file_out_id
  variable plot_sw

  variable ftr_bgcol 
  variable ftr_fgcol 
  variable calc_bgcol
  variable calc_fgcol
  variable sel_bgcol
  variable sel_fgcol
  variable act_bgcol   
  variable act_fgcol   
  variable entry_bgcol 
  variable but_abgcol  
  variable but_bgcol   
  variable scr_trough 

  #  Set several parameters.
  set bb_only 0
  set trace_only 0
  set noh 1
  set RMSDhistory ""
  set rmsd_base {top}
  set tot_rms {}

  set skip_sw 0
  set skip_ini 0
  set skip_steps 1

  set frames_sw 1
  set frames_ref 0
  set frames_all 0
  set time_sw 0
  set time_ini 0.0
  set time_step 1.0
  set file_out_sw 0
  set file_out "trajrmsd.dat"
  set plot_sw 1

  # If already initialized, just turn on
  if { [winfo exists .rmsdtt] } {
    wm deiconify $w
    return
  }

  # Create the window for the RMSD calculations.
  set w [toplevel ".rmsdtt" -bg $calc_bgcol]
  wm title $w "RMSD Trajectory Tool"
  wm resizable $w 0 0

  # Create the top part with selection windows and 
  # the "Calculate" button

  # Selection part.
  set calc_top [frame $w.top -bg $calc_bgcol]

  frame $calc_top.left -relief ridge -bd 4 -bg $calc_bgcol
  set rc_win [frame $calc_top.left.inner -bg $calc_bgcol -bd 0]
  frame $rc_win.selfr -relief sunken -bd 4 -bg $calc_bgcol

  text $rc_win.selfr.sel -bd 0 -highlightthickness 0 -insertofftime 0 \
    -bg $entry_bgcol -selectbackground $sel_bgcol -selectforeground $sel_fgcol \
    -selectborderwidth 0 -exportselection yes -height 5 -width 25 -wrap word
   $rc_win.selfr.sel insert end {resid 5 to 85}

  checkbutton $rc_win.bb -highlightthickness 0 \
    -activebackground $calc_bgcol -bg $calc_bgcol \
    -fg $calc_fgcol -activeforeground $calc_fgcol \
    -text {Backbone} -variable [namespace current]::bb_only \
    -command {::RMSDTT::ctrlbb bb}

  checkbutton $rc_win.tr -highlightthickness 0 \
    -activebackground $calc_bgcol -bg $calc_bgcol \
    -fg $calc_fgcol -activeforeground $calc_fgcol \
    -text {Trace} -variable [namespace current]::trace_only \
   -command {::RMSDTT::ctrlbb trace}

  checkbutton $rc_win.noh -highlightthickness 0 \
    -activebackground $calc_bgcol -bg $calc_bgcol \
    -fg $calc_fgcol -activeforeground $calc_fgcol \
    -text {noh} -variable [namespace current]::noh \
   -command {::RMSDTT::ctrlbb noh}

  menubutton  $rc_win.selectionhistory \
        -menu $rc_win.selectionhistory.m -padx 5 -pady 4 \
        -text {History} -relief raised \
        -direction flush


   menu $rc_win.selectionhistory.m

  # Button part.
  frame $calc_top.right -bg $calc_bgcol
  frame $calc_top.right.pushfr -relief ridge -bd 4 -bg $calc_bgcol

  button $calc_top.right.rmsd -relief raised -bd 4 -highlightthickness 0 -text {RMSD} \
    -activebackground $but_abgcol -bg $but_bgcol \
    -fg $calc_fgcol -activeforeground $calc_fgcol \
    -command {::RMSDTT::doRmsd}

  
   button $calc_top.right.align -relief raised -bd 4 -highlightthickness 0 -text {Align} \
    -activebackground $act_bgcol -bg $but_bgcol \
    -activeforeground $act_fgcol -fg $ftr_fgcol  \
    -command {::RMSDTT::doAlign}
  
  frame $calc_top.right.switch -bg $calc_bgcol -relief ridge -bd 4

  radiobutton $calc_top.right.switch.0 -highlightthickness 0 \
    -activebackground $calc_bgcol -bg $calc_bgcol \
    -fg $calc_fgcol -activeforeground $calc_fgcol \
    -text {Top} -variable [namespace current]::rmsd_base -value "top" \
    -command ::RMSDTT::ctrlgui

  radiobutton $calc_top.right.switch.1 -highlightthickness 0 \
    -activebackground $calc_bgcol -bg $calc_bgcol \
    -fg $calc_fgcol -activeforeground $calc_fgcol \
    -text {Average} -variable [namespace current]::rmsd_base -value "ave" \
    -command ::RMSDTT::ctrlgui

  radiobutton $calc_top.right.switch.2 -highlightthickness 0 \
    -activebackground $calc_bgcol -bg $calc_bgcol \
    -fg $calc_fgcol -activeforeground $calc_fgcol \
    -text {Selected} -variable [namespace current]::rmsd_base -value "selected" \
    -command ::RMSDTT::ctrlgui

  # Trajectory part.
  frame $calc_top.right.traj -bg $calc_bgcol -relief ridge -bd 4

  frame $calc_top.right.traj.frames -bg $calc_bgcol -relief ridge -bd 0

  checkbutton $calc_top.right.traj.frames.0 -highlightthickness 0 \
    -activebackground $calc_bgcol -bg $calc_bgcol \
    -fg $calc_fgcol -activeforeground $calc_fgcol \
    -text "Trajectory" -variable [namespace current]::frames_sw \
    -command ::RMSDTT::ctrlgui

  label $calc_top.right.traj.frames.reflabel -text "Frame ref:" \
    -bg $calc_bgcol -fg $calc_fgcol \
    -padx 3 -pady 3

  entry $calc_top.right.traj.frames.ref -bd 0 -highlightthickness 0 -insertofftime 0 \
    -bg $entry_bgcol -selectbackground $sel_bgcol -selectforeground $sel_fgcol \
    -selectborderwidth 0 -exportselection yes -width 5 \
    -textvariable [namespace current]::frames_ref

  checkbutton $calc_top.right.traj.frames.all -highlightthickness 0 \
    -activebackground $calc_bgcol -bg $calc_bgcol \
    -fg $calc_fgcol -activeforeground $calc_fgcol \
    -text "All" -variable [namespace current]::frames_all \
    -command ::RMSDTT::ctrlgui
  
  frame $calc_top.right.traj.skip  -bg $calc_bgcol -relief ridge -bd 0
  
  checkbutton $calc_top.right.traj.skip.0 -highlightthickness 0 \
    -activebackground $calc_bgcol -bg $calc_bgcol \
    -fg $calc_fgcol -activeforeground $calc_fgcol \
    -text "Skip" -variable [namespace current]::skip_sw \
    -command ::RMSDTT::ctrlgui

  label $calc_top.right.traj.skip.inilabel -text "Ini:" \
    -bg $calc_bgcol -fg $calc_fgcol \
    -padx 3 -pady 3

  entry $calc_top.right.traj.skip.ini -bd 0 -highlightthickness 0 -insertofftime 0 \
    -bg $entry_bgcol -selectbackground $sel_bgcol -selectforeground $sel_fgcol \
    -selectborderwidth 0 -exportselection yes -width 3 \
    -textvariable [namespace current]::skip_ini

  label $calc_top.right.traj.skip.stepslabel -text "Steps:" \
    -bg $calc_bgcol -fg $calc_fgcol \
    -padx 3 -pady 3

  entry $calc_top.right.traj.skip.steps -bd 0 -highlightthickness 0 -insertofftime 0 \
    -bg $entry_bgcol -selectbackground $sel_bgcol -selectforeground $sel_fgcol \
    -selectborderwidth 0 -exportselection yes -width 3 \
    -textvariable [namespace current]::skip_steps

  frame $calc_top.right.traj.time -bg $calc_bgcol -relief ridge -bd 0

  checkbutton $calc_top.right.traj.time.0 -highlightthickness 0 \
    -activebackground $calc_bgcol -bg $calc_bgcol \
    -fg $calc_fgcol -activeforeground $calc_fgcol \
    -text "Time (ps)" -variable [namespace current]::time_sw \
    -command ::RMSDTT::ctrlgui

  label $calc_top.right.traj.time.inilabel -text "Ini:" \
    -bg $calc_bgcol -fg $calc_fgcol \
    -padx 3 -pady 3

  entry $calc_top.right.traj.time.inival -bd 0 -highlightthickness 0 -insertofftime 0 \
    -bg $entry_bgcol -selectbackground $sel_bgcol -selectforeground $sel_fgcol \
    -selectborderwidth 0 -exportselection yes -width 6 \
    -textvariable [namespace current]::time_ini

  label $calc_top.right.traj.time.steplabel -text "Step:" \
    -bg $calc_bgcol -fg $calc_fgcol \
    -padx 3 -pady 3

  entry $calc_top.right.traj.time.stepval -bd 0 -highlightthickness 0 -insertofftime 0 \
    -bg $entry_bgcol -selectbackground $sel_bgcol -selectforeground $sel_fgcol \
    -selectborderwidth 0 -exportselection yes -width 6 \
    -textvariable [namespace current]::time_step

  frame $calc_top.right.traj.file -bg $calc_bgcol -relief ridge -bd 0

  checkbutton $calc_top.right.traj.file.plot -highlightthickness 0 \
    -activebackground $calc_bgcol -bg $calc_bgcol \
    -fg $calc_fgcol -activeforeground $calc_fgcol \
    -text "Plot" -variable [namespace current]::plot_sw

  checkbutton $calc_top.right.traj.file.0 -highlightthickness 0 \
    -activebackground $calc_bgcol -bg $calc_bgcol \
    -fg $calc_fgcol -activeforeground $calc_fgcol \
    -text "Save to file:" -variable [namespace current]::file_out_sw \
    -command [namespace code {
      if {$file_out_sw} {
	$w.top.right.traj.file.name config -state normal
      } else {
	$w.top.right.traj.file.name config -state disable
      }
    }]

  entry $calc_top.right.traj.file.name -bd 0 -highlightthickness 0 -insertofftime 0 \
    -bg $entry_bgcol -selectbackground $sel_bgcol -selectforeground $sel_fgcol \
    -selectborderwidth 0 -exportselection yes -width 15 \
    -textvariable [namespace current]::file_out -state disable

  # Pack the top part widgets.
  pack $calc_top -side top -fill both -expand 1 

  ### ...selection...
  pack $calc_top.left -side left -fill both -expand 1
  pack $rc_win -side left -fill both -expand 1
  pack $rc_win.selfr -side top -fill both -expand 1 -padx 4 -pady 4
  pack $rc_win.bb -side left 
  pack $rc_win.tr -side left 
  pack $rc_win.noh -side left 
  pack $rc_win.selectionhistory -side right
  pack $rc_win.selfr.sel -side top -fill both -expand 1
  pack $rc_win.selfr.sel -side top -fill both -expand 1
  
  ### ...button...
  pack $calc_top.right -side left -fill both
  pack $calc_top.right.pushfr -side top -fill x -expand 1
  pack $calc_top.right.rmsd $calc_top.right.align -side left -fill x -expand 1 -in $calc_top.right.pushfr
  pack $calc_top.right.switch -side top -fill both
  #pack $calc_top.right.switch.from $calc_top.right.switch.frame -side top -fill both
  #pack $calc_top.right.switch.frame -side top -fill both
  pack $calc_top.right.switch.0 $calc_top.right.switch.1 $calc_top.right.switch.2\
    -side left -fill both -padx 2 -pady 2
  
  ### ...trajectory...
  pack $calc_top.right.traj -side top -fill both
  pack $calc_top.right.traj.frames -side top -fill both
  pack $calc_top.right.traj.frames.0 $calc_top.right.traj.frames.reflabel $calc_top.right.traj.frames.ref $calc_top.right.traj.frames.all\
    -side left -fill both -padx 2 -pady 2
  pack $calc_top.right.traj.skip -side top -fill both
  pack $calc_top.right.traj.skip.0 $calc_top.right.traj.skip.inilabel $calc_top.right.traj.skip.ini\
    $calc_top.right.traj.skip.stepslabel $calc_top.right.traj.skip.steps\
    -side left -fill both -padx 2 -pady 2
  pack $calc_top.right.traj.time -side top -fill both
  pack $calc_top.right.traj.time.0 $calc_top.right.traj.time.inilabel $calc_top.right.traj.time.inival\
    $calc_top.right.traj.time.steplabel $calc_top.right.traj.time.stepval\
    -side left -fill both -padx 2 -pady 2
  pack $calc_top.right.traj.file -side top -fill both
  pack $calc_top.right.traj.file.plot $calc_top.right.traj.file.0 $calc_top.right.traj.file.name\
    -side left -fill both -padx 2 -pady 2


  # Create the bottom part, displaying the RMSDs of
  # all the molecules from the average structure.
  set calc_mid [frame $w.middle -relief ridge -bd 4 -bg $calc_bgcol]

frame $calc_mid.titlebar -bg $calc_bgcol -relief raised -bd 2
frame $calc_mid.titlebar.left -bg $calc_bgcol

grid rowconfigure    $calc_mid.titlebar	0 -pad 0  -weight 100
grid columnconfigure $calc_mid.titlebar 0 -pad 0  -weight 100
grid $calc_mid.titlebar.left -in $calc_mid.titlebar -column 0 -row 0 -columnspan 1 -rowspan 1 -ipadx 0 -ipady 0 -padx 0 -pady 0 -sticky nesw

frame $calc_mid.titlebar.right -bg $calc_bgcol
grid rowconfigure    $calc_mid.titlebar 0	 -pad 0  -weight 100
grid columnconfigure $calc_mid.titlebar 1  -pad 0  -weight 100
grid $calc_mid.titlebar.right -in $calc_mid.titlebar -column 1 -row 0 -columnspan 1 -rowspan 1 -ipadx 0 -ipady 0 -padx 0 -pady 0 -sticky nesw

 label $calc_mid.titlebar.left.label -text {ID   Molecule}  \
    -bg $calc_bgcol -fg $calc_fgcol \
    -padx 3 -pady 3
    
    label $calc_mid.titlebar.right.label -text {RMSD Average:} \
    -bg $calc_bgcol -fg $calc_fgcol \
    -padx 3 -pady 3
    
  set rms_disp [frame $w.middle.body -bg $calc_bgcol]
  frame $rms_disp.left -bg $calc_bgcol -bd 0
  frame $rms_disp.right -bg $calc_bgcol -bd 0

  set pdb_list [listbox $rms_disp.pdb_names -relief raised -bd 2 -height 10 \
      -bg $calc_bgcol -fg $calc_fgcol \
      -highlightthickness 0 -highlightbackground $calc_bgcol \
      -yscrollcommand [namespace code {$rms_disp.scrbar set}] \
			]

  set rms_list [listbox $rms_disp.rms_values -relief raised -bd 2 -height 10 \
      -bg $calc_bgcol -fg $calc_fgcol  \
      -highlightthickness 0 -highlightbackground $calc_bgcol \
      -yscrollcommand [namespace code {$rms_disp.scrbar set}] \
			]

  set all_mols [molactive]
  foreach i $all_mols {
    $pdb_list insert end [format "%-4s%-10s" [molinfo index $i] [molinfo $i get name]]
    set rms_ave($i) {}
    $rms_list insert end $rms_ave($i)
  }

  set rmstot_lbl [label $rms_disp.rmstot_lbl -text {Overall Average RMSD:} -anchor w -pady 3 \
      -bg $calc_bgcol -relief raised -bd 2 -fg $calc_fgcol]
  set rmstot_val [label $rms_disp.rmstot_val -textvariable [namespace current]::tot_rms -anchor w -pady 3 \
      -bg $calc_bgcol -relief raised -bd 2 -fg $calc_fgcol]


  scrollbar $rms_disp.scrbar \
    -relief raised -activerelief raised -bd 2 -elementborderwidth 2 \
    -bg $calc_bgcol -troughcolor $scr_trough -highlightthickness 0 \
    -activebackground $but_abgcol \
    -orient vert \
    -command  {RMSDTT::two_scroll}
   

  # Pack the bottom part widgets.
  pack $calc_mid.titlebar.left.label  -side left -in  $calc_mid.titlebar.left
  pack $calc_mid.titlebar.right.label -side left -in  $calc_mid.titlebar.right
  
  pack $calc_mid.titlebar -side top -fill x -in  $calc_mid
  pack $calc_mid -side top -fill both -expand 1
  pack $rms_disp -side top -expand 1 -fill both
  pack $rms_disp.left $rms_disp.right -side left -fill both -expand 1
  pack $rms_disp.scrbar -side left -fill y
  pack $pdb_list -side top -fill both -expand 1 -in $rms_disp.left
  pack $rmstot_lbl -side top -fill x -in $rms_disp.left

  pack $rms_list -side top -fill both -expand 1 -in $rms_disp.right
  pack $rmstot_val -side bottom -fill x -in $rms_disp.right

  frame $w.bottom -bg white -height 100
  button $w.bottom.erase -relief raised -bd 4 -highlightthickness 0 -text {Erase from list} \
    -command [namespace code {
    	set indexnum [$pdb_list index active] 
    	$pdb_list delete active 
    	$rms_list delete $indexnum 
    	}] \
    -activebackground $but_abgcol -bg $but_bgcol 
   
  button $w.bottom.allinmem -relief raised -bd 4 -highlightthickness 0 -text {All in memory} \
    -activebackground $but_abgcol -bg $but_bgcol \
    -command [namespace code {
        $pdb_list delete 0 end
        for {set i 0} {$i<[molinfo num]} {incr i} {
          set molid [molinfo index $i]
          $pdb_list insert end [format "%-4s%-10s" [molinfo index $i] [molinfo $molid get name]]      
        }
      } ]

  button $w.bottom.onlyactive -relief raised -bd 4 -highlightthickness 0 -text {Only active} \
    -activebackground $but_abgcol -bg $but_bgcol \
    -command {RMSDTT::onlyactive}

  menubutton $w.bottom.assemblymenu \
        -menu $w.bottom.assemblymenu.m -padx 5 -pady 4 \
        -text {Assembly} -relief raised 
    
  pack $w.bottom.erase $w.bottom.allinmem $w.bottom.onlyactive $w.bottom.assemblymenu -in $w.bottom  -side left -fill x -expand 1
  pack $w.bottom -in $w -side bottom -fill x -expand 1    


  menu $w.bottom.assemblymenu.m -postcommand [namespace code {
      set menu $w.bottom.assemblymenu.m
      $menu delete 0 end
      if {[array exist assembly]} {
        set nameofassemblies [array names assembly]
        for {set i 0} {$i<[llength $nameofassemblies]} {incr i} {
          set nameofcurrentassembly [lindex $nameofassemblies $i]

          $menu add  radiobutton -label $nameofcurrentassembly -variable [namespace current]::selectedassembly -value $nameofcurrentassembly \
            -command [namespace code {
                $pdb_list delete 0 end
                $rms_list delete 0 end
          
                foreach pdb $assembly($selectedassembly) {
                  $pdb_list insert end $pdb
                }
              }] 
        }
      }
    } ]


# update the History menu
 ListHisotryPullDownMenu
 ctrlgui
}

