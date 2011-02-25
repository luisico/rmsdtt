#
#             RMSD Trajectory Tool v3.0
#
# Iterative Fitting by Residue addon to RMSD Trajectory Tool
#
# http://physiology.med.cornell.edu/faculty/hweinstein/vmdplugins/rmsdtt

# Author
# ------
#      Luis Gracia, PhD
#      lug2002@med.cornell.edu

#      In collaboration with Joshua A. Speidel

#      Weill Cornell Medical College
#      New York, NY

# Documentation
# -------------
# See the index.html file distributed with this file. For more update documentation see
# http://physiology.med.cornell.edu/faculty/hweinstein/vmdplugins/rmsdtt



proc rmsdtt::iterativeFitGUI {} {
  variable w
  
  variable ifit_niter      3
  variable ifit_type       expmin
  variable ifit_factor     1.0
  variable ifit_convsw     1
  variable ifit_conv       0.01
  variable ifit_plot       0
  variable ifit_scale      RWB
  variable ifit_update     0
  variable ifit_repre      0
  variable ifit_sel2       all
  variable ifit_style      NewRibbons
  variable ifit_frames_all 1
  variable ifit_replace    1
  variable ifit_save       1
  variable ifit_file       ifitbyres.dat
  variable ifit_cluster      0
  variable ifit_cluster_only 0
  variable ifit_cluster_fit  0
  variable ifit_cluster_weights_from 0
  
  labelframe $w.ifit -text "Iterative Fitting by residue" -relief ridge -bd 2
  pack $w.ifit -side top -fill x
  
  frame $w.ifit.row1
  pack $w.ifit.row1 -side top -fill x
  
  button $w.ifit.row1.button -relief raised -bd 2 -text "Fit" -command [namespace current]::iterativeFit
  label $w.ifit.row1.niterlabel -text "Iters:"
  entry $w.ifit.row1.niter -width 4 -textvariable [namespace current]::ifit_niter
  label $w.ifit.row1.factorlabel -text "Factor:"
  menubutton $w.ifit.row1.type -textvariable [namespace current]::ifit_type -menu $w.ifit.row1.type.menu -relief raised -direction flush -width 6
  menu $w.ifit.row1.type.menu -tearoff no
  foreach type [list linear exp expmin minmax gaussian] {
    $w.ifit.row1.type.menu add radiobutton -label $type -variable [namespace current]::ifit_type -value $type
  }
  entry $w.ifit.row1.factor -width 4 -textvariable [namespace current]::ifit_factor
  checkbutton $w.ifit.row1.convsw -text "Convergence:" -variable [namespace current]::ifit_convsw
  entry $w.ifit.row1.conv -width 5 -textvariable [namespace current]::ifit_conv
  
  pack $w.ifit.row1.button $w.ifit.row1.niterlabel $w.ifit.row1.niter $w.ifit.row1.factorlabel $w.ifit.row1.type $w.ifit.row1.factor $w.ifit.row1.convsw $w.ifit.row1.conv -side left -anchor w
  
  frame $w.ifit.row2
  pack $w.ifit.row2 -side top -fill x
  
  checkbutton $w.ifit.row2.plot -text "Plot" -variable [namespace current]::ifit_plot
  checkbutton $w.ifit.row2.save -text "Save" -variable [namespace current]::ifit_save
  entry $w.ifit.row2.file -width 20 -textvariable [namespace current]::ifit_file
  
  pack $w.ifit.row2.plot $w.ifit.row2.save $w.ifit.row2.file -side left -anchor w
  
  frame $w.ifit.row3
  pack $w.ifit.row3 -side top -fill x
  
  checkbutton $w.ifit.row3.repre -text "Rep" -variable [namespace current]::ifit_repre
  entry $w.ifit.row3.sel2 -width 10 -textvariable [namespace current]::ifit_sel2
  checkbutton $w.ifit.row3.replace -text "replace" -variable [namespace current]::ifit_replace
  menubutton $w.ifit.row3.style -text "Style" -menu $w.ifit.row3.style.menu -relief raised -direction flush
  menu $w.ifit.row3.style.menu -tearoff no
  foreach style [list Lines Bonds DynamicBonds HBonds Points VDW CPK Licorice Trace Tube Ribbons NewRibbons Cartoon NewCartoon MSMS Surf VolumeSlice Isosurface Dotted Solvent] {
    $w.ifit.row3.style.menu add radiobutton -label $style -variable [namespace current]::ifit_style -value $style
  }
  menubutton $w.ifit.row3.scale -text "Scale" -menu $w.ifit.row3.scale.menu -relief raised -direction flush
  menu $w.ifit.row3.scale.menu -tearoff no
  foreach item [list unchanged RGB BGR RWB BWR RWG GWR GWB BWG BlkW WBlk] {
    $w.ifit.row3.scale.menu add radiobutton -label $item -variable [namespace current]::ifit_scale -value $item
  }
  checkbutton $w.ifit.row3.framesall -text "all frames" -variable [namespace current]::ifit_frames_all
  checkbutton $w.ifit.row3.update -text "Update" -variable [namespace current]::ifit_update
  
  pack $w.ifit.row3.repre $w.ifit.row3.sel2 $w.ifit.row3.replace $w.ifit.row3.style $w.ifit.row3.scale $w.ifit.row3.framesall $w.ifit.row3.update -side left -anchor w
  
  frame $w.ifit.row4
  pack $w.ifit.row4 -side top -fill x
  
  checkbutton $w.ifit.row4.cluster -text "cluster" -variable [namespace current]::ifit_cluster
  checkbutton $w.ifit.row4.cluster_only -text "cluster only" -variable [namespace current]::ifit_cluster_only
  checkbutton $w.ifit.row4.cluster_fit -text "cluster fit" -variable [namespace current]::ifit_cluster_fit
  label $w.ifit.row4.cluster_weights_l -text "weights from:"
  entry $w.ifit.row4.cluster_weights -width 2 -textvariable [namespace current]::ifit_cluster_weights_from
  pack $w.ifit.row4.cluster $w.ifit.row4.cluster_only $w.ifit.row4.cluster_fit $w.ifit.row4.cluster_weights_l $w.ifit.row4.cluster_weights -side left -anchor w
}


proc rmsdtt::iterativeFit {} {
  variable w
  variable datalist
  variable ifit_niter
  variable ifit_type
  variable ifit_factor
  variable ifit_convsw
  variable ifit_conv
  variable ifit_plot
  variable ifit_scale
  variable ifit_update
  variable ifit_repre
  variable ifit_sel2
  variable ifit_style
  variable ifit_frames_all
  variable ifit_replace
  variable ifit_rep
  variable ifit_save
  variable ifit_file
  variable ifit_cluster
  variable ifit_cluster_only
  variable ifit_cluster_fit
  variable ifit_cluster_weights_from

  set sel1 [set_sel]
  if {$sel1 == ""} {
    showMessage "Selection is empty selection!"
    return -code return
  }
  set sel2 $ifit_sel2
  #puts $sel1
  #puts $sel2
  
  set target_mol [$datalist(id) get 0 end]
  set nmols [llength $target_mol]

  # Check number of atoms
  set message ""
  for {set i 0} {$i < $nmols} {incr i} {
    set mol1 [lindex $target_mol $i]
    for {set j [expr {$i+1}]} {$j < $nmols} {incr j} {
      set mol2 [lindex $target_mol $j]
      if {[[atomselect $mol1 $sel1] num] != [[atomselect $mol2 $sel1] num]} {
	append message "$mol1 ([[atomselect $mol1 $sel1] num])\t\t$mol2 ([[atomselect $mol2 $sel1] num])\n"
      }
    }
  }
  if {$message != ""} {
    set message "Number of atoms selected differ for molecules:\n$message"
    showMessage $message
    return -code return
  }
  

  set divide 0
  if {$divide > 0 && $nmols > 1} {
    showMessage "divide by frames can only be used with one molecule\n"
    return -code return
  }  

  set fast 1
  set nsteps 0

  # Initialize objects and weights
  foreach mol $target_mol  {
    set nframes($mol) [molinfo $mol get numframes]

    # Make objects for molecule selections
    set sel_ref($mol) [atomselect $mol $sel1]
    set sel_current($mol) [atomselect $mol $sel1]
    set sel_move($mol) [atomselect $mol "all"]
    
    # Count number of comparison to be done per iteration
    foreach mol2 $target_mol {
      if {$mol == $mol2} {
	for {set i 1} {$i < $nframes($mol)} {incr i} {
	  set nsteps [expr {$nsteps + $i}]
	}
      } elseif {$mol < $mol2} {
	set nsteps [expr {$nsteps + $nframes($mol)*[molinfo $mol2 get numframes]}]
      }
    }

    # Make objects for residue sel
    set residues($mol) [lsort -unique -integer [[atomselect $mol $sel1] get residue]]
    for {set i 0} {$i < [llength $residues($mol)]} {incr i} {
      set sel_res_ref($mol:$i) [atomselect $mol "residue [lindex $residues($mol) $i] and $sel1"]
      set sel_res($mol:$i)     [atomselect $mol "residue [lindex $residues($mol) $i] and $sel1"]
      if {[$sel_res_ref($mol:$i) num] > 1} {
	set fast 0
      }
      #puts "$mol - [lindex $residues($mol) $i] - [$sel_res_ref($mol:$i) get index]"
    }

    # Set the user field of all frames = 1 for use in initial weighting scheme
    if {!$ifit_cluster_only} {
      set sel2_atoms [atomselect $mol $sel2]
      for {set i 0} {$i < $nframes($mol)} {incr i} {
	$sel2_atoms frame $i
	$sel2_atoms set user 1
      }
    }

    # Division of frames in contigous windows
    if {$divide > 0} {
      set divide_frames 0
      for {set i $divide} {$i < $nframes($mol)} {set i [expr {$i + $divide}]} {
	lappend divide_frames $i
      }
      puts $divide_frames
      return
    }
  }

  set signalstep [expr {$nsteps/double(20)}]
  puts "Number of steps per iteration: $nsteps"
  puts "Signal every $signalstep steps"

  if {$fast} {
    puts "Testing for fast algorithm..."
    if [catch { set ret [measure rmsd $sel_ref([lindex $target_mol 0]) $sel_ref([lindex $target_mol 0]) byatom] } msg] {
      puts "   :-( Hacked VMD not found! Please contact the RMSDTT developer"
      set fast 0
    } else {
      puts "   :-) Hacked VMD found"
    }

    if [catch {package require BLT} msg] {
      puts "   :-( Package BLT not found. Please install BLT for you Tcl/Tk version first"
      set fast 0
    } else {
      puts "   :-) Package BLT found"
    }

    if {$fast} {
      puts "   Nice, you got it all to run at light speed!"
      
      # Delete residue objects
      foreach mol $target_mol  {
	for {set i 0} {$i < [llength $residues($mol)]} {incr i} {
	  $sel_res_ref($mol:$i) delete
	  $sel_res($mol:$i) delete
	}
      }
    } else {
      set ifit_convsw 0
      puts "Switching Convergence off."
    }
  } else {
    set ifit_convsw 0
    puts "Switching Convergence off."
  }

  # Check number of residues
  set message ""
  for {set i 0} {$i < $nmols} {incr i} {
    set mol1 [lindex $target_mol $i]
    for {set j [expr {$i+1}]} {$j < $nmols} {incr j} {
      set mol2 [lindex $target_mol $j]
      if {[llength $residues($mol1)] != [llength $residues($mol1)]} {
	append message "$mol1 ([llength $residues($mol1)])\t\t$mol2 ([llength $residues($mol1)])\n"
      }
    }
  }
  if {$message != ""} {
    set message "Number of residues selected differ for molecules:\n$message"
    showMessage $message
    return -code return
  }
  set nresidues [llength $residues([lindex $target_mol 0])]

  # Initialize plot
  set plot_use 0
  if {$ifit_plot} {
    if [catch {package require multiplot} msg] {
      showMessage "Plotting in Multiplot not available: package multiplot not installed!\nDo you have the latest VMD version?"
    } else {
      set plot_use 1
      set title "Rmsd vs Residue"
      set xlab "Residue"
      set ylab "Rmsd (A)"
      set plothandle [multiplot -title $title -xlabel $xlab -ylabel $ylab -nostats]
    }
  }
  
  # Header for file
  if {$ifit_save && !$ifit_cluster_only} {
    set fid [open $ifit_file w]
    puts $fid [format "%4s %7s %5s %4s %5s %7s %5s" "iter" "residue" "resid" "name" "chain" "mean" "w"]
  }

  # Create representation
  if {$ifit_repre && !$ifit_cluster_only} {
    foreach mol $target_mol {
      mol rep $ifit_style
      mol color User
      mol selection $sel2
      set mol [expr {$mol+0}]
      set add 1
      if {$ifit_replace} {
	if {[info exists ifit_rep($mol)]} {
	  for {set i 0} {$i < [molinfo $mol get numreps]} { incr i} {
	    #puts "$i [mol repname $mol $i]"
	    if {$ifit_rep($mol) eq [mol repname $mol $i]} {
	      mol modrep [mol repindex $mol $ifit_rep($mol)] $mol
	      set add 0
	      break
	    }
	  }
	}
      }
      if {$add} {
	mol addrep $mol
	set ifit_rep($mol) [mol repname $mol [expr {[molinfo $mol get numreps]-1} ]]
      }
      if {$ifit_frames_all} {
	mol drawframes $mol [mol repindex $mol $ifit_rep($mol)] 0:[molinfo $mol get numframes]
      } else {
	mol drawframes $mol [mol repindex $mol $ifit_rep($mol)] now
      }
      #      mol modstyle [mol repindex $mol $ifit_rep($mol)] $mol $ifit_style
      #      mol modcolor [mol repindex $mol $ifit_rep($mol)] $mol User
      #      mol modselect [mol repindex $mol $ifit_rep($mol)] $mol $sel2
    }
    
    # Change scale
    if {$ifit_scale != "unchanged"} {
      color scale method $ifit_scale
    }
  }

  if {$ifit_update} {display update}


  # Initizalize rmsd by residue and weights
  if {$fast} {
    ::blt::vector create zeros($nresidues)
    ::blt::vector create weights($nresidues)
    weights expr {zeros + 1.0}
    zeros dup temp
  } else {
    for {set res 0} {$res < $nresidues} {incr res} {
      lappend rmsd_mean 0.0
    }
  }

  # Iterate over fitting and weighting niter times
  if {!$ifit_cluster_only} {
    set iter 1
    if {$ifit_convsw && $fast} {
	::blt::vector create conv_vector
	conv_vector expr {zeros + 1.0}
    }
    while {$iter <= $ifit_niter} {
      puts -nonewline "Iteration: $iter "
      
      # Reset rmsd by residue
      if {$fast} {
	zeros dup rmsd_mean
      } else {
	for {set res 0} {$res < $nresidues} {incr res} {
	  lset rmsd_mean $res 0.0
	}
      }
      
      set count 0
      set count2 0
      for {set i 0} {$i < $nmols} {incr i} {
	set mol1 [lindex $target_mol $i]
	for {set j 0} {$j < $nframes($mol1)} {incr j} {
	  $sel_ref($mol1) frame $j
	  set sel_reference $sel_ref($mol1)
	  if {!$fast} {
	    for {set res 0} {$res < $nresidues} {incr res} {
	      $sel_res_ref($mol1:$res) frame $j
	    }
	  }
	  
	  for {set k 0} {$k < $nmols} {incr k} {
	    set mol2 [lindex $target_mol $k]
	    for {set l 0} {$l < $nframes($mol2)} {incr l} {
	      if {$i == $k && $j >= $l   ||   $i > $k} {continue}
	      
	      incr count
	      incr count2
	      if {$count2 > $signalstep} {
		puts -nonewline "."
		set count2 0
	      }

	      $sel_current($mol2) frame $l
	      $sel_move($mol2) frame $l
	      
	      if {$fast} {
		$sel_move($mol2) move [measure fit $sel_current($mol2) $sel_reference weight [weights range 0 end]]
		lassign [measure rmsd $sel_reference $sel_current($mol2) byatom] global_rmsd ifit_rmsd
		#puts "$j $l $ifit_rmsd"
		temp set $ifit_rmsd
		rmsd_mean set [rmsd_mean + temp]
		
	      } else {
		$sel_move($mol2) move [measure fit $sel_current($mol2) $sel_reference weight user]
		for {set res 0} {$res < $nresidues} {incr res} {
		  $sel_res($mol2:$res) frame $l
		  lset rmsd_mean $res [expr {[lindex $rmsd_mean $res] + [measure rmsd $sel_res_ref($mol1:$res) $sel_res($mol2:$res)]}]
		}
	      }
	    }
	  }
	}
      }
      
      if {$fast} {
	# Compute mean, mix and max
	rmsd_mean set [rmsd_mean / $count]
	set rmsd_min [set [namespace current]::rmsd_mean(min)]
	set rmsd_max [set [namespace current]::rmsd_mean(max)]
	
	# Compute weights
	switch $ifit_type {
          linear {
            weights expr { $ifit_factor * rmsd_mean }
          }
	  exp {
	    weights expr { exp(-$ifit_factor * rmsd_mean) }
	  }
	  expmin {
	    weights expr { exp(-$ifit_factor * ( rmsd_mean - $rmsd_min)) }
	  }
	  minmax {
	    if {$rmsd_max == $rmsd_min} {#!!!!!!!!!!!!!!!!
	      set weight 1
	    } else { 
	      weights expr { ($rmsd_max - rmsd_mean) / ($rmsd_max - $rmsd_min) }
	    }
	  }
	  gaussian {
	    weights expr { exp(-(rmsd_mean * rmsd_mean) / $ifit_factor) }
	  }
	}
	
	# Update display
	if {$ifit_update} {
	  foreach mol $target_mol  {
	    for {set i 0} {$i < [molinfo $mol get numframes]} {incr i} {
	      $sel_ref($mol) frame $i
	      $sel_ref($mol) set user [weights range 0 end]
	      $sel_ref($mol) set beta [rmsd_mean range 0 end]
	    }
	  }
	}
	
	# Plot
	if {$plot_use} {
	  set color [index2rgb $iter]
	  set legend "Iter $iter"
	  $plothandle add $residues([lindex $target_mol 0]) [rmsd_mean range 0 end] -marker point -radius 2 -fillcolor $color -linecolor $color -nostats -legend $legend
	}
	
	# Save
	if {$ifit_save} {
	  set data [$sel_ref([lindex $target_mol 0]) get {residue resid resname chain}]
	  for {set res 0} {$res < $nresidues} {incr res} {
	    lassign [lindex $data $res] d_residue d_resid d_resname d_chain
	    puts $fid [format "%4s %7d %5d %4s %5s %7.3f %5.3f" $iter $d_residue $d_resid $d_resname $d_chain [rmsd_mean index $res] [weights index $res]]
	  }
	}

	# Convergence
	if {$ifit_convsw} {
	  if {$iter == 1} {
	    conv_vector expr { weights }
	    set t [::blt::vector expr { sqrt(sum( weights * weights )) }]
	    set conv 0
	  } else {
	    set t [::blt::vector expr { sqrt(sum( (weights - conv_vector) * (weights - conv_vector) )) }]
	  }
	  set conv [expr abs($conv - $t)]
	  puts -nonewline [format " %8.5f" $conv]
	  if {$conv < $ifit_conv} {
	    puts "\nConvergence achieved"
	    break
	  }
	  set conv $t
	}
	
      } else { # not fast
	if {$plot_use} {
	  set y {}
	  set x {}
	}
      
	# Compute mean, mix and max
	set rmsd_min [expr {[lindex $rmsd_mean 0] / $count}]
	set rmsd_max $rmsd_min
	for {set res 0} {$res < [llength $rmsd_mean]} {incr res} {
	  lset rmsd_mean $res [expr {[lindex $rmsd_mean $res] / $count}]
	  if {[lindex $rmsd_mean $res] < $rmsd_min} {
	    set rmsd_min [lindex $rmsd_mean $res]
	    continue
	  }
	  if {[lindex $rmsd_mean $res] > $rmsd_max} {
	    set rmsd_max [lindex $rmsd_mean $res]
	  }
	}
	
	# Compute weights
	for {set res 0} {$res < $nresidues} {incr res} {
	  set r [lindex $rmsd_mean $res]
	  switch $ifit_type {
            linear {
              set weight [expr {$ifit_factor*$r}]
            }
	    exp {
	      set weight [expr {exp(-$ifit_factor*$r)}]
	    }
	    expmin {
	      set weight [expr {exp(-$ifit_factor*($r - $rmsd_min))}]
	    }
	    minmax {
	      if {$rmsd_max == $rmsd_min} {
		set weight 1
	      } else { 
		set weight [expr {($rmsd_max-$r) / ($rmsd_max-$rmsd_min)}]
	      }
	    }
	    gaussian {
	      set weight [expr {exp(-($r*$r)/$ifit_factor)}]
	    }
	  }
	  #puts [format "\tRes: %5d %6.3f %6.3f %6.3f %6.2f" $res $rmsd_mean($res) $rmsd_min $rmsd_max $weight]
	  
	  # Update
	  foreach mol $target_mol  {
	    for {set i 0} {$i < [molinfo $mol get numframes]} {incr i} {
	      $sel_res($mol:$res) frame $i
	      $sel_res($mol:$res) set user $weight
	      $sel_res($mol:$res) set beta $r
	    }
	  }
	  
	  # Plot
	  if {$plot_use} {
	    lappend x $res
	    lappend y $r
	    set color [index2rgb $iter]
	    set legend "Iter $iter"
	    $plothandle add $x $y -marker point -radius 2 -fillcolor $color -linecolor $color -nostats -legend $legend
	  }
	  
	  # Save
	  if {$ifit_save} {
	    set data [lindex [$sel_res([lindex $target_mol 0]:$res) get {residue resid resname chain user}] 0]
	    puts $fid [format "%4s %7d %5d %4s %5s %7.3f %5.3f" $iter [lindex $data 0] [lindex $data 1] [lindex $data 2] [lindex $data 3] $r [lindex $data 4]]
	  }
	}

	# convergence not available if not fast
      }
      
      incr iter
      puts ""
    }

    puts "Done with fitting"
    
    # Update atom properties
    if {$fast} {
      foreach mol $target_mol  {
	for {set i 0} {$i < $nframes($mol)} {incr i} {
	  $sel_ref($mol) frame $i
	  $sel_ref($mol) set user [weights range 0 end]
	  $sel_ref($mol) set beta [rmsd_mean range 0 end]
	}
      }
    }
    if {$ifit_update} {display update}
    
    if {$plot_use} {$plothandle replot}
    
    if {$ifit_save} {close $fid}
    
  }
 
  # clustering
  if {$ifit_cluster} {
    puts "Clustering"
    set fid1 [open "$ifit_file.cluster1" w]
#     set fid2 [open "$ifit_file.cluster2" w]
#     puts -nonewline $fid2 [format "%6s %6s" "str1" "str2"]
#     for {set res 0} {$res < $nresidues} {incr res} {
#       puts -nonewline $fid2 [format " %7d" [lindex $residues([lindex $target_mol 0]) $res]]
#     }
#     puts $fid2 ""
    
    if {$ifit_cluster_only} {
      weights set [[atomselect $ifit_cluster_weights_from $sel1] get user]
    }
    if {$fast} {
      set weight_sum [::blt::vector expr { sum(weights) }]
    }
    
    for {set i 0} {$i < $nmols} {incr i} {
      set mol1 [lindex $target_mol $i]
      for {set j 0} {$j < $nframes($mol1)} {incr j} {
	$sel_ref($mol1) frame $j
	set sel_reference $sel_ref($mol1)
	if {!$fast} {
	  for {set res 0} {$res < $nresidues} {incr res} {
	    $sel_res_ref($mol1:$res) frame $j
	  }
	}
	
	for {set k 0} {$k < $nmols} {incr k} {
	  set mol2 [lindex $target_mol $k]
	  for {set l 0} {$l < $nframes($mol2)} {incr l} {
	    if {$i == $k && $j >= $l   ||   $i > $k} {continue}

	    $sel_current($mol2) frame $l
	    if {$ifit_cluster_fit} {
	      $sel_move($mol2) frame $l
	      if {$fast} {
		$sel_move($mol2) move [measure fit $sel_current($mol2) $sel_reference weight [weights range 0 end]]
	      } else {
		$sel_move($mol2) move [measure fit $sel_current($mol2) $sel_reference weight user]
	      }
	    }

	    set cluster_rms [measure rmsd $sel_current($mol2) $sel_reference]

#	    puts -nonewline $fid2 [format "%2d:%-3d %2d:%-3d" $mol1 $j $mol2 $l]
	    if {$fast} {
	      lassign [measure rmsd $sel_current($mol2) $sel_reference byatom weight [weights range 0 end]] cluster_rmsw cluster_ifit
#	      puts $fid2 $cluster_ifit
	      temp set $cluster_ifit
	      set cluster_ifit [expr {[::blt::vector expr {sum(temp)}] / $weight_sum }]

	    } else {   # not fast
	      set cluster_rmsw [measure rmsd $sel_current($mol2) $sel_reference weight user]
	      set cluster_ifit 0.0
	      set weight_tot 0.0
	      for {set res 0} {$res < $nresidues} {incr res} {
		$sel_res($mol2:$res) frame $l
		set rmsd [measure rmsd $sel_res_ref($mol1:$res) $sel_res($mol2:$res)]
		set weight [lindex [$sel_res_ref($mol1:$res) get user] 0]
		set weight_tot [expr {$weight_tot + $weight}]
		set cluster_ifit [expr {$cluster_ifit + $weight * $rmsd}]
#		puts -nonewline $fid2 [format " %7.3f" $rmsd]
	      }
	      set cluster_ifit [expr {$cluster_ifit / $weight_tot}]
#	      puts $fid2 ""
	    }
	    
	    puts $fid1 [format "%2d:%-3d %2d:%-3d %7.3f %7.3f %7.3f" $mol1 $j $mol2 $l $cluster_rms $cluster_rmsw $cluster_ifit]

	  }
	}
      }
    }
    close $fid1
#    close $fid2

  }

  
  #puts [array get rmsd_mean]
#  return [array get rmsd_mean]

}
