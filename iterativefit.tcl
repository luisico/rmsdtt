#
#             RMSD Trajectory Tool v4.0
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
  variable ifit_convsw     0
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
  checkbutton $w.ifit.row4.cluster_fit -text "cluster fit" -variable [namespace current]::ifit_cluster_fit
  checkbutton $w.ifit.row4.cluster_only -text "cluster only" -variable [namespace current]::ifit_cluster_only
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
  variable ifit_hack
  variable ifit_tensor

  set sel1 [set_sel]
  if {$sel1 == ""} {
    showMessage "Selection is empty selection!"
    return -code return
  }
  set sel2 $ifit_sel2

  set target_mol [$datalist(id) get 0 end]
  set nmols [llength $target_mol]

  # Check number of atoms
  [namespace current]::iterativeFit_check_atoms

  # Test fast algorithm
  [namespace current]::iterativeFit_test_fast

  # Initialize objects and weights
  set nsteps 0
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
    set ifit_natoms 0
    set residues($mol) [lsort -unique -integer [[atomselect $mol $sel1] get residue]]
    for {set i 0} {$i < [llength $residues($mol)]} {incr i} {
      set sel_res_ref($mol:$i) [atomselect $mol "residue [lindex $residues($mol) $i] and $sel1"]
      set sel_res($mol:$i)     [atomselect $mol "residue [lindex $residues($mol) $i] and $sel1"]
      if {!$ifit_natoms && [$sel_res_ref($mol:$i) num] > 1} {
        # More than one atom selected per residue
        puts "Warning: more than one atom selected per residue"
        set ifit_natoms 1
      }
      # puts "$mol - [lindex $residues($mol) $i] - [$sel_res_ref($mol:$i) get index]"
    }

    # Set the user field of all frames = 1 for use in initial weighting scheme
    if {!$ifit_cluster_only} {
      set sel2_atoms [atomselect $mol $sel2]
      for {set i 0} {$i < $nframes($mol)} {incr i} {
        $sel2_atoms frame $i
        $sel2_atoms set user 1
      }
    }
  }

  set ifit_convsw 0
  puts "Switching Convergence off."

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

  # Only clustering
  if {$ifit_cluster_only} {
    set ifit_plot 0
    set ifit_save 0
    set ifit_repre 0
  }

  # Initialize plot
  if {$ifit_plot} {
    set plothandle [[namespace current]::iterativeFit_plot_header]
  }

  # Header for file
  if {$ifit_save} {
    set fid [open $ifit_file w]
    puts -nonewline $fid [[namespace current]::iterativeFit_save_header]
  }

  # Create representation
  if {$ifit_repre} {
    [namespace current]::iterativeFit_representation
  }

  if {$ifit_update} {display update}

  # Initizalize rmsd by residue and weights
  # TODO: another way to initialize lists?
  for {set res 0} {$res < $nresidues} {incr res} {
    lappend rmsd 0.0
    lappend weights 1.0
  }
  if {$ifit_tensor} {
    if {! [::tensor::exists weightsT]} { ::tensor::create weightsT -size $nresidues }
    weightsT = scalar 1.0
    if {! [::tensor::exists rmsdT]} { ::tensor::create rmsdT -size $nresidues }
    if {! [::tensor::exists tmp]} { ::tensor::create tmp -size $nresidues }
    if {! [::tensor::exists tmp_sum]} { ::tensor::create tmp_sum -size 1 }
    if {! [::tensor::exists conv_vector]} { ::tensor::create conv_vector -size $nresidues }
    if {!$ifit_hack} {
      for {set res 0} {$res < $nresidues} {incr res} {
        lappend ifit_rmsd 0.0
      }
    }
  }

  # Iterate over fitting and weighting niter times
  set iter 1
  if {$ifit_convsw} {
    if {$ifit_tensor} {
      conv_vector = scalar 1.0
    } else {
      # TODO: convergence?
    }
  }

  set signalstep [expr {$nsteps/double(20)}]
  puts "Number of steps per iteration: $nsteps"
  puts "Signal every $signalstep steps"

  while {$iter <= $ifit_niter} {
    puts -nonewline "Iteration: $iter "

    # Reset rmsd by residue
    if {$ifit_tensor} {
      rmsdT = scalar 0.0
    } else {
      # TODO: convert to foreach
      for {set res 0} {$res < $nresidues} {incr res} {
        lset rmsd $res 0.0
      }
    }

    set count 0
    for {set i 0} {$i < $nmols} {incr i} {
      set mol1 [lindex $target_mol $i]
      for {set j 0} {$j < $nframes($mol1)} {incr j} {
        $sel_ref($mol1) frame $j
        set sel_reference $sel_ref($mol1)
        if {!$ifit_hack} {
          for {set res 0} {$res < $nresidues} {incr res} {
            $sel_res_ref($mol1:$res) frame $j
          }
        }

        for {set k 0} {$k < $nmols} {incr k} {
          if {$i > $k} {continue}
          set mol2 [lindex $target_mol $k]
          if {$i == $k} {
            set l0 [expr {$j+1}]
          } else {
            set l0 0
          }
          # TODO: put above in $l0  below as an expression
          for {set l $l0} {$l < $nframes($mol2)} {incr l} {

            incr count
            if {$count > $signalstep} {
              puts -nonewline "."
              flush stdout
              set count 0
            }

            $sel_current($mol2) frame $l
            $sel_move($mol2) frame $l

            $sel_move($mol2) move [measure fit $sel_current($mol2) $sel_reference weight user]
            # TODO: this could be tidy up
            if {$ifit_tensor} {
              if {$ifit_hack} {
                lassign [measure rmsd $sel_reference $sel_current($mol2) byres] global_rmsd ifit_rmsd
              } else {
                for {set res 0} {$res < $nresidues} {incr res} {
                  $sel_res($mol2:$res) frame $l
                  lset ifit_rmsd $res [measure rmsd $sel_res_ref($mol1:$res) $sel_res($mol2:$res)]
                }
              }
              rmsdT += array $ifit_rmsd
            } else {
              if {$ifit_hack} {
                lassign [measure rmsd $sel_reference $sel_current($mol2) byres] global_rmsd ifit_rmsd
                for {set res 0} {$res < $nresidues} {incr res} {
                  lset rmsd $res [expr {[lindex $rmsd $res] + [lindex $ifit_rmsd $res]}]
                }
              } else {
                for {set res 0} {$res < $nresidues} {incr res} {
                  $sel_res($mol2:$res) frame $l
                  lset rmsd $res [expr {[lindex $rmsd $res] + [measure rmsd $sel_res_ref($mol1:$res) $sel_res($mol2:$res)]}]
                }
              }
            }
          }
        }
      }
    }

    # Compute mean, mix and max
    if {$ifit_tensor} {
      rmsdT /= scalar $nsteps
      lassign [rmsdT minmax] rmsd_min rmsd_max

      # For convineance later
      set rmsd [rmsdT]
    } else {
      set rmsd_min [expr {[lindex $rmsd 0] / $nsteps}]
      set rmsd_max $rmsd_min
      for {set res 0} {$res < [llength $rmsd]} {incr res} {
        lset rmsd $res [expr {[lindex $rmsd $res] / $nsteps}]
        if {[lindex $rmsd $res] < $rmsd_min} {
          set rmsd_min [lindex $rmsd $res]
          continue
        }
        if {[lindex $rmsd $res] > $rmsd_max} {
          set rmsd_max [lindex $rmsd $res]
        }
      }
    }

    # Compute weights
    if {$ifit_tensor} {
      switch $ifit_type {
        linear {
          weightsT = tensor rmsdT
          weightsT *= scalar $ifit_factor
        }
        exp {
          weightsT = tensor rmsdT
          weightsT *= scalar -$ifit_factor
          weightsT apply exp
        }
        expmin {
          weightsT = tensor rmsdT
          weightsT -= scalar $rmsd_min
          weightsT *= scalar -$ifit_factor
          weightsT apply exp
        }
        minmax {
          if {$rmsd_max == $rmsd_min} {#TODO
            set weight 1
          } else {
            weightsT expr { ($rmsd_max - rmsdT) / ($rmsd_max - $rmsd_min) }
            weightsT = scalar $rmsd_max
            weightsT -= tensor rmsdT
            weightsT /= scalar ($rmsd_max - $rmsd_min)
          }
        }
        gaussian {
          weightsT expr { exp(-(rmsdT * rmsdT) / $ifit_factor) }
          weightsT = tensor rmsdT
          weightsT *= tensor rmsdT
          weightsT /= scalar -$ifit_factor
          weightsT apply exp
        }
      }

      # For convineance later
      set weights [weightsT]
    } else {
      for {set res 0} {$res < $nresidues} {incr res} {
        set r [lindex $rmsd $res]
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
        #puts [format "\tRes: %5d %6.3f %6.3f %6.3f %6.2f" $res [lindex $rmsd $res] $rmsd_min $rmsd_max $weight]
        lset weights $res $weight
      }
    }

    # Update molecules
    foreach mol $target_mol {
      for {set i 0} {$i < [molinfo $mol get numframes]} {incr i} {
        if {$ifit_natoms} {
          for {set res 0} {$res < $nresidues} {incr res} {
            $sel_res($mol:$res) frame $i
            $sel_res($mol:$res) set user [lindex $weights $res]
            $sel_res($mol:$res) set beta [lindex $rmsd $res]
          }
        } else {
          $sel_ref($mol) frame $i
          $sel_ref($mol) set user $weights
          $sel_ref($mol) set beta $rmsd
        }
      }
    }

    # Plot
    if {$ifit_plot} {
      [namespace current]::iterativeFit_plot $iter $plothandle $residues([lindex $target_mol 0]) $rmsd
    }

    # Save
    if {$ifit_save} {
      puts -nonewline $fid [[namespace current]::iterativeFit_save $iter $sel_ref([lindex $target_mol 0])]
    }

    # Convergence
    if {$ifit_convsw} {
      if {$ifit_tensor} {
        if {$iter == 1} {
          conv_vector = tensor weightsT
          tmp = tensor weightsT
          tmp *= tensor weightsT
          tmp_sum = contraction tmp 0
          set t [expr {sqrt([tmp_sum])}]
          set conv 0
        } else {
          tmp = tensor weightsT
          tmp -= tensor conv_vector
          tmp *= tensor tmp
          tmp_sum = contraction tmp 0
          set t [expr {sqrt([tmp_sum])}]
        }
        set conv [expr abs($conv - $t)]
        puts -nonewline [format " %8.5f" $conv]
        if {$conv < $ifit_conv} {
          puts "\nConvergence achieved"
          break
        }
        set conv $t
      } else {
        # TODO: convergence?
      }
    }

    if {$ifit_repre} {[namespace current]::iterativeFit_representation_update}

    if {$ifit_update} {display update}

    incr iter
    puts ""
  }

  puts "Done fitting"

  # Update atom properties
  # TODO: needed?
  #  foreach mol $target_mol  {
  #    for {set i 0} {$i < [molinfo $mol get numframes]} {incr i} {
  #      $sel_ref($mol) frame $i
  #      $sel_ref($mol) set user $weights
  #      $sel_ref($mol) set beta $rmsd
  #    }
  #  }

  if {$ifit_repre} {[namespace current]::iterativeFit_representation_update}

  if {$ifit_update} {display update}

  if {$ifit_plot} {$plothandle replot}

  if {$ifit_save} {close $fid}


  # Clustering
  if {$ifit_cluster} {
    puts "Clustering"
    set fid1 [open "$ifit_file.cluster" w]
    puts $fid1 [format "%6s %6s %7s %7s %7s" "pair1" "pair2" "rms" "rmsw" "lrmsw"]

    if {$ifit_cluster_only} {
      set weights [[atomselect $ifit_cluster_weights_from $sel1] get user]
      if {$ifit_tensor} {
        weightsT = array $weights
      }
    }
    if {$ifit_tensor} {
      tmp_sum = contraction weightsT 0
      set weight_sum [tmp_sum]
    }

    for {set i 0} {$i < $nmols} {incr i} {
      set mol1 [lindex $target_mol $i]
      for {set j 0} {$j < $nframes($mol1)} {incr j} {
        $sel_ref($mol1) frame $j
        set sel_reference $sel_ref($mol1)
        if {!$ifit_hack} {
          for {set res 0} {$res < $nresidues} {incr res} {
            $sel_res_ref($mol1:$res) frame $j
          }
        }

        for {set k 0} {$k < $nmols} {incr k} {
          if {$i > $k} {continue}
          set mol2 [lindex $target_mol $k]
          if {$i == $k} {
            set l0 [expr {$j+1}]
          } else {
            set l0 0
          }
          # TODO: put above in $l0  below as an expression
          for {set l $l0} {$l < $nframes($mol2)} {incr l} {

            $sel_current($mol2) frame $l
            if {$ifit_cluster_fit} {
              $sel_move($mol2) frame $l
              if {$ifit_tensor} {
                $sel_move($mol2) move [measure fit $sel_current($mol2) $sel_reference weight [weightsT]]
              } else {
                $sel_move($mol2) move [measure fit $sel_current($mol2) $sel_reference weight user]
              }
            }

            set cluster_rms [measure rmsd $sel_current($mol2) $sel_reference]

            if {$ifit_hack} {
              lassign [measure rmsd $sel_current($mol2) $sel_reference byres] foo rmsd_by_res
              if {$ifit_tensor} {
                set cluster_rmsw [measure rmsd $sel_current($mol2) $sel_reference weight [weightsT]]

                tmp = array $rmsd_by_res
                tmp_sum = multensor weightsT tmp {0 0}
                set cluster_ifit [expr {[tmp_sum] / $weight_sum}]
              } else {
                set cluster_rmsw [measure rmsd $sel_current($mol2) $sel_reference weight user]
                set cluster_ifit 0.0
                set weight_tot 0.0
                for {set res 0} {$res < $nresidues} {incr res} {
                  $sel_res($mol2:$res) frame $l
                  set rmsd [lindex $rmsd_by_res $res]
                  set weight [lindex [$sel_res_ref($mol1:$res) get user] 0]
                  set weight_tot [expr {$weight_tot + $weight}]
                  set cluster_ifit [expr {$cluster_ifit + $weight * $rmsd}]
                }
                set cluster_ifit [expr {$cluster_ifit / $weight_tot}]
              }

            } else {
              set cluster_rmsw [measure rmsd $sel_current($mol2) $sel_reference weight user]
              set cluster_ifit 0.0
              set weight_tot 0.0
              for {set res 0} {$res < $nresidues} {incr res} {
                $sel_res($mol2:$res) frame $l
                set rmsd [measure rmsd $sel_res_ref($mol1:$res) $sel_res($mol2:$res)]
                set weight [lindex [$sel_res_ref($mol1:$res) get user] 0]
                set weight_tot [expr {$weight_tot + $weight}]
                set cluster_ifit [expr {$cluster_ifit + $weight * $rmsd}]
              }
              set cluster_ifit [expr {$cluster_ifit / $weight_tot}]
            }

            puts $fid1 [format "%2d:%-3d %2d:%-3d %7.3f %7.3f %7.3f" $mol1 $j $mol2 $l $cluster_rms $cluster_rmsw $cluster_ifit]
          }
        }
      }
    }

    close $fid1
  }
}


proc rmsdtt::iterativeFit_test_fast {} {
  variable ifit_hack 0
  variable ifit_tensor 0

  puts "Testing for fast algorithm..."
  if [catch { set ret [measure rmsd [atomselect top "residue 1 to 2"] [atomselect top "residue 1 to 2"] byatom] } msg] {
    puts "   :-( Hacked VMD not found! Please contact the RMSDTT developer"
  } else {
    puts "   :-) Hacked VMD found"
    set ifit_hack 1
  }

  if [catch {package require Tensor} msg] {
    puts "   :-( Package Tensor not found. Please install Tensor for you Tcl/Tk version first"
  } else {
    puts "   :-) Package Tensor found"
    set ifit_tensor 1
  }
}

proc rmsdtt::iterativeFit_check_atoms {} {
  variable datalist

  set target_mol [$datalist(id) get 0 end]
  set nmols [llength $target_mol]
  set sel1 [set_sel]

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
}

proc rmsdtt::iterativeFit_representation {} {
  variable datalist
  variable ifit_scale
  variable ifit_repre
  variable ifit_sel2
  variable ifit_style
  variable ifit_frames_all
  variable ifit_replace
  variable ifit_rep

  set target_mol [$datalist(id) get 0 end]

  mol rep $ifit_style
  mol color User
  mol selection $ifit_sel2
  foreach mol $target_mol {
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
  }

  # Change scale
  if {$ifit_scale != "unchanged"} {
    color scale method $ifit_scale
  }
}

proc rmsdtt::iterativeFit_representation_update {} {
  variable datalist
  variable ifit_rep

  set target_mol [$datalist(id) get 0 end]
  foreach mol $target_mol {
    mol modcolor [mol repindex $mol $ifit_rep($mol)] $mol User
  }

}

proc rmsdtt::iterativeFit_save_header {} {
  return [format "%4s %7s %5s %4s %5s %7s %5s\n" "iter" "residue" "resid" "name" "chain" "mean" "w"]
}

proc rmsdtt::iterativeFit_save {iter sel} {
  set result {}
  set data [$sel get {residue resid resname chain beta user}]
  for {set res 0} {$res < [llength $data]} {incr res} {
    lassign [lindex $data $res] residue resid resname chain rmsd weight
    append result [format "%4s %7d %5d %4s %5s %7.3f %5.3f\n" $iter $residue $resid $resname $chain $rmsd $weight]
  }
  return $result
}

proc rmsdtt::iterativeFit_plot_header {} {
  variable ifit_plot

  if [catch {package require multiplot} msg] {
    showMessage "Plotting in Multiplot not available: package multiplot not installed!\nDo you have the latest VMD version?"
    set ifit_plot 0
    return 0
  } else {
    set title "Rmsd vs Residue"
    set xlab "Residue"
    set ylab "Rmsd (A)"
    set plothandle [multiplot -title $title -xlabel $xlab -ylabel $ylab -nostats]
    return $plothandle
  }
}

proc rmsdtt::iterativeFit_plot {iter plothandle residues rmsd} {
  set color [index2rgb $iter]
  set legend "Iter $iter"
  $plothandle add $residues $rmsd -marker point -radius 2 -fillcolor $color -linecolor $color -nostats -legend $legend
}
