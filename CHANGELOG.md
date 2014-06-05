### 4.0.0 (2014-06-04)

* Publicly release iterative fitting

### 3.0.1 (2014-06-04)

* Move code to github

### 3.0 (2010-08-05)

* Included in VMD 1.8.8
* New options to highlight equivalent atoms in the selection
* New options to do weighted rmsd and alignment

### 2.01 (2006-04-24)

* Using Average as the reference molecule was not working.
* Syntax of VMD command measure avpos updated (thanks to Christopher Lanci for reporting this bug)

### 2.0 (2006-04-03)

* Included in VMD 1.8.4
* Complete rewrite to speed up the calculations
* Add Multiplot support and redesign the gui

### 1.10 (2006-01-23)

Bug in alignment code solved

### 1.9.2.2 (2005-08-05)

Bug in alignment code solved

### 1.9.2.1 (2005-07-27)

Added option for exclude hydrogens from the atoms selection

### 1.9 (2005-04-06)

* Added option to calculate pairwise rmsd for all frames.
  * If 'all' is selected, each frame of the selected (top,...) molecule will be used as reference
  * Changed procs compute_rms, saveData, doRmsd, ctrltraj and rmsdtt
  * Added checkbutton 'all'
  * Added proc saveDataAll
* Variable names changed: frames_ref to frame_ref, base to rmsd_base
* Other minor changes and bug-solving

### 1.8 (2005-01-18)

* Change to use the VMD "measure rmsd" command instead of calculating the rmsd from scratch, except when using the average as reference
* Average option now works with trajectories for the rmsd button (not align)
* Change order of some subroutines in the code
* Added option to show the results in time instead of frame
* Added option to plot the results in Excel when using Windows
  * Needs the tcom tcl library.
  * VMD will complain upon exit if Excel is still open, but nothing else happens

### 1.7 (2004-10-24)

* 'Backbone' option change to select atoms 'C CA N', which differs from VMD definition, and deactivated by default
* Added 'Trace' option. Will add 'and name CA' to the actual selection
* Multi-line selection text, with added support for comments. Selection widget changed from entry to text
* New warning windows instead of the console warning messages. Added procedure 'showMessage'
* More and improved warnings
* Code split in procedures: onlyactive, doAlign, doRmsd
* Minor bugs

### 1.4 (2004-08-20)

* Initial version
