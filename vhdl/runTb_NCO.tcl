# To start ModelSim and source this script from the command line, type this:
#       vsim -do runTb_NCO.tcl
# Or, if ModelSim is already running with the correct working directory, type this in the ModelSim main window:
#       source runTb_NCO.tcl



set library_file_list {
                        design {pck_myhdl_011.vhd
                                NCO.vhd}
                        test   {tb_NCO.vhd}
}

set top_level   work.tb_NCO
set wave_patterns {
                    /*
}
set wave_radices {
                        decimal {o_cos o_sin}
                        unsigned {i_fw i_ph0}
                        
}

proc r  {} {uplevel #0 source runTb_NCO.tcl}
proc rr {} {global last_compile_time
            set last_compile_time 0
            r                            }
proc q  {} {quit -force                  }


#Does this installation support Tk?
set tk_ok 1
if [catch {package require Tk}] {set tk_ok 0}

# Prefer a fixed point font for the transcript
set PrefMain(font) {Courier 10 roman normal}

# Compile out of date files
set time_now [clock seconds]
if [catch {set last_compile_time}] {
  set last_compile_time 0
}
vlib work

foreach {library file_list} $library_file_list {
  vlib $library
  vmap work $library
  if {$library == {test}} {
    set vhdl_ver -2008
  } else {
    set vhdl_ver -93
  }  
  foreach file $file_list {
    if { $last_compile_time < [file mtime $file] } {
      if [regexp {.vhdl?$} $file] {
        vcom $vhdl_ver $file
      } else {
        vlog $file
      }
      set last_compile_time 0
    }
  }
}
set last_compile_time $time_now

# Load the simulation
eval vsim $top_level -onfinish stop

# If waves are required
if [llength $wave_patterns] {
  noview wave
  foreach pattern $wave_patterns {
    add wave $pattern
  }
  configure wave -signalnamewidth 1
  foreach {radix signals} $wave_radices {
    foreach signal $signals {
      catch {property wave -radix $radix $signal}
    }
  }
}

# Run the simulation
run -all

# If waves are required
if [llength $wave_patterns] {
  if $tk_ok {wave zoomfull}
}

puts {
  Script commands are:

  r = Recompile changed and dependent files
 rr = Recompile everything
  q = Quit without confirmation
}

# How long since project began?
if {[file isfile start_time.txt] == 0} {
  set f [open start_time.txt w]
  puts $f "Start time was [clock seconds]"
  close $f
} else {
  set f [open start_time.txt r]
  set line [gets $f]
  close $f
  regexp {\d+} $line start_time
  set total_time [expr ([clock seconds]-$start_time)/60]
  puts "Project time is $total_time minutes"
}