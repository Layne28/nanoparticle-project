#!/usr/bin/tclsh

proc initial_config {} {
	mol new init_config.xyz waitfor all autobonds off
	mol delrep 0 top
	# add new representations
	mol color Element
	mol representation points 15	
	mol selection {name Au}
	mol addrep top
}

proc load_conf {fname} {
	mol new $fname waitfor all autobonds off
	mol delrep 0 top
	# add new representations
	mol color Element
	mol representation points 15	
	mol selection {name Au}
	mol addrep top
}

proc load_ana {fname} {
	mol new $fname waitfor all autobonds off
	mol delrep 0 top
	# add new representations
	mol color Name
	mol representation points 15	
	mol selection {name Ag}
	mol addrep top
}
proc start {} {
	mol new traj.xyz waitfor all autobonds off
	mol delrep 0 top
	# add new representations
	mol color Element
	mol representation points 15	
	mol selection {name Au}
	mol addrep top
}

proc makemovie {numframes} {
  set i 0
  set fcount 0
  while {$fcount < $numframes} {
    mol new [format "conf_%04d.xyz" $i] waitfor -1 autobonds off
    mol delrep 0 top
    # add new representations
    mol color Element
    mol representation QuickSurf 
#    mol material EdgyShiny
    mol selection {name Au}
    mol addrep top
    axes location off
    color Display Background white
#    display depthcue off
    display resize 800 800
#    display resetview
    scale by 0.8
    render TachyonInternal [format "movie_%04d.tga" $i]
    set i [expr $i + 1]
    incr fcount 
    mol delete all
  }
}
    
