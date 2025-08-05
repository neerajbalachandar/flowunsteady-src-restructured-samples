#-----------------------------------------------------------NEERAJ BALACHANDAR-----------------------------------------------------------------------------------------------------
__authors__ = "NEERAJ BALACHANDAR"
__contact__ = "neerajbalachandar@gmail.com"

"""
    UVLMVehicle{N, M, R}(system; tilting_systems, vlm_system, wake_system, grids)

Type handling all geometries and subsystems that define a vehicle made
out of FLOWVLM components (Wing, WingSystem).

# ARGUMENTS WE HAVE CONSIDERED INCLUDE:
* `system::FLOWVLM.WingSystem`:        System of all FLOWVLM objects. This system
                                    is considered as the entire vehicle. Not all
                                    components in this system will be solved,
                                    but they will all be rotated and translated
                                    according to the maneuver.

* `tilting_systems::NTuple{N, FLOWVLM.WingSystem}`:   Tuple of all FLOWVLM
                                    tilting objects, where `tilting_systems[i]`
                                    contains the i-th FLOWVLM system of lifting
                                    surfaces and rotors that tilt together.

* `vlm_system::FLOWVLM.WingSystem`:    System of all FLOWVLM objects to be solved
                                    through the VLM solver.
* `wake_system::FLOWVLM.WingSystem`:   System of all FLOWVLM objects that will
                                    shed a VPM wake.
* `grids::Array{gt.GridTypes}`:         Array of grids that will be translated and
                                    rotated along with `system`.

# State variables
* `V::Vector`                   : Current vehicle velocity
* `W::Vector`                   : Current vehicle angular velocity
* `prev_data::Array{Any}`       : Information about previous step
* `grid_O::Vector{Vector}`       : Origin of every grid
"""

#-----------------------------------------------------------FLOWUnsteady_openvsp.jl-----------------------------
function wing_geometry()
    
    # Read degenerate geometry components
    components = uns.read_degengeom(joinpath(@__DIR__, "Wing_DegenGeom.csv"))

    wing_component = components[findfirst(c -> c.name == "Wing", components)]
    
    wing_left = uns.import_vsp(wing_component; geomType="wing")
    wing_right = uns.import_vsp(wing_component; geomType="wing", flip_y=true)
  
    println("wing geometry imported")

#-------------------------------------------------------------FLOWUnsteady_vehicle_vlm_unsteady.jl--------------------
    
    # Build comprehensive wing system with proper hierarchy
    system = uns.vlm.WingSystem()
    
    # Add wings to main vehicle system
    uns.vlm.addwing(system, "Wing_L", wing_left)
    uns.vlm.addwing(system, "Wing_R", wing_right)
    
    # Configure tilting systems as WingSystem subgroups
    left_tilting_system = uns.vlm.WingSystem()
    right_tilting_system = uns.vlm.WingSystem()
    
    # Reference wings in tilting systems
    uns.vlm.addwing(left_tilting_system, "Left", wing_left)
    uns.vlm.addwing(right_tilting_system, "Right", wing_right)
    
    # Initialize vehicle with proper kinematic hierarchy
    vehicle = uns.UVLMVehicle(
        system;
        tilting_systems = (left_tilting_system, right_tilting_system),
        vlm_system = system,  # Solve entire system with VLM
        wake_system = system, # Shed wake from all components
        V = zeros(3),  # Initial linear velocity
        W = zeros(3),  # Initial angular velocity
        prev_data = [
            deepcopy(system),  # Previous vlm_system state
            deepcopy(system),  # Previous wake_system state
            ()  # Empty rotor systems
        ]
    )
    
    # Validate panel counts with proper indexing
    for (i, wing) in enumerate(vehicle.system.wings)
        m_panels = uns.vlm.get_m(wing)
        println("Wing $i: $m_panels panels")
        if m_panels < 50
            @warn "Low panel count ($m_panels) detected in wing $i"
        end
    end
    
    println("UVLM vehicle configuration successful")
    return vehicle

    #EXPORTING---------------
    # Define output directory
    save_geom_path = "/home/dysco/FLOWUnsteady/FLOWUnsteady/flapping_wing_journal/pitching plunging wing/geometry/"

    # Clean/create directory
    if isdir(save_geom_path)
        rm(save_geom_path; recursive=true, force=true)
    end
    mkdir(save_geom_path)

    # Set freestream velocity for visualization context
    uns.vlm.setVinf(system, Vinf)

    # Save VTK files and get master PVD file name
    pvd_file = uns.save_vtk(vehicle, "p_p_flapping"; 
                            path=save_geom_path,
                            save_wopwopin=false)  # Set true if needing acoustic inputs

    # Open in ParaView using the master PVD file
    paraview_file = joinpath(save_geom_path, pvd_file)
    run(`paraview --data=$paraview_file`, wait=false)

end
#----------------------------------------------------------------------------------------------------------------------