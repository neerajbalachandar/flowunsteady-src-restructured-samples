#-----------------------------------------------------------NEERAJ BALACHANDAR-----------------------------------------------------------------------------------------------------
__authors__ = "NEERAJ BALACHANDAR"
__contact__ = "neerajbalachandar@gmail.com"

"""

"""

#-----------------------------------------------------------FLOWUnsteady_openvsp.jl--------------------------------------------------------------------------------
function import_dragonfly_geometry()
"""
    import_dragonfly_geometry() -> UVLMVehicle
  
Import and process the dragonfly geometry from VSP degenerate geometry file.
Returns a configured UVLM vehicle with all wing systems.
"""
    
    # Read degenerate geometry components
    components = uns.read_degengeom(joinpath(@__DIR__, "dragonfly_DegenGeom.csv"))
    
    # Extract individual components by name
    body_component = components[findfirst(c -> c.name == "Body", components)]
    fore_wing_component = components[findfirst(c -> c.name == "Fore Wing", components)]
    hind_wing_component = components[findfirst(c -> c.name == "Hind Wing", components)]
    
    # Generate body geometry
    body = uns.import_vsp(body_component; geomType="fuselage")
    # Process wing geometries with proper orientation
    fore_wing_left = uns.import_vsp(fore_wing_component; geomType="wing", symmetric=true)
    fore_wing_right = uns.import_vsp(fore_wing_component; geomType="wing")
    hind_wing_left = uns.import_vsp(hind_wing_component; geomType="wing", symmetric=true)
    hind_wing_right = uns.import_vsp(hind_wing_component; geomType="wing")
  
    println("dragonfly geometry imported")

#-------------------------------------------------------------FLOWUnsteady_vehicle_vlm_unsteady.jl--------------------------------------------------------------------------
    # Create body grid
    fuselage_grid = uns.gt.MultiGrid(3)  # Fixed variable name consistency
    uns.gt.addgrid(fuselage_grid, "Body", body)

    # Build comprehensive wing system
    system = uns.vlm.WingSystem()
    uns.vlm.addwing(system, "ForeWing_L", fore_wing_left)
    uns.vlm.addwing(system, "ForeWing_R", fore_wing_right)
    uns.vlm.addwing(system, "HindWing_L", hind_wing_left)
    uns.vlm.addwing(system, "HindWing_R", hind_wing_right)
    
    # Configure tilting systems properly - each tilting system should be a WingSystem
    # containing wings that tilt together
    fore_tilting_system = uns.vlm.WingSystem()
    uns.vlm.addwing(fore_tilting_system, "ForeWing_L", fore_wing_left)
    uns.vlm.addwing(fore_tilting_system, "ForeWing_R", fore_wing_right)
    
    hind_tilting_system = uns.vlm.WingSystem()
    uns.vlm.addwing(hind_tilting_system, "HindWing_L", hind_wing_left)
    uns.vlm.addwing(hind_tilting_system, "HindWing_R", hind_wing_right)
    
    # Define tilting systems as tuple of WingSystems
    tilting_systems = (fore_tilting_system, hind_tilting_system)
    
    # Surface regularization should be handled differently
    # Remove the direct sol access approach as it's not standard
    
    grids = [fuselage_grid]  # Fixed variable name
    vlm_system = system
    wake_system = system
    
    # Set freestream velocity
    uns.vlm.setVinf(system, freestream_velocity)
    
    # Create UVLM vehicle with proper tilting systems
    vehicle = uns.VLMVehicle(system;
        tilting_systems=tilting_systems,  # Add tilting systems
        grids=grids,
        vlm_system=vlm_system,
        wake_system=wake_system
    );
    
    # Check panel counts with proper wing access
    for (i, wing) in enumerate(system.wings)
        m_panels = uns.vlm.get_m(wing)
        wing_name = "Wing_$i"  # Use index-based naming or store names separately
        println("$wing_name: $m_panels panels")
        if m_panels < 50
            @warn "Low panel count detected for $wing_name"
        end
    end
    println("Configuring UVLM vehicle successful")
    
    return vehicle

end
