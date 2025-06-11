#--------------------------------------------------------------------------------NEERAJ BALACHANDAR-----------------------------------------------------------------------------------------------------
__authors__ = "NEERAJ BALACHANDAR"
__contact__ = "neerajbalachandar@gmail.com"

function generate_monitor_dragonfly(vehicle, rho, Vinf, nsteps, save_path;
                                    add_wings=true,
                                    wingmonitor_optargs=[])
    
    # Create output directory
    isdir(save_path) || mkpath(save_path)

    monitors = []

    # -------------------- WING MONITORS ---------------------------------------
    if add_wings
    
        system = vehicle.system
        #Dimensionalised
        b_ref = 1.0
        ar_ref = 1.0
        qinf = 1.0

        # Force component directions (Z-up, X-forward)
        L_dir = [0, 0, 1]  
        D_dir = [-1, 0, 0] 

        # Aerodynamic force calculation
        calc_aerodynamicforce_fun = uns.generate_calc_aerodynamicforce(;
            add_parasiticdrag=true,
            add_skinfriction=true,
            airfoilpolar="xf-n0012-il-500000-n5.csv" # Replace with NACA0003
        )

        fore_system = vehicle.tilting_systems[1]
        monitor_fore = uns.generate_monitor_wing(
            fore_system, Vinf, b_ref, ar_ref, rho, qinf, nsteps;
            calc_aerodynamicforce_fun=calc_aerodynamicforce_fun,
            save_path=save_path,
            run_name="fore_wings",
            CL_lbl="Fore Wing Lift (N)",
            CD_lbl="Fore Wing Drag (N)",
            L_dir=L_dir,
            D_dir=D_dir,
            wingmonitor_optargs...
        )

        hind_system = vehicle.tilting_systems[2]
        monitor_hind = uns.generate_monitor_wing(
            hind_system, Vinf, b_ref, ar_ref, rho, qinf, nsteps;
            calc_aerodynamicforce_fun=calc_aerodynamicforce_fun,
            save_path=save_path,
            run_name="hind_wings",
            CL_lbl="Hind Wing Lift (N)",
            CD_lbl="Hind Wing Drag (N)",
            L_dir=L_dir,
            D_dir=D_dir,
            wingmonitor_optargs...
        )

        push!(monitors, monitor_fore, monitor_hind)
    end

    # -------------------- FLAPPING KINEMATICS --------------------------------
    function monitor_flapping(pfield, vehicle, t, dt)
        # Track flapping angles
        fore_angle = vehicle.tilting_systems[1].angle(t)
        hind_angle = vehicle.tilting_systems[2].angle(t)
        
        # Store in global data structure
        global flapping_data = get!(flapping_data, t, Dict())
        flapping_data[t] = (fore=fore_angle, hind=hind_angle)
        
        return false  # Don't stop simulation
    end

    # -------------------- VORTEX METRICS ----------------------------------
    function monitor_vortex_formation(pfield, vehicle, t, dt)
        # Calculate LEV characteristics
        global vortex_data = get!(vortex_data, t, Dict())
        
        # Implementation for vortex tracking would go here
        # (Requires access to flow field particles)
        
        return false
    end

    # -------------------- STANDARD MONITORS -----------------------------------
    # Vehicle state variables
    state_monitor = uns.generate_monitor_statevariables(;
        save_path=save_path,
        figname="vehicle_states"
    )

    # Flow field enstrophy
    enstrophy_monitor = uns.generate_monitor_enstrophy(;
        save_path=save_path,
        figname="flow_energy"
    )

    # Combine all monitors
    all_monitors = uns.concatenate(
        state_monitor,
        enstrophy_monitor,
        monitor_flapping,
        monitor_vortex_formation,
        monitors...
    )

    return all_monitors
end

