#-----------------------------------------------------------NEERAJ BALACHANDAR-----------------------------------------------------------------------------------------------------
__authors__ = "NEERAJ BALACHANDAR"
__contact__ = "neerajbalachandar@gmail.com"

#-------------------------------------------MID FIDELITY WITH PARTICLE FIELD AND FLUID DOMAIN CALCULATION. UVLM,rVPM simulation-------------------------------------------------------------------
#-------------------------------------------SIMULATION DEFINITION WITH POST PROCESSING--------------------------------------------------------------------------------------

import FLOWUnsteady as uns
import FLOWUnsteady: vlm, vpm, gt, Im, dot, norm

include("maneuver_definition.jl")     
include("vehicle_definition.jl")    
include("monitor_definition.jl")   

# =============================================================================
# SIMULATION PARAMETERS
# =============================================================================

# ----------------- GENERAL SIMULATION SETUP ----------------------------------
run_name        = "dragonfly"    # Name of this simulation
save_path       = ""  # Where to save results
paraview        = true                              # Whether to visualize with Paraview
compute_fluid_domain = true                         # Whether to compute fluid domain post-simulation

# ----------------- GEOMETRY PARAMETERS ----------------------------------------
n_factor        = 2                         # Discretization factor
add_wings       = true                      # Whether to include wings
add_rotors      = false                     # No rotors for dragonfly

# Reference lengths (typical dragonfly dimensions)
wingspan        = 1.0                      # (m) dragonfly wingspan ~10cm

# ----------------- FLIGHT PARAMETERS ------------------------------------------
# Flapping kinematics
flap_amplitude  = 60.0                      # (degrees) flapping amplitude
flap_frequency  = 30.0                      # (Hz) flapping frequency
vehicle_velocity = 1.0                      # (m/s) forward flight velocity
angle_of_attack = 5.0                       # (degrees) vehicle angle of attack

# Flight conditions
Vinf(X,t)       = 20.0*[1.0, 0.0, 0.0]  # (m/s) freestream velocity
rho             = 1.225                     # (kg/m^3) air density
mu              = 1.81e-5                   # (kg/ms) air dynamic viscosity

# Time parameters
flap_cycles     = 3                         # Number of flapping cycles to simulate
ttot            = flap_cycles / flap_frequency  # (s) total simulation time
nsteps          = 150 * flap_cycles         # Time steps (150 per cycle for efficiency)
dt              = ttot / nsteps             # (s) time step

# Simulation time window
tstart          = 0.00 * ttot               # (s) start time
tquit           = 1.00 * ttot               # (s) end time

# ----------------- SOLVER PARAMETERS ------------------------------------------

# Aerodynamic solver
VehicleType     = uns.UVLMVehicle           # Unsteady solver for flapping wings

# VPM particle shedding
p_per_step      = 2                         # Particles shed per time step (reduced for efficiency)
shed_starting   = false                     # No starting vortex needed
shed_unsteady   = true                      # Shed vorticity from unsteady loading
unsteady_shedcrit = 0.01                    # Circulation fluctuation threshold

# Regularization parameters
sigma_vlm_surf  = wingspan/150              # VLM-on-VPM smoothing radius
lambda_vpm      = 2.125                     # VPM core overlap
sigma_vpm_overwrite = lambda_vpm * flap_frequency * flap_amplitude * pi/180 * dt / p_per_step
sigmafactor_vpmonvlm = 1                    # VPM-on-VLM velocity factor

# Wing solver parameters
vlm_rlx         = 0.4                       # VLM relaxation factor
vlm_vortexsheet = false                     # Actuator surface model
vlm_vortexsheet_overlap = 2.125
vlm_vortexsheet_distribution = uns.g_pressure
vlm_vortexsheet_sigma_tbv = sigma_vpm_overwrite
vlm_vortexsheet_maxstaticparticle = vlm_vortexsheet ? 1000000 : nothing

# Force calculation parameters
KJforce_type    = "regular"                 # Kutta-Joukowski force type
include_trailingboundvortex = false
include_unsteadyforce = true                # Critical for flapping wing aerodynamics
add_unsteadyforce = false
include_parasiticdrag = true
add_skinfriction = true
calc_cd_from_cl = false
wing_polar_file = "xf-n0012-il-500000-n5.csv"  # Thin airfoil for dragonfly wings

# VPM solver settings
vpm_integration = vpm.euler                 # Temporal integration scheme
vpm_viscous     = vpm.Inviscid()           # Viscous diffusion

# Subfilter-scale model
vpm_SFS = vpm.DynamicSFS(vpm.Estr_fmm, vpm.pseudo3level_positive;
                         alpha=0.999, maxC=1.0,
                         clippings=[vpm.clipping_backscatter],
                         controls=[vpm.control_directional, vpm.control_magnitude])

# ----------------- FLUID DOMAIN PARAMETERS --------------------------------
# Grid parameters for fluid domain computation
L_domain        = wingspan                  # (m) reference length for domain
dx_domain       = L_domain/40               # (m) cell size in each direction
dy_domain       = L_domain/40
dz_domain       = L_domain/40

# Domain bounds (adjusted for dragonfly)
Pmin_domain     = L_domain*[-1.0, -1.5, -1.5]   # (m) minimum bounds
Pmax_domain     = L_domain*[ 3.0,  1.5,  1.5]   # (m) maximum bounds
NDIVS_domain    = ceil.(Int, (Pmax_domain .- Pmin_domain)./[dx_domain, dy_domain, dz_domain])
nnodes_domain   = prod(NDIVS_domain .+ 1)        # Total number of nodes

# VPM settings for fluid domain
maxparticles_domain = Int(1.0e6 + nnodes_domain)
fmm_domain      = vpm.FMM(; p=4, ncrit=50, theta=0.4, phi=0.3)
f_sigma_domain  = 0.5                       # Smoothing of node particles
maxsigma_domain = L_domain/15               # Maximum particle size
maxmagGamma_domain = Inf                    # Maximum vortex strength

# File naming for fluid domain
pfield_prefix   = "dragonfly_pfield"        # Prefix of particle field files
staticpfield_prefix = "dragonfly_staticpfield"  # Prefix of static particle field files
fdom_prefix     = "dragonfly"               # Prefix of fluid domain files

# =============================================================================
# SIMULATION EXECUTION
# =============================================================================

println("="^80)
println("DRAGONFLY FLAPPING WING SIMULATION WITH FLUID DOMAIN ANALYSIS")
println("="^80)

# ----------------- 1) VEHICLE DEFINITION --------------------------------------
println("\n1) Generating dragonfly vehicle geometry...")

vehicle = dragonfly_geometry()

println("   Vehicle generated with $(length(vehicle.system.wings)) wings")
println("   Total panels: $(sum(uns.vlm.get_m(wing) for wing in vehicle.system.wings))")

# ----------------- 2) MANEUVER DEFINITION -------------------------------------
println("\n2) Generating flapping maneuver...")

maneuver = generate_dragonfly_maneuver(;
    disp_plot = false,  # Disable plot for automated run
    add_wings = add_wings,
    vehicle_velocity = vehicle_velocity,
    angle_of_attack = angle_of_attack
)

println("   Maneuver configured: $(flap_cycles) cycles at $(flap_frequency) Hz")

# ----------------- 3) SIMULATION DEFINITION -----------------------------------
println("\n3) Setting up simulation...")

# Reference parameters
Vref = vehicle_velocity
RPMref = flap_frequency * 60                # Reference RPM equivalent

# Initial conditions
t0 = tstart / ttot
Vinit = Vref * maneuver.Vvehicle(t0)
Winit = pi/180 * (maneuver.anglevehicle(t0 + 1e-12) - maneuver.anglevehicle(t0)) / (ttot * 1e-12)

# Create simulation object
simulation = uns.Simulation(vehicle, maneuver, Vref, RPMref, ttot;
                           Vinit=Vinit, Winit=Winit, t=tstart)

# Memory allocation for particles
max_particles = ceil(Int, (nsteps + 2) * (2 * vlm.get_m(vehicle.wake_system) * (p_per_step + 1) + p_per_step))
max_particles = min(3000000, max_particles)  # Reasonable limit for dragonfly

println("   Simulation configured: $nsteps steps, max $max_particles particles")

# ----------------- 4) FORCE CALCULATION SETUP --------------------------------
println("\n4) Setting up force calculations...")

forces = []

# Kutta-Joukowski force (primary lift mechanism)
kuttajoukowski = uns.generate_aerodynamicforce_kuttajoukowski(
    KJforce_type,
    sigma_vlm_surf, 0.0,  # No rotor surface
    vlm_vortexsheet, vlm_vortexsheet_overlap,
    vlm_vortexsheet_distribution,
    vlm_vortexsheet_sigma_tbv;
    vehicle=vehicle
)
push!(forces, kuttajoukowski)

# Unsteady force (critical for flapping aerodynamics)
if include_unsteadyforce
    unsteady(args...; optargs...) = uns.calc_aerodynamicforce_unsteady(
        args...; add_to_Ftot=add_unsteadyforce, optargs...
    )
    push!(forces, unsteady)
end

# Parasitic drag force
if include_parasiticdrag
    parasiticdrag = uns.generate_aerodynamicforce_parasiticdrag(
        wing_polar_file;
        calc_cd_from_cl=calc_cd_from_cl,
        add_skinfriction=add_skinfriction
    )
    push!(forces, parasiticdrag)
end

# Combined force calculation function
function calc_aerodynamicforce_fun(vlm_system, args...; per_unit_span=false, optargs...)
    fieldname = per_unit_span ? "ftot" : "Ftot"
    if fieldname in keys(vlm_system.sol)
        pop!(vlm_system.sol, fieldname)
    end
    
    Ftot = nothing
    for force in forces
        Ftot = force(vlm_system, args...; per_unit_span=per_unit_span, optargs...)
    end
    
    return Ftot
end

# ----------------- 5) MONITORS SETUP ------------------------------------------
println("\n5) Setting up monitoring functions...")

# Wing monitor options
wingmonitor_optargs = (
    include_trailingboundvortex=include_trailingboundvortex,
    calc_aerodynamicforce_fun=calc_aerodynamicforce_fun
)

# Generate monitoring functions
monitors = generate_monitor_dragonfly(
    vehicle, rho, Vinf, nsteps, save_path;
    add_wings=add_wings,
    wingmonitor_optargs=wingmonitor_optargs
)

# ----------------- 6) WAKE TREATMENT ------------------------------------------
println("\n6) Configuring wake treatment...")

# Remove weak particles
rmv_strength = 2 * 2 / p_per_step * dt / (1 / flap_frequency)
minmaxGamma = rmv_strength * [0.0001, 0.15]
wake_treatment_strength = uns.remove_particles_strength(
    minmaxGamma[1]^2, minmaxGamma[2]^2; every_nsteps=10
)

# Remove oversized particles
minmaxsigma = sigma_vpm_overwrite * [0.1, 6]
wake_treatment_sigma = uns.remove_particles_sigma(
    minmaxsigma[1], minmaxsigma[2]; every_nsteps=10
)

# Remove distant particles
wake_treatment_sphere = uns.remove_particles_sphere(
    (2.5 * wingspan)^2, 1; Xoff=[0.5 * wingspan, 0, 0]
)

# Combine wake treatments
wake_treatment = uns.concatenate(
    wake_treatment_sphere, 
    wake_treatment_strength, 
    wake_treatment_sigma
)

# Combined runtime function
runtime_function = uns.concatenate(wake_treatment, monitors)

# ----------------- 7) RUN SIMULATION ------------------------------------------
println("\n7) Running simulation...")
println("   Simulating $flap_cycles flapping cycles over $(round(ttot, digits=3)) seconds")
println("   Forward velocity: $(vehicle_velocity) m/s")
println("   Angle of attack: $(angle_of_attack)°")

# Execute simulation
uns.run_simulation(simulation, nsteps;
    # Simulation options
    Vinf=Vinf,
    rho=rho, mu=mu, 
    tquit=tquit,
    
    # Solver options
    p_per_step=p_per_step,
    max_particles=max_particles,
    max_static_particles=vlm_vortexsheet_maxstaticparticle,
    vpm_integration=vpm_integration,
    vpm_viscous=vpm_viscous,
    vpm_SFS=vpm_SFS,
    sigma_vlm_surf=sigma_vlm_surf,
    sigma_rotor_surf=0.0,  # No rotors
    sigma_vpm_overwrite=sigma_vpm_overwrite,
    sigmafactor_vpmonvlm=sigmafactor_vpmonvlm,
    vlm_vortexsheet=vlm_vortexsheet,
    vlm_vortexsheet_overlap=vlm_vortexsheet_overlap,
    vlm_vortexsheet_distribution=vlm_vortexsheet_distribution,
    vlm_vortexsheet_sigma_tbv=vlm_vortexsheet_sigma_tbv,
    vlm_rlx=vlm_rlx,
    hubtiploss_correction=nothing,  # No rotors
    shed_starting=shed_starting,
    shed_unsteady=shed_unsteady,
    unsteady_shedcrit=unsteady_shedcrit,
    extra_runtime_function=runtime_function,
    
    # Output options
    save_path=save_path,
    run_name=run_name,
    save_wopwopin=false,
)

println("   Simulation completed successfully!")

# =============================================================================
# FLUID DOMAIN COMPUTATION
# =============================================================================

if compute_fluid_domain
    println("\n" * "="^80)
    println("FLUID DOMAIN COMPUTATION")
    println("="^80)
    
    # ----------------- FLUID DOMAIN SETUP -------------------------------------
    println("\n8) Setting up fluid domain computation...")
    
    # Time steps to process (select representative steps)
    # Process every 10th step for the last cycle
    last_cycle_start = nsteps - 150  # Last cycle
    fluid_domain_nums = collect(last_cycle_start:10:nsteps)
    
    println("   Processing $(length(fluid_domain_nums)) time steps for fluid domain")
    println("   Domain grid: $(NDIVS_domain) cells, $(nnodes_domain) nodes")
    
    # Fluid domain save path
    fdom_save_path = joinpath(save_path, "fluid_domain")
    
    # Create save path
    gt.create_path(fdom_save_path, false)  # No prompt for automated execution
    
    # Grid orientation (aligned with vehicle)
    Oaxis_domain = gt.rotation_matrix2(0, 0, angle_of_attack)
    
    # Include static particles for complete analysis
    include_staticparticles = true
    other_file_prefs = include_staticparticles ? [staticpfield_prefix] : []
    other_read_paths = [save_path for i in 1:length(other_file_prefs)]
    
    # Generate preprocessing function for particle field
    preprocessing_pfield = uns.generate_preprocessing_fluiddomain_pfield(
        maxsigma_domain, maxmagGamma_domain;
        verbose=true, v_lvl=1
    )
    
    # ----------------- COMPUTE FLUID DOMAIN -----------------------------------
    println("\n9) Computing fluid domain...")
    
    # Process fluid domain computation
    uns.computefluiddomain(
        Pmin_domain, Pmax_domain, NDIVS_domain,
        maxparticles_domain,
        fluid_domain_nums, save_path, pfield_prefix;
        Oaxis=Oaxis_domain,
        fmm=fmm_domain,
        f_sigma=f_sigma_domain,
        save_path=fdom_save_path,
        file_pref=fdom_prefix, 
        grid_names=["_fdom"],
        other_file_prefs=other_file_prefs,
        other_read_paths=other_read_paths,
        userfunction_pfield=preprocessing_pfield,
        verbose=true, 
        v_lvl=0
    )
    
    println("   Fluid domain computation completed!")
    println("   Results saved to: $fdom_save_path")
    
else
    println("\n   Skipping fluid domain computation (compute_fluid_domain = false)")
end

# =============================================================================
# POST-PROCESSING AND VISUALIZATION
# =============================================================================

println("\n" * "="^80)
println("POST-PROCESSING")
println("="^80)

# ----------------- RESULTS SUMMARY ----------------------------------------
println("\n10) Simulation summary:")
println("    • Main results: $save_path")
if compute_fluid_domain
    println("    • Fluid domain: $(joinpath(save_path, "fluid_domain"))")
end
println("    • VTK files generated for:")
println("      - Vehicle geometry and wake")
println("      - Particle field (pfield)")
println("      - Static particle field (staticpfield)")
if compute_fluid_domain
    println("      - Fluid domain (fdom)")
end

# ----------------- PARAVIEW VISUALIZATION ---------------------------------
if paraview
    println("\n11) Opening results in ParaView...")
    
    # Find the main VTK file
    vtk_files = filter(x -> endswith(x, ".vtk"), readdir(save_path))
    if !isempty(vtk_files)
        main_vtk = joinpath(save_path, vtk_files[1])
        println("    Opening: $main_vtk")
        try
            run(`paraview --data=$main_vtk`, wait=false)
        catch e
            println("    Could not launch ParaView automatically: $e")
            println("    Please open $main_vtk manually in ParaView")
        end
    end
    
    # Also mention fluid domain files
    if compute_fluid_domain
        fdom_path = joinpath(save_path, "fluid_domain")
        fdom_files = filter(x -> endswith(x, ".vtk"), readdir(fdom_path))
        if !isempty(fdom_files)
            println("    Fluid domain files available in: $fdom_path")
        end
    end
end

# ----------------- ANALYSIS NOTES -----------------------------------------
println("\n12) Analysis notes:")
println("    • Force time histories are saved in monitor output files")
println("    • Wing performance can be analyzed separately for fore and hind wings")
println("    • Vortex formation data is tracked during simulation")
if compute_fluid_domain
    println("    • Fluid domain provides detailed velocity and vorticity fields")
    println("    • NOTE: Freestream velocity must be added manually in ParaView")
end

println("\n" * "="^80)
println("DRAGONFLY SIMULATION COMPLETE")
println("="^80)

