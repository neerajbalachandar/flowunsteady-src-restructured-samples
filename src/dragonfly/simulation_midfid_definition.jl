
#-----------------------------------------------------------NEERAJ BALACHANDAR-----------------------------------------------------------------------------------------------------
__authors__ = "NEERAJ BALACHANDAR"
__contact__ = "neerajbalachandar@gmail.com"

#-------------------------------------------MID FIDELITY WITH PARTICLE FIELD AND FLUID DOMAIN CALCULATION-------------------------------------------------------------------
#-------------------------------------------SIMULATION DEFINITION WITH POST PROCESSING--------------------------------------------------------------------------------------

#=##############################################################################
# DESCRIPTION
    Complete dragonfly flapping wing simulation with integrated fluid domain
    computation. This implementation runs the UVLM/rVPM simulation and then
    processes the particle field to generate volumetric flow field data for
    detailed aerodynamic analysis.

# AUTHORSHIP
  * Author          : NEERAJ BALACHANDAR
  * Email           : neerajbalachandar@gmail.com
  * Created         : Based on FLOWUnsteady framework
  * License         : MIT
=###############################################################################

import FLOWUnsteady as uns
import FLOWUnsteady: vlm, vpm, gt, Im, dot, norm

# Include dragonfly-specific modules
include("maneuver_definition.jl")      # Kinematic maneuver generation
include("vehicle_definition.jl")      # Vehicle geometry construction
include("monitor_definition.jl")       # Monitoring and analysis functions




# =============================================================================
# SIMULATION PARAMETERS
# =============================================================================

# ----------------- GENERAL SIMULATION SETUP ----------------------------------
run_name        = "dragonfly"    # Name of this simulation
save_path       = "/home/dysco/Neeraj/Dragonfly flowunsteady/modular/results/"  # Where to save results
paraview        = true                              # Whether to visualize with Paraview
compute_fluid_domain = true                         # Whether to compute fluid domain post-simulation




# ----------------- GEOMETRY PARAMETERS ----------------------------------------
n_factor        = 4                         # Discretization factor
add_wings       = true                      # Whether to include wings
add_rotors      = false                     # No rotors for dragonfly

# Reference lengths (typical dragonfly dimensions)
wingspan        = 1.0                      # (m) dragonfly wingspan





# ----------------- FLIGHT PARAMETERS ------------------------------------------
# Flapping kinematics
flap_amplitude  = 40.0                      # (degrees) flapping amplitude
flap_frequency  = 4.0                      # (Hz) flapping frequency
vehicle_velocity = 0.001                      # (m/s) forward flight velocity
angle_of_attack = 3.0                       # (degrees) vehicle angle of attack

# Flight conditions
rho             = 1.225                     # (kg/m^3) air density
mu              = 1.81e-5                   # (kg/ms) air dynamic viscosity

# Time parameters
flap_cycles     = 10                         # Number of flapping cycles to simulate
ttot            = flap_cycles / flap_frequency  # (s) total simulation time
nsteps          = 150 * flap_cycles         # Time steps (150 per cycle for efficiency)
dt              = ttot / nsteps             # (s) time step

# Simulation time window
tstart          = 0.00 * ttot               # (s) start time
tquit           = 1.00 * ttot               # (s) end time

Vinf(X, t)      = 2.0*[1.0, 0.0, 0.0]  # Freestream function





# ----------------- SOLVER PARAMETERS ------------------------------------------

# Aerodynamic solver
VehicleType     = uns.VLMVehicle           # Unsteady solver for flapping wings

# VPM particle shedding
p_per_step      = 4                         # Particles shed per time step (reduced for efficiency)
shed_starting   = true                     # No starting vortex needed
shed_unsteady   = true                      # Shed vorticity from unsteady loading
unsteady_shedcrit = 0.01                    # Circulation fluctuation threshold

# Regularization parameters
sigma_vlm_surf  = wingspan/50              # VLM-on-VPM smoothing radius
sigma_vlm_solver = 0.002
lambda_vpm      = 2.125                     # VPM core overlap
# sigma_vpm_overwrite = lambda_vpm * flap_frequency * flap_amplitude * pi/180 * dt / p_per_step #------------------is this correct----------------
sigma_vpm_overwrite = 0.1
sigmafactor_vpmonvlm = 1                    # VPM-on-VLM velocity factor

# Wing solver parameters
vlm_rlx         = 0.2                       # VLM relaxation factor
vlm_vortexsheet = false                     # Actuator surface model
vlm_vortexsheet_overlap = 2.125
vlm_vortexsheet_distribution = uns.g_pressure
vlm_vortexsheet_sigma_tbv = sigma_vpm_overwrite
vlm_vortexsheet_maxstaticparticle = vlm_vortexsheet ? 1000000 : nothing

# Force calculation parameters
KJforce_type    = "regular"                 # Kutta-Joukowski force type
include_trailingboundvortex = false
include_unsteadyforce = false                
add_unsteadyforce = false
include_parasiticdrag = false
add_skinfriction = false
calc_cd_from_cl = false
wing_polar_file = "xf-n0012-il-500000-n5.csv" 

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
fdom_prefix     = "dragonfly_fdom"               # Prefix of fluid domain files





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
RPMref = flap_frequency * 1                # Reference RPM equivalent

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


# ---------------------------------------------------------------REDEFINED SCRIPT OF UNS-------------------------------------------------------------------

function Vvpm_on_Xs(pfield::vpm.ParticleField, Xs::Array{T, 1}; static_particles_fun=(args...)->nothing, dt=0, fsgm=1) where {T}

    if length(Xs)!=0 && vpm.get_np(pfield)!=0
        # Omit freestream
        Uinf = pfield.Uinf
        pfield.Uinf = (t)->zeros(3)

        org_np = vpm.get_np(pfield)             # Original particles

        # Singularize particles to correct tip loss
        # NOTE: This doesn't include static particles, but there shouldn't be
        #       any in the field at this point anyways
        if abs(fsgm) != 1
            for P in vpm.iterator(pfield)
                P.sigma .*= fsgm
            end
        end

        # Add static particles
        static_particles_fun(pfield, pfield.t, dt)

        sta_np = vpm.get_np(pfield)             # Original + static particles

        # Add probes
        for X in Xs
            add_probe(pfield, X)
        end

        # Evaluate velocity field
        pfield.UJ(pfield)

        # Retrieve velocity at probes
        Vvpm = [Array(P.U) for P in vpm.iterator(pfield; start_i=sta_np+1)]

        # Remove static particles and probes
        for pi in vpm.get_np(pfield):-1:(org_np+1)
            vpm.remove_particle(pfield, pi)
        end

        # De-singularize particles
        if abs(fsgm) != 1
            for P in vpm.iterator(pfield)
                P.sigma ./= fsgm
            end
        end

        # Restore freestream
        pfield.Uinf = Uinf
    else
        Vvpm = [zeros(3) for i in 1:length(Xs)]
    end

    return Vvpm
end

add_probe(pfield::vpm.ParticleField, X) = vpm.add_particle(pfield, X, zeros(3), 1e-6; vol=0)



function generate_calc_aerodynamicforce(; add_parasiticdrag=false,
                                          add_skinfriction=false,
                                          airfoilpolar="xf-n0012-il-500000-n5.csv",
                                          parasiticdrag_args=(),
                                          )


    # Aerodynamic force from Kutta-Joukowski's theorem
    kuttajoukowski = generate_aerodynamicforce_kuttajoukowski("regular",
                                                                nothing, nothing,
                                                                false, nothing,
                                                                nothing, vehicle)

    if add_parasiticdrag
        parasiticdrag = generate_aerodynamicforce_parasiticdrag(airfoilpolar;
                                                                add_skinfriction=add_skinfriction,
                                                                parasiticdrag_args...)
    end

    # Force due to unsteady circulation
    unsteady(args...; optargs...) = calc_aerodynamicforce_unsteady(args...; add_to_Ftot=false, optargs...)


    function calc_aerodynamicforce(vlm_system, args...; per_unit_span=false, optargs...)

        # Delete any previous force field
        fieldname = per_unit_span ? "ftot" : "Ftot"
        if fieldname in keys(vlm_system.sol)
            pop!(vlm_system.sol, fieldname)
        end

        # Calculate unsteady force
        Ftot = unsteady(vlm_system, args...; per_unit_span=per_unit_span, optargs...)

        # Calculate Kutta-Joukowski force
        Ftot = kuttajoukowski(vlm_system, args...; per_unit_span=per_unit_span, optargs...)
        
        if add_parasiticdrag
            Ftot = parasiticdrag(vlm_system, args...;
                                    per_unit_span=per_unit_span, optargs...)
        end

        return Ftot
    end

    return calc_aerodynamicforce
end


function generate_aerodynamicforce_kuttajoukowski(KJforce_type::String,
                                sigma_vlm_surf, sigma_rotor_surf,
                                vlm_vortexsheet, vlm_vortexsheet_overlap,
                                vlm_vortexsheet_distribution,
                                vlm_vortexsheet_sigma_tbv;
                                vehicle=vehicle
                                )

    function calc_aerodynamicforce_kuttajoukowski(vlm_system::Union{vlm.Wing, vlm.WingSystem},
                                    prev_vlm_system, pfield, Vinf, dt, rho; t=0.0,
                                    per_unit_span=false, spandir=[0, 1, 0],
                                    include_trailingboundvortex=false,
                                    Vout=nothing, lenout=nothing,
                                    lencrit=-1, debug=false)

        m = vlm.get_m(vlm_system)    # Number of horseshoes



        # === START OF FIX ===
        if vehicle != nothing && !("Gamma" in keys(vlm_system.sol))
            # Get the global Gamma solution from the overall vehicle system
            global_gamma_solution = vehicle.system.sol["Gamma"]

            # Determine the starting index for the current vlm_system within the global Gamma array.
            idx_start = 1

            if vlm_system === vehicle.tilting_systems[1]  # fore_left
                idx_start = 1
            elseif vlm_system === vehicle.tilting_systems[2]  # fore_right
                m_fore_left = vlm.get_m(vehicle.tilting_systems[1])
                idx_start = m_fore_left + 1
            elseif vlm_system === vehicle.tilting_systems[3]  # hind_left
                m_fore_left = vlm.get_m(vehicle.tilting_systems[1])
                m_fore_right = vlm.get_m(vehicle.tilting_systems[2])
                idx_start = m_fore_left + m_fore_right + 1
            elseif vlm_system === vehicle.tilting_systems[4]  # hind_right
                m_fore_left = vlm.get_m(vehicle.tilting_systems[1])
                m_fore_right = vlm.get_m(vehicle.tilting_systems[2])
                m_hind_left = vlm.get_m(vehicle.tilting_systems[3])
                idx_start = m_fore_left + m_fore_right + m_hind_left + 1
            elseif vlm_system === vehicle.system
                # If it's the full system, it should already have Gamma
                @warn "Full vehicle system passed but Gamma is missing. This is unexpected."
                return nothing
            else
                # If for some reason vlm_system is not identified
                error("Cannot determine global Gamma slice for the given vlm_system. " *
                      "It is not identified as any of the four wing systems, and 'Gamma' is missing.")
            end

            # Calculate the end index for the slice.
            m = vlm.get_m(vlm_system)
            idx_end = idx_start + m - 1

            # Basic sanity check to prevent out-of-bounds access
            if idx_end > length(global_gamma_solution)
                error("Calculated Gamma slice indices ($idx_start to $idx_end) exceed the bounds of the global Gamma solution (length: $(length(global_gamma_solution))).")
            end

            # Manually assign the relevant slice of Gamma to the current vlm_system's solution dictionary.
            vlm_system.sol["Gamma"] = global_gamma_solution[idx_start:idx_end]

        elseif vehicle == nothing && !("Gamma" in keys(vlm_system.sol))
            error("Gamma solution not found in vlm_system.sol and the vehicle object is not provided for global Gamma access. Cannot proceed with Kutta-Joukowski calculation.")
        end
        # === END OF FIX ===



        # Nodes of every horseshoe
        Ap = uns._get_Xs(vlm_system, "Ap")
        A = uns._get_Xs(vlm_system, "A")
        B = uns._get_Xs(vlm_system, "B")
        Bp = uns._get_Xs(vlm_system, "Bp")

        # Midpoints of bound vortices
        ApA = (Ap .+ A)/2
        AB = (A .+ B)/2
        BBp = (B .+ Bp)/2

        # Evaluate VPM on each midpoint
        Xs = vcat(ApA, AB, BBp)

        # Evaluate VLM on each midpoint
        Vvlm = vlm.Vind.(Ref(vlm_system), Xs; t=t, ign_col=true, ign_infvortex=true)

        # Evaluate Vinf on each midpoint
        Vinfs = Vinf.(Xs, t)

        # Evaluate kinematic velocity on each node
        Vtran = uns._Vkinematic(vlm_system, prev_vlm_system, dt; t=t,
                                                    targetX=["Ap", "A", "B", "Bp"])

        # Pre-allocate memory
        Ftot = [zeros(3) for i in 1:m]
        if debug
            Vvpms, Vvlms, Vtranmids = ([zeros(3) for i in 1:m] for j in 1:3)
        end

        # Bound vortices to include in force calculation
        vortices = include_trailingboundvortex ? (1:3) : (2:2)  # 1==ApA, 2==AB, 3==BBp



        if !vlm_vortexsheet || KJforce_type=="regular"

            ## NOTE: Instead of calling the VPM, we use what was calculated
            ## by `solve()`, which includes Rotor-on-VLM velocities
            Vvpm = Vvpm_on_Xs(pfield, Xs; dt=dt)
            #Vvpm = vcat(vlm_system.sol["Vvpm_ApA"], vlm_system.sol["Vvpm_AB"], vlm_system.sol["Vvpm_BBp"])

        elseif vehicle == nothing

            error("KJ force calculation based on vortex sheet has been
                    requested, but vehicle has not been provided.")

        else

            # ----- Estimate VPM velocity (including rotors) on vortex sheet

            # Generate function of static particles
            static_particles_function = generate_static_particle_fun(pfield, pfield,
                                            vehicle,
                                            sigma_vlm_surf, sigma_rotor_surf;
                                            vlm_vortexsheet=vlm_vortexsheet,
                                            vlm_vortexsheet_overlap=vlm_vortexsheet_overlap,
                                            vlm_vortexsheet_distribution=vlm_vortexsheet_distribution,
                                            vlm_vortexsheet_sigma_tbv=vlm_vortexsheet_sigma_tbv)

            # Pre-calculate direction of lifting bound vortices
            BVdir = [zeros(3) for i in 1:m]

            for i in 1:m
                BVdir[i] .= B[i]
                BVdir[i] .-= A[i]
                BVdir[i] ./= norm(BVdir[i])
            end

            # Add static particles representing the vortex sheet
            org_np = vpm.get_np(pfield)
            static_particles_function(pfield)

            # Nullify strength of vortex sheet particles to avoid double accounting for VLM velocity
            for P in vpm.iterate(pfield; include_static=true)

                if P.static[1] && 1 <= P.index[1] <= m

                    # Save strength for future reference
                    for i in 1:3; P.M[i] = P.Gamma[i]; end;

                    # Save index before it gets overwritten by the FMM
                    P.vol[1] = P.index[1]

                    # Zero it out
                    P.Gamma .*= 1e-14
                end
            end

            # Evaluate velocity in entire particle field
            vpm._reset_particles(pfield)
            pfield.UJ(pfield)

            # Collect velocities
            weights = zeros(m)
            Vvpm = [zeros(3) for i in 1:m]

            for P in vpm.iterate(pfield; include_static=true)

                # if P.static[1] && 1 < P.index[1] <= m
                if P.static[1] && 1 <= P.vol[1] <= m

                    # ind = P.index[1]
                    ind = Int(P.vol[1])

                    # Calculate weight of this static particle as projection to lifting BV
                    weight = KJforce_type=="weighted" ? abs(dot(view(P.M, 1:3), BVdir[ind])) :
                             KJforce_type=="averaged" ? 1.0 :
                             error("Invalid KJforce_type. Options are `\"regular\"`, `\"weighted\"`, or `\"averaged\"`; got $(KJforce_type).")

                    # Add weighted velocity
                    for i in 1:3
                        Vvpm[ind][i] += weight*P.U[i]
                    end
                    weights[ind] += weight

                end
            end

            # Divide by total weight to get the final weighted velocity
            Vvpm ./= weights

            Vvpm = vcat(Vvpm, Vvpm, Vvpm)
            # ------
        end


        # Calculate KJ force
        for i in 1:m                 # Iterate over horseshoes
            for j in vortices        # Iterate over bound vortices (1==ApA, 2==AB, 3==BBp)

                # Bound vortex' length vector
                if j==1
                    l = (A[i]-Ap[i])
                elseif j==2
                    l = (B[i]-A[i])
                else
                    l = (Bp[i]-B[i])
                end

                # Bound vortex' midpoint
                # X = Xs[i + m*(j-1)]

                # Kinematic velocity at midpoint
                Vtranmid = (Vtran[i + m*(j-1)] + Vtran[i + m*j])/2

                # Effective velocity at midpoint
                V = Vvpm[i + m*(j-1)] + Vvlm[i + m*(j-1)] + Vinfs[i + m*(j-1)] + Vtranmid
                if Vout != nothing
                    push!(Vout, [Vvpm[i + m*(j-1)], Vvlm[i + m*(j-1)], Vinfs[i + m*(j-1)], Vtranmid])
                end

                # Circulation
                Gamma = vlm_system.sol["Gamma"][i]

                # Normalization factor
                if per_unit_span
                    len = abs((B[i][1]-A[i][1])*spandir[1] + (B[i][2]-A[i][2])*spandir[2] + (B[i][3]-A[i][3])*spandir[3])
                else
                    len = 1
                end
                if lenout != nothing
                    push!(lenout, per_unit_span==false || len>lencrit ? len : -10*lencrit)
                end

                # Kutta–Joukowski theorem: F = rho * V × vecGamma
                if per_unit_span==false || len > lencrit # NOTE: Use lencrit to avoid dividing by zero
                    Ftot[i][1] += rho * Gamma * (V[2]*l[3] - V[3]*l[2]) / len
                    Ftot[i][2] += rho * Gamma * (V[3]*l[1] - V[1]*l[3]) / len
                    Ftot[i][3] += rho * Gamma * (V[1]*l[2] - V[2]*l[1]) / len

                    if debug && j==2
                        Vtranmids[i] .+= Vtranmid
                        Vvpms[i] .+= Vvpm[i + m*(j-1)]
                        Vvlms[i] .+= Vvlm[i + m*(j-1)]
                    end
                end

            end
        end


        if vlm_vortexsheet && KJforce_type!="regular"
            # Remove static particles
            for pi in vpm.get_np(pfield):-1:org_np+1
                vpm.remove_particle(pfield, pi)
            end
        end

        if debug
            # Save Kutta-Joukowski force as a solution field
            vlm._addsolution(vlm_system, (per_unit_span ? "f" : "F")*"rk-vector", deepcopy(Ftot); t=t)
            vlm._addsolution(vlm_system, "Vtranmid-vector", Vtranmids; t=t)
            vlm._addsolution(vlm_system, "Vvpm-vector", Vvpms; t=t)
            vlm._addsolution(vlm_system, "Vvlm-vector", Vvlms; t=t)
        end

        # Add Kutta-Joukowski force to any existing force calculation
        fieldname = per_unit_span ? "ftot" : "Ftot"
        if fieldname in keys(vlm_system.sol)
            Ftot .+= vlm_system.sol[fieldname]
        end

        # Save total force as a solution field
        vlm._addsolution(vlm_system, fieldname, Ftot; t=t)

        return Ftot
    end

    return calc_aerodynamicforce_kuttajoukowski
end



# ----------------- 4) FORCE CALCULATION SETUP --------------------------------
println("\n4) Setting up force calculations...")

# Create unified force calculation function that matches monitor requirements
calc_aerodynamicforce_fun = generate_calc_aerodynamicforce(;
    add_parasiticdrag=include_parasiticdrag,
    add_skinfriction=add_skinfriction,
    airfoilpolar=wing_polar_file,
)

# Initialize forces array for simulation
forces = []

# Kutta-Joukowski force (primary lift mechanism) - simplified approach
kuttajoukowski = generate_aerodynamicforce_kuttajoukowski(
    "regular",            
    sigma_vlm_surf, 0.000000001,         # No rotor surface
    false, 2.125,                
    uns.g_pressure, sigma_vpm_overwrite;
    vehicle=vehicle,
)
push!(forces, kuttajoukowski)


# Combined force calculation function for simulation runtime
function simulation_calc_aerodynamicforce_fun(vlm_system, args...; per_unit_span=false, optargs...)
    fieldname = per_unit_span ? "ftot" : "Ftot"
    
    # Clear previous force calculations
    if fieldname in keys(vlm_system.sol)
        pop!(vlm_system.sol, fieldname)
    end
    
    Ftot = nothing
    
    # Apply each force component sequentially
    for (i, force) in enumerate(forces)
        try
            Ftot = force(vlm_system, args...; per_unit_span=per_unit_span, optargs...)
        catch e
            println("Warning: Force component $i failed: $e")
        end
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
# println("   Reynolds number: $(round(Re, digits=0))")

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
    # vpm_fmm = vpm.FMM(; p=4, ncrit=100, theta=0.4, phi=0.5),
    sigma_vlm_surf=sigma_vlm_surf,
    sigma_vlm_solver = sigma_vlm_solver, 
    sigma_rotor_surf=0.000000001,  # No rotors
    sigma_vpm_overwrite=sigma_vpm_overwrite,
    sigmafactor_vpmonvlm=sigmafactor_vpmonvlm,
    vlm_vortexsheet=vlm_vortexsheet,
    vlm_vortexsheet_overlap=vlm_vortexsheet_overlap,
    vlm_vortexsheet_distribution=vlm_vortexsheet_distribution,
    vlm_vortexsheet_sigma_tbv=vlm_vortexsheet_sigma_tbv,
    vlm_rlx=vlm_rlx,
    vlm_init=true,
    hubtiploss_correction=nothing,  # No rotors
    shed_starting=shed_starting,
    shed_unsteady=shed_unsteady,
    unsteady_shedcrit=unsteady_shedcrit,
    # extra_runtime_function=runtime_function,
    extra_runtime_function=monitors,
    
    # Output options
    save_path=save_path,
    run_name=run_name,
    save_wopwopin=false,
    # save_pfield=true,     # Required for fluid domain computation
    # save_staticpfield=true # Required for complete flow analysis
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
