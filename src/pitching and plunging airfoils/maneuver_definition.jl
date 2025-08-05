#--------------------------------------------------------------------------------NEERAJ BALACHANDAR-----------------------------------------------------------------------------------------------------
__authors__ = "NEERAJ BALACHANDAR"
__contact__ = "neerajbalachandar@gmail.com"

"""
    wing_maneuver(; disp_plot = true, add_wings = true, vehicle_velocity::Real=0.001, 
                    reduced_frequency::Real=4.0, h_amplitude::Real=0.375, 
                    pitch_amplitude::Real=15.0, phase_offset::Real=pi/2)

# ARGUMENTS
* `reduced_frequency::Real`: Non-dimensional frequency k = 2piyfc/U∞
* `h_amplitude::Real`: Non-dimensional plunge amplitude h = h₀/c  
* `pitch_amplitude::Real`: Pitch amplitude θ₀ in degrees
* `phase_offset::Real`: Phase offset φ between plunge and pitch motions (radians)
* `vehicle_velocity::Real`: Forward flight velocity (normalized)

The kinematic equations follow:
* Plunge motion: y(t) = -h₀ cos(kt)  
* Pitch motion: θ(t) = -θ₀ cos(kt + φ)

Where t is non-dimensionalized time from 0 to 1.
"""

function wing_maneuver(; disp_plot = true, 
                         add_wings = true, 
                         vehicle_velocity::Real=0.0,
                         angle_of_attack::Real=0.0)


    chord = 1.50
    wingspan = 8.25

    reduced_frequency=3.0
    h_amplitude=0.10
    pitch_amplitude=15.0
    phase_offset=π/2

    freq = (reduced_frequency)/(2*pi*chord)
    plunge_amp = h_amplitude*chord
    pitch_amp = pitch_amplitude * (π/180)


    vehicle_velocity_func(t) = [vehicle_velocity, 0, 0]
    vehicle_angle_func(t) = [0, angle_of_attack, 0]  


    # Left wing motion (assuming symmetric flapping)
    left_wing_angle(t) = begin
    # Plunge motion affects the "flapping" angle (Ax rotation)
    plunge_angle_left = plunge_amp * sin(freq * 2 * pi * t * ttot) # Ay rotation for plunge

    # Pitch motion (Az rotation around span axis)
    pitch_angle_left = -pitch_amp * sin(freq * 2pi * t * ttot + phase_offset)

    # Return [Ax, Ay, Az] in degrees
    [plunge_angle_left * (180/pi), pitch_angle_left * (180/pi), 0.0]
    end

    # Right wing motion (symmetric)
    right_wing_angle(t) = begin
    # Plunge motion (opposite direction for symmetric flapping)
    plunge_angle_right = -plunge_amp * sin(freq * 2 * pi * t * ttot)

    # Pitch motion (same as left wing)
    pitch_angle_right = -pitch_amp * sin(freq * 2 * pi * t * ttot + phase_offset)

    # Return [Ax, Ay, Az] in degrees
    [plunge_angle_right * (180/pi), pitch_angle_right * (180/pi), 0.0]
    end
    
    
    wing_angles = (left_wing_angle, right_wing_angle)
    rotor_rpms = ()

    maneuver = uns.KinematicManeuver(wing_angles, rotor_rpms, 
                                     vehicle_velocity_func, vehicle_angle_func)
    
    if disp_plot
        uns.plot_maneuver(maneuver)
    end
    
    return maneuver
end










    

    # function asymmetric_phase(t, downstroke_fraction)
    #     if t < downstroke_fraction
    #         return (t / downstroke_fraction) * pi
    #     else
    #         return pi + ((t - downstroke_fraction) / (1 - downstroke_fraction)) * pi
    #     end
    # end

    # downstroke_fraction = 0.4

    # function left_wing_angle(t)
    #     t_dim = t * ttot  # Dimensionalized time
    #     t_cycle = mod(freq * t_dim, 1.0)
    #     phase = asymmetric_phase(t_cycle, downstroke_fraction)
    #     plunge_angle_left = plunge_amp * cos(phase)
    #     pitch_angle_left = -pitch_amp * cos(phase + phase_offset)
    #     [plunge_angle_left * (180/pi), pitch_angle_left * (180/pi), 0.0]
    # end

    # function right_wing_angle(t)
    #     t_dim = t * ttot
    #     t_cycle = mod(freq * t_dim, 1.0)
    #     phase = asymmetric_phase(t_cycle, downstroke_fraction)
    #     plunge_angle_right = -plunge_amp * cos(phase)
    #     pitch_angle_right = -pitch_amp * cos(phase + phase_offset)
    #     [plunge_angle_right * (180/pi), pitch_angle_right * (180/pi), 0.0]
    # end


#     function wing_maneuver(; disp_plot=true, 
#                          add_wings=true, 
#                          vehicle_velocity::Real=0.0,
#                          angle_of_attack::Real=0.0,
#                          ttot::Real = 1.0)  # ttot should be passed in

#     chord = 1.50
#     wingspan = 8.25

#     reduced_frequency = 3.0
#     h_amplitude = 0.10
#     pitch_amplitude = 15.0
#     phase_offset = pi/2

#     freq = reduced_frequency / (2 * pi * chord)
#     plunge_amp = h_amplitude * chord
#     pitch_amp = pitch_amplitude * (pi / 180)  # degrees to radians

#     vehicle_velocity_func(t) = [vehicle_velocity, 0, 0]
#     vehicle_angle_func(t) = [0, angle_of_attack, 0]


#     function left_wing_angle(t)
#         theta = 2*pi * freq * t * ttot
#         phase = mod(theta, 2*pi)
#         plunge = plunge_amp * sin(theta)
#         if phase < pi
#             pitch = 0.0
#         else
#             pitch = -pitch_amp * sin(theta + phase_offset)
#         end
#         return [plunge * (180/pi), pitch * (180/pi), 0.0]
#     end

#     function right_wing_angle(t)
#         theta = 2*pi * freq * t * ttot
#         phase = mod(theta, 2*pi)
#         plunge = -plunge_amp * sin(theta)
#         if phase < pi
#             pitch = 0.0
#         else
#             pitch = -pitch_amp * sin(theta + phase_offset)
#         end
#         return [plunge * (180/pi), pitch * (180/pi), 0.0]
#     end





#     wing_angles = (left_wing_angle, right_wing_angle)
#     rotor_rpms = ()
#     maneuver = uns.KinematicManeuver(wing_angles, rotor_rpms,
#                                      vehicle_velocity_func, vehicle_angle_func)

#     if disp_plot
#         uns.plot_maneuver(maneuver)
#     end
#     return maneuver
# end
