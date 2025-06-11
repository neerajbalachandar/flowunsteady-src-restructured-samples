#--------------------------------------------------------------------------------NEERAJ BALACHANDAR-----------------------------------------------------------------------------------------------------
__authors__ = "NEERAJ BALACHANDAR"
__contact__ = "neerajbalachandar@gmail.com"

#--------------WILL BE UPDATED BASED ON SOLVING THE VEHICLE KINEMATIC EQUATIONS TO COUPLE THE CONTROL INPUT WITH KINEMATICS, EITHER BY EXPERIMENTAL VALIDATION FOR HOVERING CONDITION AND FORWARD FLIGHT OR THROUGH SIMULATION-----
#-----------------------------------------------------------FIXED VEHICLE WITH SYNCHRONISED PHASE FLAPPING OF SAME AMPLITUDE AND FREQUENCY FOR EACH WING----------------------

function generate_dragonfly_maneuver(; vehicle_velocity::Real=0.0, angle_of_attack::Real=0.0)
    
    amp = 60 * (pi/180) 
    freq = 30 # Hz      
    
    vehicle_velocity_func(t) = [vehicle_velocity, 0, 0]
    vehicle_angle_func(t) = [0, angle_of_attack, 0]
    
    fore_wing_left_angle(t) = [0, amp * sin(freq * t), 0]
    fore_wing_right_angle(t) = [0, -amp * sin(freq * t), 0]
    hind_wing_left_angle(t) = [0, amp * sin(freq * t), 0]
    hind_wing_right_angle(t) = [0, -amp * sin(freq * t), 0]
    
    wing_angles = (
        fore_wing_left_angle,
        fore_wing_right_angle,
        hind_wing_left_angle,
        hind_wing_right_angle
    )
  
    rotor_rpms = ()
    
    # kinematic maneuver
    maneuver = uns.KinematicManeuver(
        wing_angles, 
        rotor_rpms, 
        vehicle_velocity_func, 
        vehicle_angle_func
    )
    return maneuver
end
