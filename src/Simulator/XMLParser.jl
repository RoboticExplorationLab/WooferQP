include("../Common/Config.jl")

function ParseXML(woofer::WooferConfig; μ = 1.5, Δt = 0.0001)
###### ROBOT PARAMETERS #####

    ## Solver params ##
    woofer_timestep = Δt # timestep
    woofer_joint_solref = "0.001 1" # time constant and damping ratio for joints
    woofer_joint_solimp = "0.9 0.95 0.001" # joint constraint parameters

    woofer_geom_solref = "0.005 2" # time constant and damping ratio for geom contacts
    woofer_geom_solimp = "0.9 0.95 0.001" # geometry contact parameters

    woofer_armature = 0.0024 # armature for joints [kgm2]

    ## Geometry params ##
    woofer_thigh_length = woofer.geometry.upper_link_length / 2 # half length of thigh
    woofer_shin_length = woofer.geometry.lower_link_length / 2 # half length of shin
    woofer_shin_pos = 0.5 * ((2 * woofer_shin_length)^2 - (2 * woofer_thigh_length)^2)^(0.5)
    woofer_shin_angle = asin(woofer_thigh_length / woofer_shin_length)
    woofer_leg_radius = woofer.geometry.foot_radius # radius of leg capsule
    woofer_friction = μ # friction between legs and ground
    woofer_half_size = "$(woofer.geometry.body_length/2) $(woofer.geometry.body_width/2) $(woofer.geometry.body_height/2)" # half-size of body box

    woofer_start_position = "0 0 $(2*woofer_shin_pos + woofer_leg_radius)"# Initial position of the robot torso

    woofer_force_geom = "0 0 -0.34"

    ## Mass/Inertia Params ##
    woofer_frame_mass = woofer.inertial.frame_mass
    woofer_module_mass = woofer.inertial.module_mass
    woofer_shin_mass = woofer.inertial.lower_link_mass
    woofer_thigh_mass = woofer.inertial.upper_link_mass

	woofer_abduction_offset = woofer.geometry.abduction_offset

    woofer_frame_inertia = "0.0065733 0.074011 0.077763"
    woofer_module_inertia = "0.002449 0.005043 0.006616 -0.001784 -.00002 -0.000007"
    woofer_leg_inertia = "0.003575 0.006356 0.002973 -0.0001326 -0.0001079 -0.0002538"

    ## Joint params ##
    woofer_joint_range = "$(-woofer.actuator.max_joint_range) $(woofer.actuator.max_joint_range)"# joint range in rads for angular joints
    woofer_joint_force_range = "$(-woofer.actuator.max_joint_torque) $(woofer.actuator.max_joint_torque)" # force range for ab/ad and forward/back angular joints
    woofer_joint_damping = 0.2 # damping on ab/ad and f/b angular joints [Nm/rad/s]

    ## Sensor Noise Parameters ##
    woofer_accel_noise = 0.01
    woofer_encoder_noise = 0.001
    woofer_gyro_noise = 0.02
    woofer_encoder_vel_noise = 0.01
    woofer_force_noise = 0

    ## FILE PATHS  ##
    in_file = joinpath(@__DIR__, "woofer.xml")
    out_file = joinpath(@__DIR__, "woofer_out.xml")

    ## Parse the xml ##
    print("Parsing MuJoCo XML file:")
    print("Input xml: $(in_file)")
    print("Output xml: $(out_file)")

    open(in_file) do file
        filedata = read(file, String)

        ## Replace variable names with values ##
        # Solver specs
        filedata = replace(filedata, "woofer_timestep" => woofer_timestep)
        filedata = replace(filedata, "woofer_joint_solref" => woofer_joint_solref)
        filedata = replace(filedata, "woofer_geom_solref" => woofer_geom_solref)
        filedata = replace(filedata, "woofer_friction" => woofer_friction)
        filedata = replace(filedata, "woofer_armature" => woofer_armature)
        filedata = replace(filedata, "woofer_joint_solimp" => woofer_joint_solimp)
        filedata = replace(filedata, "woofer_geom_solimp" => woofer_geom_solimp)

        # Joint specs
        filedata = replace(filedata, "woofer_joint_range" => woofer_joint_range)
        filedata = replace(filedata, "woofer_joint_force_range" => woofer_joint_force_range)
        filedata = replace(filedata, "woofer_joint_damping" => woofer_joint_damping)

        # Geometry specs
        filedata = replace(filedata, "woofer_frame_mass" => woofer_frame_mass)
        filedata = replace(filedata, "woofer_module_mass" => woofer_module_mass)
        filedata = replace(filedata, "woofer_thigh_mass" => woofer_thigh_mass)
        filedata = replace(filedata, "woofer_shin_mass" => woofer_shin_mass)

        filedata = replace(filedata, "woofer_frame_inertia" => woofer_frame_inertia)
        filedata = replace(filedata, "woofer_module_inertia" => woofer_module_inertia)
        filedata = replace(filedata, "woofer_leg_inertia" => woofer_leg_inertia)
        filedata = replace(filedata, "woofer_leg_radius" => woofer_leg_radius)
        filedata = replace(filedata, "woofer_half_size" => woofer_half_size)
        filedata = replace(filedata, "woofer_leg_fb" => woofer.geometry.hip_center_x)
        filedata = replace(filedata, "woofer_leg_lr" => woofer.geometry.hip_center_y)
        filedata = replace(filedata, "woofer_start_position" => woofer_start_position)
        filedata = replace(filedata, "woofer_force_geom" => woofer_force_geom)
        filedata = replace(filedata, "woofer_thigh_length" => woofer_thigh_length)
        filedata = replace(filedata, "woofer_shin_length" => woofer_shin_length)
        filedata = replace(filedata, "woofer_shin_pos" => woofer_shin_pos)
        filedata = replace(filedata, "woofer_shin_angle" => woofer_shin_angle)
        filedata = replace(filedata, "woofer_abduction_offset" => woofer_abduction_offset)

        # Sensor noise
        filedata = replace(filedata, "woofer_accel_noise" => woofer_accel_noise)
        filedata = replace(filedata, "woofer_gyro_noise" => woofer_gyro_noise)
        filedata = replace(filedata, "woofer_encoder_noise" => woofer_encoder_noise)
        filedata = replace(filedata, "woofer_encoder_vel_noise" => woofer_encoder_vel_noise)
        filedata = replace(filedata, "woofer_force_noise" => woofer_force_noise)

        # Write the xml file
        open(out_file, "w") do f
            write(f, filedata)
        end
    end
end
