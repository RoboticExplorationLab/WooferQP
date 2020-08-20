using Revise
using StaticArrays, Rotations
import StaticArrays: SUnitRange

includet("../../src/QuadrupedDynamics.jl")
includet("../../src/MPCControl/MPCControl.jl")
includet("../../src/Utilities.jl")
includet("../../src/Config.jl")
includet("XMLParser.jl")

import .MPCControl
using MuJoCoSimulator

function run_sim()
    yaml_path = joinpath(@__DIR__, "Simulator.yaml")
    simulator_yaml = YAML.load(open(joinpath(yaml_path)))

    woofer = WooferConfig()
    ParseXML(woofer)
    xml_path = joinpath(@__DIR__, "woofer_out.xml")

    # Initialize simulator object
    control_timestep = simulator_yaml["simulator"]["sim_data_dt"]
    buffer_timestep = simulator_yaml["simulator"]["buffer_update_dt"]
    sim = MuJoCoSimulator.Simulator{Float64, Int64}(xml_path, 26, 12, control_timestep, buffer_timestep)
    
    param = MPCControl.ControllerParams(Float64, Int64)

    control! = let param = param
        function (torques, q, q̇, t)
            rot = UnitQuaternion(q[4], q[5], q[6], q[7])
            mrp = MRP(rot)
            ω = rot \ q̇[SUnitRange(4, 6)]

            x_est = [   q[SUnitRange(1, 3)]; 
                        Rotations.params(mrp); 
                        q̇[SUnitRange(1, 3)]; 
                        ω   ]

            # annoying way to get rid of knee joint measurements
            joint_pos = q[@SVector [8,9,11,13,14,16,18,19,21,23,24,26]]
            joint_vel = q̇[@SVector [7,8,10,12,13,15,17,18,20,22,23,25]]

            MPCControl.control!(torques, x_est, joint_pos, joint_vel, t, param)
        end
    end

    MuJoCoSimulator.simulate(sim, control!)
end

run_sim()