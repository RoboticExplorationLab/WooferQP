using StaticArrays

include("../Common/QuadrupedDynamics.jl")
include("../Common/MPCControl/MPCControl.jl")
include("../Common/Utilities.jl")
include("../Common/Config.jl")
include("XMLParser.jl")

import .MPCControl
using MuJoCoSimulator

function run_sim()
    yaml_path = joinpath(@__DIR__, "Simulator.yaml")
    simulator_yaml = YAML.load(open(joinpath(yaml_path)))

    woofer = WooferConfig()
    ParseXML(woofer)
    xml_path = yaml_path = joinpath(@__DIR__, "woofer_out.xml")

    # Initialize simulator object
    control_timestep = simulator_yaml["simulator"]["sim_data_dt"]
    buffer_timestep = simulator_yaml["simulator"]["buffer_update_dt"]
    sim = MuJoCoSimulator.Simulator{Float64, Int64}(xml_path, 26, 12, control_timestep, buffer_timestep)
    
    param = MPCControl.ControllerParams(Float64, Int64)

    control! = let param = param
        function (torques, q, q̇, t)
            MPCControl.control!(torques, q, q̇, t, param)
        end
    end

    MuJoCoSimulator.simulate(sim, control!)
end

run_sim()