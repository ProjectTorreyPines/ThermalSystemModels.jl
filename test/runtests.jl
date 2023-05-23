using ThermalSystem_Models
using Test

TSM = ThermalSystem_Models
TSMode = TSM.Dynamics
TSMdt = TSM.DeadTime

@testset "ThermalSystem_Models.jl" begin
    
    TSM = ThermalSystem_Models
    TSM.DeadTime.init_test()
    TSM.Dynamics.test_full_system_rankine()
    TSM.Dynamics.test_full_system_brayton()
end
