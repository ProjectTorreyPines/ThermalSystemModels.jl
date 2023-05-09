using .ThermalSystem_Models


function init_g()
TSM = ThermalSystem_Models
sys = TSM.rankine_with_loops()
run_sys(sys)
h=TSM.network2graph2(sys)
gdawg, Edict, Cdict, name_list = network2graph2(sys)
return gdawg, Edict, name_list
end

