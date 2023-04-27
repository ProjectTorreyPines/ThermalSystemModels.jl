#===============================================================#
#               COMPRESSION COMPONENTS                          #
#   Compressors
#   Pumps
#
#===============================================================#
#                    EQUATIONS
#   V̄ = velocity
#   v = specific volume
#   V = volume
#
#   e = P/ρ + V̄^2 / 2 + gz  =>  u + ke + pe
#   E = m * e
#   Ė = ṁ * e
#
#   Wflow = pressure * volume = Pv
#   Flow total energy
#        θ = Pv + e
#           =>  (Pv + u) + (V^2/2) + (gz)    
#           =>  h + ke + pe
#           =>  h + V^2 / 2 + gz 
#
#   Conservation of Mass => always applies
#
#       ṁ_in = ṁ_out
#
#   Conservation of energy
#   
#       Ėin     = Q̇in  + Ẇin  + ∑ (ṁ_in  * θ_in)  
#       Ėout    = Q̇out + Ẇout + ∑ (ṁ_out * θ_out)
#
#   Steady State ΔĖsys = 0
#                       Ėin = Ėin
#           (Q̇in - Q̇out) + (Ẇin - Ẇout) = ∑ (ṁ_out * θ_out) - ∑ (ṁ_in  * θ_in) 
#   SIGN CONVENTIONS
#       Q, W 
#       
#    Always true
#       Ėin - Ėout = ΔĖsys
#===============================================================#
using Graphs

abstract type thermal_system    end 

abstract type model_component     <: Thermal_System end       # "Blocks"
abstract type fluid_component     <: model_component end      # functinoal element - pump, compressor, turbine, boiler ...
abstract type source_element      <: model_component end      # Flow source, heat sources
abstract type reservoir           <: model_component end      # heat and work reservoirs

abstract type model_nodes            <: thermal_system   end
abstract type fluid_nodes            <: model_nodes      end
abstract type work_nodes             <: model_nodes      end
abstract type heat_nodes            <: model_nodes  end


abstract type fluid <: thermal_system  end
abstract type component_parameters end
abstract type external_interface <: model_component end     # for tracking heat and work loads between cycle

mutable struct no_connector <: connector end

#ISENTROPIC RELATIONSHIPS - for creating wrapper functions in the other files
ideal_gas_Δu(cv,T2,T1) = cv * (T2 - T1)
ideal_gas_Δh(cp,T2,T1) = cp * (T2 - T1)









