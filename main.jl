"""
Initial draft of calibrating an LPM-based model
It is assumed that LoneParentsModel.jl is located in the following relative path ../LoneParentsModel.jl
This version works with LoneParentsModel.jl version TAG-Verion 0.6.0
                                                    =================
later version of LoneParentsModel may be not compatible (simply switch to TAG-Version) 
can be directly executing as 
# julia main.jl 
or within REPL as 
julia> include("main.jl")
""" 


# Step I - OK 
# ===========  
#       adding paths to the LPM and other nonstandard packages 
#

# Viewing LPM as a library    
module LPMLib
    # export VERSION
    export LPMPATH 
    export loadAndSeedSimulationPars, loadModelParameters

    # const VERSION = r"V0.6.0"   
    
    include("./utilities.jl")
    const LPMPATH = "../LoneParentsModel.jl"
    addToLoadPath!(LPMPATH) 
    addToLoadPath!("$(LPMPATH)/src")
    addToLoadPath!("$(LPMPATH)/src/generic")

    using LPM.ParamTypes: SimulationPars, seed!, DemographyPars

    # TODO some significant (simulation & model) parameters 
    #   could be made as flags but not necessarily all (e.g. matrices)
    #   to make it possible to load from an input file 
    function loadAndSeedSimulationPars() 
        simPars = SimulationPars() # default values  
        simPars.seed = 0           # Every execution will be different 
        seed!(simPars)             
        simPars 
    end 

    function loadModelParameters()
        pars = DemographyPars() 
        # adhoc fix: todo improve
        pars.datapars.datadir =  "$(LPMPATH)/$(pars.datapars.datadir)"  
        pars.poppars.initialPop = 500 
        pars 
    end

# Step II-  OK 
# ============
#       establish model definition to be used for calibration 
#       as well as simulation parameters 
#       every time the model is calibrated a deep copy instance is employed  

    module ModelDef
        using ..LPMLib: LPMPATH, loadModelParameters, addToLoadPath!, loadAndSeedSimulationPars
        include("$(LPMPATH)/mainHelpers.jl")

        export MODEL, MODPARS, SIMPARS, Model, setupModel 

        # simulation parameters defining the 
        const SIMPARS =  loadAndSeedSimulationPars() 

        # model parameter definitions that will be variated from a simulation to another 
        const MODPARS  = loadModelParameters()

        # This is a model definition that will be deep copied for every simulation instance 
        const MODEL = setupModel(MODPARS) 
    end # ModelDefinition
    
end # LPMLib

# if the above package is moved to LoneParentsModel.jl 
# @assert LoneParentsModel.VERSION == r"V0.6.0"  


# Step III 
# ========
#       determine the active parameters of the LPM , their lower bounds and upper bouns 
#       i.e. which parameters are of interest and shall be variated from a simulation to another 
#       a. initially this could be done manually as a proof of concept (OK)
#       b. later: active parameters to be set in from an input file with flags  

mutable struct ActiveParameter{ValType} 
    lowerbound::ValType 
    upperbound::ValType  
    group::Symbol
    name::Symbol         
    
    function ActiveParameter{ValType}(low,upp,gr,id) where ValType 
        @assert low <= upp 
        new(low,upp,gr,id)
    end
end 

# TODO Base.show methods 

# as a hint @Umberto
# The following is how you can find out the group names and parameter names 
# from REPL 
# julia> include("main.jl")
# julia> using .LPMLib: DemographyPars
# julia> using .LPMLib.ModelDef: MODPARS 
# julia> println(fieldnames(DemographyPars))
# julia> println(fieldnames(typeof(MODPARS.poppars))) 
#  

# adhoc choosing 7 active parameters 
# TODO: choices and values to be loaded from a file 
const femaleAgeScaling      = ActiveParameter{Float64}(10.0,20,:poppars,:femaleAgeScaling) 
const maleAgeScaling        = ActiveParameter{Float64}(5.0,15.0,:poppars,:maleAgeScaling)

const femaleMortalityBias   = ActiveParameter{Float64}(0.77,0.89,:poppars,:femaleMortalityBias)
const maleMortalityBias     = ActiveParameter{Float64}(0.71,0.83,:poppars,:maleMortalityBias)

const variableDivorce       = ActiveParameter{Float64}(0.01,0.14,:divorcepars,:variableDivorce)
const basicDivorceRate      = ActiveParameter{Float64}(0.01,0.15,:divorcepars,:basicDivorceRate)

const basicMaleMarraigeProb = ActiveParameter{Float64}(0.5,0.9,:marriagepars,:basicMaleMarriageProb)

# TODO: to add meaningful constraints, e.g. 
#       femaleAgeScaling > maleAgeScaling + something
#       maleMortalityBias < femaleMortalityBias - something  
 
const activePars = [femaleAgeScaling, maleAgeScaling, femaleMortalityBias, maleMortalityBias,
                    variableDivorce, basicDivorceRate, 
                    basicMaleMarraigeProb]


# Step III 
# ========
# generate a set of models and parameters for multiple simulations 

# TODO configuring using flags / input files 
const NUMSIMS = 100    # number of simulations 

function setParValueExp(pname::String, active::ActiveParameter{ValType}, val::ValType) where ValType  
    @assert val >= active.lowerbound && val <= active.upperbound
    gr = active.group 
    id = active.name  
    Meta.parse("$(pname).$(gr).$(id) = $(val)")
end 

using Distributions: Uniform 
function setRandParValueExp(pname,active) 
    val = rand(Uniform(active.lowerbound,active.upperbound))
    setParValueExp(pname,active,val)
end

# TODO subject to improvement, to make the first argument actual rather than string
setRandParValue!(pname,active) = eval(setRandParValueExp(pname,active))

# establish model parameter sets according to active parameters 

using .LPMLib.ModelDef: MODPARS, setupModel, Model 

# TODO establish a function that generates random parameters  
#   need to take care that this makes eval(exp) become local to a function 

using .LPMLib: DemographyPars, SimulationPars

const parameters =  DemographyPars[] 
# TODO No need to have an array of models if model construction is merged with simulation       
const models     =  Model[]              
# const simpars    =  SimulationPars[]   # probably not needed, it is the same across all simulations 

# TODO : this is now only a proof of concept. 
#        it makes sense to merge this loop with the simulation for loop
for index in 1:NUMSIMS 
    push!(parameters,deepcopy(MODPARS)) 
    for active in activePars 
        setRandParValue!("parameters[$index]", active)
    end
    # TODO: deepcopy could be fragile, in this case, model to be constructed here
    mod = setupModel(parameters[index]) # instead of uncertain deepcopy(MODEL)
    push!(models,mod) 
    #println(parameters[index].poppars.femaleAgeScaling)
end

 
# Step IV
# =======
# Establish / load data to which the model definition is going to be calibrated 
#   & placing the imperical data as well as other related stuffs (OK)
#   & which model variables do they correspond to (N.A.)
# a. initially hard-coded 
# b. initially Population Pyramid data 
# c. later from input files / flags and other data & fitness indices  


# loading data/202006PopulationPyramid.csv
#       columns correspnds to male numbers vs. female numbers 
#       last row corresponds to population number of age 0 
#       2nd  row corresponds to population number of age 89
#       1st  row corresponds to population number of age 90+ 
#       the male (female) population should sum to 33145709	(33935525) 

using CSV 
using Tables 

const malePopPyramid2020 = reverse!(CSV.File("./data/202006PopulationPyramid.male.csv", header=0) |> Tables.matrix) 
@assert sum(malePopPyramid2020) == 33145709 
@assert size(malePopPyramid2020) == (91,1)  

const femalePopPyramid2020 = reverse!(CSV.File("./data/202006PopulationPyramid.female.csv", header=0) |> Tables.matrix) 
@assert sum(femalePopPyramid2020) == 33935525 
@assert size(femalePopPyramid2020) == (91,1) 


# Step V
# initially define a cost function (e.g. sum of least squares) 
# can be also a vector rather than a single value depending on the imperical data 

"
Compute population ratio for each age class
    assuming that arguments correspond to one dimensional matrix of identical lengths 
"  
function computePPRatio(agemale,agefemale)
    @assert( length(agefemale) == length(agemale))
    res = Vector{Float64}(undef,length(agemale))
    npop = sum(agefemale) + sum(agemale)
    for index in 1:length(agemale) 
       res[index] =  agemale[index] / npop 
    end
    res 
end 

malePPRatio = computePPRatio(malePopPyramid2020,femalePopPyramid2020)

@assert sum(malePPRatio) < 
    sum(malePopPyramid2020) /( sum(malePopPyramid2020) + sum(femalePopPyramid2020) ) + eps()

# Assumption: 
# data corresponds to year 2020-07
# model is simulated till  2020-07 
function evaluatePPIndex(ppRatio1,ppRatio2::Vector{Float64}) 
    # compute Population pyramid 
end 

# Step VI 
# conduct calibration multiple simulation
# a. initially brute-force naive technique 
# b. Other suggested packatges, e.g. hypercube, differential evolution .. 
# ... 





