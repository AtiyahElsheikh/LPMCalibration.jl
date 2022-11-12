"""
Initial draft of calibrating an LPM-based model
It is assumed that LoneParentsModel.jl is located in the following relative path ../LoneParentsModel.jl
This version works with LoneParentsModel.jl version TAG 0.6.1
                                                    =========
later version of LoneParentsModel is probably not compatible (simply switch to TAG-Version) 

This script can be directly executing as 
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
    
    include("./utilities.jl")
    const LPMPATH = "../LoneParentsModel.jl"
    addToLoadPath!(LPMPATH) 
    addToLoadPath!("$(LPMPATH)/src")
    addToLoadPath!("$(LPMPATH)/src/generic")

    using LPM: LPMVERSION 
    using LPM.ParamTypes: SimulationPars, seed!, DemographyPars
    using XAgents: alive, isFemale, isMale 

    # TODO some significant (simulation & model) parameters 
    #   could be made as flags but not necessarily all (e.g. matrices)
    #   to make it possible to load from an input file 
    function loadAndSeedSimulationPars() 
        simPars = SimulationPars(seed=0, finishTime = 2020 + 6 // 12) # default values  
        seed!(simPars)             
        simPars 
    end 

    function loadModelParameters()
        pars = DemographyPars() 
        # adhoc fix: todo improve
        pars.datapars.datadir =  "$(LPMPATH)/$(pars.datapars.datadir)"  
        pars.poppars.initialPop = 5000 
        pars 
    end

# Step II-  OK 
# ============
#       establish model definition to be used for calibration 
#       as well as simulation parameters 
#       every time the model is calibrated a deep copy instance is employed  

    module ModelDef
        using ..LPMLib: LPMPATH, loadModelParameters, addToLoadPath!, loadAndSeedSimulationPars, 
                        alive, isFemale, isMale
        include("$(LPMPATH)/mainHelpers.jl")

        #export MODPARS, 
        export SIMPARS, Model, setupModel, runModel!
        export alivePeople, females, yearsold  

        # simulation parameters used for all simulations 
        const SIMPARS =  loadAndSeedSimulationPars() 

        # help functions for fitness computations 
        alivePeople(model) = [person for person in model.pop if alive(person)] 
        females(model) = [fem for fem in alivePeople(model) if isFemale(fem)]
        males(model) = [m for m in alivePeople(model) if isMale(m)]
        # alive males of age a 
        males(model,a::Int) = [amale for amale in males(model) if yearsold(amale) == a]
        females(model,a::Int) = [afemale for afemale in females(model) if yearsold(afemale) == a]
        
        """
        population Pyramid of the available alive population within the model 
        return PP of males and females 
        """
        function populationPyramid(model)
            ppf = Vector{Int}(undef,91)
            ppm = Vector{Int}(undef,91); 
            for index in 1:90 
                ppm[index] = length(males(model,index - 1))    
                ppf[index] = length(females(model,index - 1))    
            end
            ppm[91] = length(males(model)) - sum(ppm[1:90])
            ppf[91] = length(females(model)) - sum(ppf[1:90])
            ppm,ppf 
        end
        
    end # ModelDefinition
    
end # LPMLib

@assert LPMLib.LPMVERSION == r"0.6.1"  

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
# generate a set of parameters for multiple simulations 
# TODO: can be merged with the simulation loop 

# TODO configuring using flags / input files 
const NUMSIMS = 10   # number of simulations 

setParValue!(parameter,activePar,val) = 
    setfield!( getfield(parameter,activePar.group), activePar.name, val  )


using Distributions: Uniform 
setRandParValue!(parameter,activePar) = 
    setParValue!(parameter,activePar, 
                    rand(Uniform(activePar.lowerbound,activePar.upperbound)))


# Step IV (OK)
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

const malePopPyramid2020 = 
    reverse!(dropdims(
                CSV.File("./data/202006PopulationPyramid.male.csv", header=0) |> Tables.matrix,
                dims=2)) 
@assert sum(malePopPyramid2020) == 33145709 
@assert length(malePopPyramid2020) == 91  

const femalePopPyramid2020 = 
    reverse!(dropdims(
                CSV.File("./data/202006PopulationPyramid.female.csv", header=0) |> Tables.matrix,
                dims = 2)) 
@assert sum(femalePopPyramid2020) == 33935525 
@assert length(femalePopPyramid2020) == 91 


# Step V (Ok / subject to validation)
# initially define a cost function (e.g. sum of least squares) 
# can be also a vector rather than a single value depending on the imperical data 
"
Compute population ratio for each age class
    assuming that arguments correspond to one dimensional matrix of identical lengths 
"  
function ppRatio(ppGender1,ppGender2)
    @assert( length(ppGender1) == length(ppGender2))
    res1 = Vector{Float64}(undef,length(ppGender1))
    res2 = Vector{Float64}(undef,length(ppGender1))
    npop = sum(ppGender1) + sum(ppGender2)
    res1 = ppGender1 / npop 
    res2 = ppGender2 / npop 
    res1,res2  
end 

const malePPRatio2020, femalePPRatio2020 = ppRatio(malePopPyramid2020,femalePopPyramid2020)

@assert  1.0 - eps() < sum(malePPRatio2020) + sum(femalePPRatio2020) < 1 + eps() 
@assert sum(malePopPyramid2020) /( sum(malePopPyramid2020) + sum(femalePopPyramid2020) ) - eps() <
            sum(malePPRatio2020) < 
            sum(malePopPyramid2020) /( sum(malePopPyramid2020) + sum(femalePopPyramid2020) ) + eps()

# TODO some optional plots / histograms

""" 
Assumption: 
    data and simulation data has the same length and 
    they correspond to the same year 
""" 
ppVecIndex(ppRatioData,ppRatioSim::Vector{Float64}) = 
    ( (ppRatioData .- ppRatioSim) ./ ppRatioData) .^ 2  

ppindex(ppVIndex1,ppVIndex2) = sum(ppVIndex1) + sum(ppVindex2)


# conduct calibration multiple simulation
# a. initially brute-force naive technique 
# b. Other suggested packatges, e.g. hypercube, differential evolution .. 
# ... 

using .LPMLib: loadModelParameters
using .LPMLib.ModelDef: SIMPARS,
            setupModel, runModel!, populationPyramid,
            alivePeople, females  

const parameters =  loadModelParameters()
const model = setupModel(parameters)

bestPPIndex = Inf
bestParSet = parameters 

# TODO merge parameter construction with fitness computations 
for index in 1:NUMSIMS 
    global bestPPIndex, bestParSet

    # new parameter set 
    for active in activePars 
        setRandParValue!(parameters,active)
    end
    # println(parameters[index].poppars.femaleAgeScaling)

    model = setupModel(parameters) 
    @time runModel!(model,SIMPARS,parameters)

    # Evaluate PP index for each run 
    #= TODO associate the following with Logging 
    npop     = length(alivePeople(model))
    nfemales = length(females(model))
    nmales   = npop - nfemales 
    
    println("# of pop: $(length(alivePeople(model))) out of which $(length(females(model))) females") 
    =# 

    malePPSim,femalePPSim = populationPyramid(model)
    femalePPRatioSim, malePPRatioSim = ppRatio(femalePPSim,malePPSim)
    malePPVIndex  =  ppVecIndex(malePPRatio2020,malePPRatioSim)
    femalePPVIndex = ppVecIndex(femalePPRatio2020, femalePPRatioSim)
    currppindex = sum(malePPVIndex) + sum(femalePPVIndex)

    println("Iteration $(index) : pp fitness: $(currppindex)")

    if currppindex < bestPPIndex
        bestPPIndex = currppindex
        bestParSet  = parameters
    end 
end

println("best parameter set obtained pp fitness of $bestPPIndex")


