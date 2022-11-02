"""
Initial draft of calibrating an LPM-based model
""" 

include("./utilities.jl")

# Step I 
# adding paths to the LPM and other nonstandard packages 
# here it will be assumed that LPM.jl is located in the following relative path ../LPM.jl 
# ... 

# Step II 
# determine the active parameters of the LPM 
# i.e. which parameters are of interest and shall be variated from a simulation to another 
# a. initially this could be done manually as a proof of concept 
# b. later: active parameters to be set in from an input file with flags  

# Step III 
# Set an upper and lower bounds of every active parameter p, p_min and p_max   
# a. initially manually hard-coded 
# b. later: from an input file and flags etc. 

# Step IV 
# set up time-stamped simulation folder for placing (if needed) the results 
#   & placing the imperical data as well as other related stuffs 
#   & which model variables do they correspond to 
# a. initially hard-coded 
# b. later from input files / flags 

# Step V 
# solution-dependent 
# initially define a cost function (e.g. sum of least squares) 
# can be also a vector rather than a single value depending on the imperical data 

# Step VI 
# conduct calibration multiple simulation
# a. initially brute-force naive technique 
# b. Other suggested packatges, e.g. hypercube 





