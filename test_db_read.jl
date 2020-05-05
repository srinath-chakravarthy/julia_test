using DataFrames, CSV
using JSON
using Plots
using PyCall
si = pyimport("scipy.interpolate")
data = CSV.read("all_data_interpolated2_new.csv", copycols=true)
