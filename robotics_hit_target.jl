using CSV, DataFrames, Statistics, Glob
using Dierckx
using Plots
using LaTeXStrings

# --- Simple to plot ballistic angle
x = collect(0.5:0.01:5)
y = 0.61
v = 6.78 # m/s
v2 = v^2
v4 = v^4
g = 9.81 # m/s^2
theta = []
# theta2 = []
# t1 = []
for x1 in x
    t = v4 - g * (g * x1^2 + 2 * v2 * y)
    if t > 0
        tt = atand((v2 + sqrt(t))/(g * x1))
        tt2 = atand((v2 - sqrt(t))/(g * x1))
    else
        tt = 0
        tt2 = 0
    end
    push!(theta, minimum([tt,tt2]))
end
# temp = similar(x)
# theta1 = similar(x)
# theta2 = similar(x)
#
# @. temp = sqrt(v4 - g * (g * x^2 + 2 * v2 * y))
# @. theta1 = atand((v2 + temp)/(g * x))
# @. theta2 = atand((v2 - temp)/(g * x))
