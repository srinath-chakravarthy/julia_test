using CSV, DataFrames, Statistics
using Dierckx
using Plots
using LaTeXStrings
## --- Offset sinusoidal curves
dw = 2.0
dh = 12.0
offset = dw/10.0
t = LinRange(0, 1, 101)
x = dw .* (1 .- t)
y = (-dh/2.0) .* (1.0 .+ cos.(pi .* (1 .- t)))
dx = -dw
dy = -(dh/2.0) .* sin.(pi .* (1 .- t))
den = sqrt.(dx.^2 .+ dy.^2)
newx = x .+ offset .* dy ./ den
newy = y .- offset .* dx ./ den
