# Find curvature of top surface
using CSV, DataFrames, Statistics
using Dierckx
cd("/home/srinath/Projects/brian_sheldon_expts/initial_expt/100um_flaw")
filebase = "test_uy_time"
df_fin = DataFrame([Float64, Float64],[:time, :curvature])
for i in collect(1:1:159)
    filename = filebase * "_" * string(i-1) *".csv"
    df = CSV.read(filename)
    x = df[!,5].*1e-6 # Convert to m
    y = df[!,3].*1e-6 # Convert to m
    time = df[1,1]/3600 # convert to hours
    spl = Spline1D(x,y)
    xx = collect(-200:0.001:200)
    xx .= xx .* 1e-6
    yy = evaluate(spl, xx)
    yp = derivative(spl, xx)
    ypp = derivative(spl, xx; nu = 2)
    curvature = ypp ./ (1.0 .+ yp .^2) .^(1.5)
    println(mean(curvature))
    push!(df_fin, [time, mean(curvature)])
end
