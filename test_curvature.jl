# Find curvature of top surface
using CSV, DataFrames, Statistics
using Dierckx
using Plots
using LaTeXStrings
## adfasdf
# cd("/home/srinath/Projects/brian_sheldon_expts/initial_expt/5um_flaw_peridoic/")
# cd("/home/srinath/data/Projects/brian_sheldon_expts/initial_expt/5um_flaw_peridoic/crack/5.0um")
workdir = "/home/srinath/data/Projects/moose_tests/supercomp_latest/curved_surface/glued/dw_2.0_dh_6.0/csv"
cd(workdir)
filebase = "full_model_quartz_top"
df_fin = DataFrame([Float64, Float64],[:time, :curvature])
# for i in collect(1:1:34)
pyplot()
p = plot(xlabel = "x" ,
        ylabel = "Curvature (1/m)",
        xmirror = false,
        framestyle = :box,
        legend = :outertopright,
        legendfontsize = 14,
        legendtitlefontsize = 16,
        tickfontsize = 14,
        guidefontsize = 16,
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        titlefontsize = 16,
        grid = false)

p2 = plot(xlabel = "x" ,
        ylabel = "Displacement",
        xmirror = false,
        framestyle = :box,
        legend = :outertopright,
        legendfontsize = 14,
        legendtitlefontsize = 16,
        tickfontsize = 14,
        guidefontsize = 16,
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        titlefontsize = 16,
        grid = false)
df2 = DataFrame(CSV.File("full_model.csv"))
for i in collect(4:1:83)
    filename = filebase * "_" * lpad(i-1,4,"0") *".csv"
    df = DataFrame!(CSV.File(filename))
    x = df[!,1].*1e-6 # Convert to m
    y = df[!,5].*1e-6 # Convert to m
    time1 = df2[i,1] # convert to hours
    spl = Spline1D(x,y, s= length(x), k = 3)
    xx = collect(0:0.001:5)
    xx .= xx .* 1e-6
    yy = evaluate(spl, xx)
    yp = derivative(spl, xx)
    ypp = derivative(spl, xx; nu = 2)
    curvature = ypp ./ (1.0 .+ yp .^2) .^(1.5)
    labstr = string(round(time1/100,digits=2))
    if (i % 10 == 0)
        plot!(p, xx.*1e6, curvature, label=labstr)
        plot!(p2, xx.*1e6, yy*1e9, label = labstr)
    end
    println(mean(curvature))
    push!(df_fin, [time1, mean(curvature)])
end
# workdir1 = "/home/srinath/data/Projects/moose_tests/supercomp_latest/curved_surface/glued/0.2mA/dw_0.1_dh_0.2/rst"
# cd(workdir1)
# for i in collect(1:1:200)
#     filename = filebase * "_" * string(i-1) *".csv"
#     df = DataFrame!(CSV.File(filename))
#     x = df[!,5].*1e-6 # Convert to m
#     y = df[!,3].*1e-6 # Convert to m
#     time1 = df[1,1] # convert to hours
#     spl = Spline1D(x,y, s= length(x), k = 5)
#     xx = collect(0:0.001:5)
#     xx .= xx .* 1e-6
#     yy = evaluate(spl, xx)
#     yp = derivative(spl, xx)
#     ypp = derivative(spl, xx; nu = 2)
#     curvature = ypp ./ (1.0 .+ yp .^2) .^(1.5)
#     labstr = string(time1/100)
#     # if (i % 10 == 0)
#     #     plot!(p, xx.*1e6, curvature, label=labstr)
#     #     plot!(p2, xx.*1e6, yy*1e6, label = labstr)
#     # end
#     println(mean(curvature))
#     push!(df_fin, [time1, mean(curvature)])
# end
# workdir1 = "/home/srinath/data/Projects/moose_tests/supercomp_latest/curved_surface/glued/0.04mA/dw_0.1_dh_0.2/rst"
# cd(workdir1)
# for i in collect(2:1:36)
#     filename = filebase * "_" * string(i-1) *".csv"
#     df = DataFrame!(CSV.File(filename))
#     x = df[!,5].*1e-6 # Convert to m
#     y = df[!,3].*1e-6 # Convert to m
#     time1 = df[1,1] # convert to hours
#     spl = Spline1D(x,y, s= length(x), k = 5)
#     xx = collect(0:0.001:5)
#     xx .= xx .* 1e-6
#     yy = evaluate(spl, xx)
#     yp = derivative(spl, xx)
#     ypp = derivative(spl, xx; nu = 2)
#     curvature = ypp ./ (1.0 .+ yp .^2) .^(1.5)
#     labstr = string(time1/100)
#     # if (i % 2 == 0)
#     #     plot!(p, xx.*1e6, curvature, label=labstr)
#     #     plot!(p2, xx.*1e6, yy*1e6, label = labstr)
#     # end
#     println(time1, mean(curvature))
#     push!(df_fin, [time1, mean(curvature)])
# end

png(p,"curvature_space_time.png")
png(p2, "uy_space_time.png")
cd(workdir)
## a
p = plot(xlabel = "Time (hrs)" ,
        ylabel = "Curvature (1/m)",
        xmirror = false,
        framestyle = :box,
        legend = :outertopright,
        legendfontsize = 14,
        legendtitlefontsize = 16,
        tickfontsize = 14,
        guidefontsize = 16,
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        titlefontsize = 16,
        grid = false)

plot!(df_fin[!,1], -df_fin[!,2])
png(p,"curvature2.png")

## asdf
p = plot(xlabel = "Position" ,
        ylabel = "Rct",
        xmirror = false,
        framestyle = :box,
        tickfontsize = 14,
        guidefontsize = 16,
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        titlefontsize = 16,
        grid = false)

f(x) = (0.5e-9*cos(2.0*3.14157 * x/10.0) + 0.5e-9 + 1e-11)*1e9
# f(x) = (1.0/(1.0 + exp(-2.0*2.0*(x-100.0)))*(1000-10) + 10)
plot!(f, -200,200)
png(p,"Rct.png")
## Stress vs time
df = DataFrame(CSV.File("stress_time.csv"))
p = plot(xlabel = "Time" ,
        ylabel = "Stress (MPa)",
        xmirror = false,
        framestyle = :box,
        tickfontsize = 14,
        guidefontsize = 16,
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        titlefontsize = 16,
        grid = false)

plot!(df[!,:Time]./100, df[!,2])
png(p, "stress_time.csv")
