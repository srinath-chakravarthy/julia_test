# Find curvature of top surface
using CSV, DataFrames, Statistics, Glob
using Dierckx
using Plots
using LaTeXStrings
## Load all files
# cd("/home/srinath/Projects/brian_sheldon_expts/initial_expt/5um_flaw_peridoic/")
# cd("/home/srinath/data/Projects/brian_sheldon_expts/initial_expt/5um_flaw_peridoic/crack/5.0um")
workdir = "/home/srinath/data/Projects/moose_tests/supercomp_latest/curved_surface/glued/new_simulations/1mA/dw_0.5_dh_0.5/csv"
cd(workdir)
filebase1 = "full_model_ceramic_surface"
# filebase2 = "full_model_interface"
filebase3 = "full_model_quartz_top"
filebase4 = "full_model_ceramic_bot"

files_ceramic_top = glob(filebase1 * "_*csv", workdir)
files_ceramic_bot = glob(filebase4 * "_*csv", workdir)
# files_interface = glob(filebase2 * "_*csv", workdir)
files_quartz_top  = glob(filebase3 * "_*csv", workdir)


post_proc_file = "full_model.csv"
dfs_ceramic_top = DataFrame.(CSV.File.(files_ceramic_top))
# dfs_interface = DataFrame.(CSV.File.(files_interface))
dfs_quartz = DataFrame.(CSV.File.(files_quartz_top))
dfs_ceramic_bot = DataFrame.(CSV.File.(files_ceramic_top))
df_post = DataFrame(CSV.File(post_proc_file))

## --- now sync all of these and compute and plot curvature

df_fin = DataFrame([Float64, Float64],[:time, :curvature])

p = plot(xlabel = L"x $(\mu m)$" ,
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

p2 = plot(xlabel = L"x $(\mu m)$" ,
        ylabel = L"uy (nm)",
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


for i in 2:length(dfs_quartz)
    dfs_quartz[i][!,:time] .= df_post[i,:time]
    x = dfs_quartz[i][!,"x"].*1e-6 # Convert to m
    y = dfs_quartz[i][!,"uy"].*1e-6 # Convert to m
    time1 = df_post[i,1] # convert to hours
    spl_top = Spline1D(x,y, s= length(x), k = 3)
    # x = dfs_ceramic_bot[i][!,"x"].*1e-6 # Convert to m
    # y = dfs_ceramic_bot[i][!,"uy"].*1e-6 # Convert to m
    time1 = df_post[i,1] # convert to hours
    # spl_bot = Spline1D(x,y, s= length(x), k = 3)


    xx = collect(0:0.001:5)
    xx .= xx .* 1e-6
    yy_top = evaluate(spl_top, xx)
    # yy_bot = evaluate(spl_bot, xx)
    yy = yy_top;
    spl = Spline1D(xx, yy, s= length(x), k = 3)
    yp = derivative(spl, xx)
    ypp = derivative(spl, xx; nu = 2)
    curvature = ypp ./ (1.0 .+ yp .^2) .^(1.5)
    labstr = string(round(time1,digits=2))
    if (i  == 10 || i == 15 || i == 20 || i == 40 || i == 100)
        plot!(p, xx.*1e6, -curvature, legend = :false)
        plot!(p2, xx.*1e6, abs.(yy*1e9), legend = :false)
        # plot!(p2, x.*1e6, y.*1e6)
    end
    println(time1, mean(curvature))
    push!(df_fin, [time1, mean(curvature)])
end
png(p,"curvature_space_time2.png")
png(p2, "uy_space_time2.png")
cd(workdir)
CSV.write("curvature_time2.csv", df_fin)

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

plot!(df_fin[!,1]/3600, -df_fin[!,2])
png(p,"curvature4.png")
# ## ----
# p = plot(xlabel = L"$\bar{x} = x/dw$" ,
#         ylabel = L"$\eta/\eta_0$",
#         xmirror = false,
#         framestyle = :box,
#         legend = :false,
#         legendfontsize = 14,
#         legendtitlefontsize = 16,
#         tickfontsize = 14,
#         guidefontsize = 16,
#         foreground_color_legend=nothing,
#         background_color_legend=nothing,
#         titlefontsize = 16,
#         grid = false)
#         df_new = dfs_quartz[40]

# plot!(p, df_new[!,:x]./2.0, df_new[!,:li_ion_V]./1.0)
# vline!(p, [1.0])
# png(p, "plating_overpotential.png")
