# Find curvature of top surface
using CSV, DataFrames, Statistics, Glob
using Dierckx
using Plots
using LaTeXstringings
## --- Now load each different directory and plot the curvature
dw = 2.0
dh = [2.0, 4.0, 6.0, 8.0]
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
p2 = plot(xlabel = "Time (hrs)" ,
        ylabel = "Max. Stress (MPa)",
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

for h in dh
    workdir = "/home/srinath/data/Projects/moose_tests/supercomp_latest/curved_surface/glued/dw_" * string(2.0) *"_dh_" * string(h) * "/csv"
    cd(workdir)
    filebase = "curvature_time.csv"
    df = DataFrame(CSV.File(filebase))
    df2 = DataFrame(CSV.File("max_stress.csv"))
    labstr = string(h)
    plot!(p,df[!,1]/3600, -df[!,2], label=labstr)
    plot!(p2,df2[!,:Time]/3600, df2[!,"avg(Maximum)"], label=labstr)
end
cd("../../")

png(p, "curvature.png")
png(p2, "stress.png")

## Plot both values for dh 6.0
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
p2 = plot(xlabel = "Time (hrs)" ,
        ylabel = "Max. Stress (MPa)",
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
workdir = "/home/srinath/data/Projects/moose_tests/supercomp_latest/curved_surface/glued/dw_2.0_dh_4.0/csv"
cd(workdir)
filebase = "curvature_time.csv"
df = DataFrame(CSV.File(filebase))
df2 = DataFrame(CSV.File("max_stress.csv"))
plot!(p,df[!,1]/3600, -df[!,2], label="0.1")
plot!(p2,df2[!,:Time]/3600, df2[!,"avg(Maximum)"], label="0.1")
workdir = "/home/srinath/data/Projects/moose_tests/supercomp_latest/curved_surface/glued/increasing_current/dw_2.0_dh_4.0/csv"
cd(workdir)
filebase = "curvature_time.csv"
df = DataFrame(CSV.File(filebase))
df2 = DataFrame(CSV.File("max_stress.csv"))
plot!(p,df[!,1]/3600, -df[!,2], label="0.01")
plot!(p2,df2[!,:Time]/3600, df2[!,"avg(Maximum)"], label="0.01")
workdir = "/home/srinath/data/Projects/moose_tests/supercomp_latest/curved_surface/glued"
cd(workdir)
png(p, "curvature_dh_4.0_compare.png")
png(p2, "stress_dh_4.0_compare.png")
