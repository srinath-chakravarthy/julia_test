# Find curvature of top surface
using CSV, DataFrames, Statistics, Glob
using Dierckx
using Plots
using LaTeXStrings
using Printf
## ---- Import equilibrium potential file first ----
basedir="/home/srinath/github/electro_chemo_mech2/problems/cathode_trial"
cd(basedir)
equil_file="nmc_equilibrium_potential.csv"
df=DataFrame(CSV.File(equil_file))
csv_dir=basedir * "/csv"
cd(csv_dir)
rct = [1.0 10.0 100.0 200.0 500.0 1000.0]
current = [3e-3 5e-3 10e-3 15e-3 20e-3]
current2 = similar(current)
rct2 = similar(rct)
@. rct2 = 1.0 / rct * 1e-2
@. current2 = current * 1000
xlabel = L"Capacity $(mAh/cm^2)$"
ylabel = L"Voltage $(V)$"
pyplot()
p = plot(xlabel = xlabel,
        ylabel = ylabel,
        xmirror = false,
        framestyle = :box,
        legend = :outertopright,
        legendfontsize = 14,
        legendtitlefontsize = 18,
        tickfontsize = 16,
        guidefontsize = 18,
        legendtitle = L"i $(mA/cm^2)$",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        grid = false,
        ylim = (2.7, 4.3),
        xlim = (0.01, 35))
for c in current2
    r = rct[1]
    st = @sprintf "%.1f" c
    filename = "rct_" * string(Int(r)) * "_i_" * st *".csv"
    xx = string(c/10.0)
    labstr = xx
    dfx = DataFrame(CSV.File(filename))
    plot!(dfx[!,"time"].*c ./3600.0, dfx[!,"Voltage_Cathode"] ./ 1000, label=labstr)
    println(filename)
end

png(p, "test.png")
