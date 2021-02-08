using DataFrames, DataFramesMeta, CSV, Statistics
# using Colors, ColorSchemes, ColorSchemeTools
using Plots, Plots.PlotMeasures
using PyCall, Conda, LaTeXStrings
using LsqFit
# workdir = "C:\\Users\\s1.chakravar\\Documents\\Projects\\Li_swelling_moose\\stripping2"
workdir = "/home/srinath/Projects/moose_tests/supercomp_latest/stripping/0.4mA"
datadir = workdir
cd(datadir)
dirs = ["2.0","4.0","./", "10.0"]
currents = [2.0,4.0,6.0,10.0] # ma/cm^2
Faraday = 96485.3329 # C/mol
rho = 1e-5 # m^3/mol
# dirs = ["0.4mA"]
pyplot()
p = plot(xlabel = L" time (a.u)",
        ylabel = L"$\frac{E_{WE}-E_0}{E_0}$",
        xmirror = false,
        framestyle = :box,
        legend = :topright,
        legendfontsize = 14,
        legendtitlefontsize = 16,
        tickfontsize = 14,
        guidefontsize = 16,
        legendtitle = L"$ \Delta_{def} (\mu m)$",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        title = "",
        titlefontsize = 16,
        grid = false,
        right_margin = 20mm,
        ylim = (0,5))

for (i,d) in enumerate(dirs)
    if (i != 3)
        cd("size_effects/" * d)
    end
    filename = "Rct_low_0.4mA_csv.csv"
    println(filename)
    df = DataFrame(CSV.File(filename))
    df11 = @where(df, :time .== 90.0)
    df2 = @where(df,:time .> 90)
    E0 = df2[1,:over_potential]
    if (i != 3)
        cd("../../")
    end

    labstr = string(currents[i])
    y = (df2[!,:over_potential] .- E0) ./ E0 * 5
    t = df2[!,:time]
    plated_Li = t * 10.0 * 10.0 * currents[i]/ Faraday * rho * 1e6 # um
    # println(t, y)
    println(E0)
    plot!(t ,y, label=labstr, ls=:dash)
    # if (i == 1)
    #     plot!(t ,y, label=labstr, ls=:dash)
    # elseif ( i == 2)
    #     plot!(plated_Li,y, label=labstr, ls=:dash)
    # elseif (i == 3)
    #     plot!(plated_Li,y, label=labstr, ls=:dash)
    # else
    #     plot!(plated_Li,y, label=labstr, ls=:dash)
    # end
    if ( i == 3)
        plot!(twinx(), t, -df2[!,:ext_pressure]*1e6,
                color=:black, guidefontsize=16, tickfontsize = 14,
                label="", ylabel=L"$\sigma (MPa)$", ylim = (0.0, 3.2))
    end
    # println(last(df2,10))
end
# annotate!()
# p.o[:legend](bbox_to_anchor = (1.05,1), loc=2, borderaxespad = 1.0)
cd(workdir)
png(p,"test_over_potential.png")
