using DataFrames, DataFramesMeta, CSV, Statistics
# using Colors, ColorSchemes, ColorSchemeTools
using Plots, Plots.PlotMeasures
using PyCall, Conda, LaTeXStrings
using LsqFit
# workdir = "C:\\Users\\s1.chakravar\\Documents\\Projects\\Li_swelling_moose\\stripping2"
workdir = "/home/srinath/repo/Projects/moose_tests/supercomp/stripping/new"
datadir = workdir
cd(datadir)
dirs = ["0.2mA","0.4mA","0.8mA", "1.2mA"]
labels = [0.2, 0.4, 0.8, 1.2]
# dirs = ["0.4mA"]
pyplot()
p = plot(xlabel = "time (a.u)",
        ylabel = L"$\frac{E_{WE}-E_0}{E_0}$",
        xmirror = false,
        framestyle = :box,
        legend = :bottomright,
        legendfontsize = 14,
        legendtitlefontsize = 16,
        tickfontsize = 14,
        guidefontsize = 16,
        legendtitle = L"$i (\frac{mA}{cm^2})$",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        title = "",
        titlefontsize = 16,
        grid = false,
        right_margin = 20mm,
        xlim = (0,500), ylim = (0,5))

for (i,d) in enumerate(dirs)
    cd(d)
    # if ( i == 1)
    #     cd("redo")
    # end
    filename = "Rct_low_" * d * "_csv.csv"
    println(filename)
    df = DataFrame(CSV.File(filename))
    # if (i == 1)
    df11 = @where(df, :time .== 70.0)
    # else
    #     df11 = @where(df, :time .== 150.025)
    # end
    E0 = df11[1,:over_potential]
    df2 = @where(df,:time .> 70)
    cd("../")
    # if (i == 1)
    #     cd("../")
    # end
    labstr = string(labels[i])
    y = (df2[!,:over_potential] .- E0) ./ E0
    t = df2[!,:time]
    println(t, y)
    # println(E0)
    if (i == 1)
        plot!(t ,y, label=labstr, ls=:dash)
    elseif ( i == 2)
        plot!(t,y, label=labstr, ls=:dash)
    elseif (i == 3)
        plot!(t,y, label=labstr, ls=:dash)
    else
        plot!(t,y, label=labstr, ls=:dash)
    end
    if ( i == 3)
        plot!(twinx(), t, -df2[!,:ext_pressure]*1e6,
                color=:black, guidefontsize=16, tickfontsize = 14,
                label="", ylabel=L"$\sigma (MPa)$", xlim = (0,500), ylim = (0.5, 3.2))
    end
end
# annotate!()
# p.o[:legend](bbox_to_anchor = (1.05,1), loc=2, borderaxespad = 1.0)
cd(workdir)
png(p,"test_over_potential.png")
