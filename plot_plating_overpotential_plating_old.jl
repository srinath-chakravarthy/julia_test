using DataFrames, DataFramesMeta, CSV, Statistics
# using Colors, ColorSchemes, ColorSchemeTools
using Plots, Plots.PlotMeasures
using PyCall, Conda, LaTeXStrings
using LsqFit
# workdir = "C:\\Users\\s1.chakravar\\Documents\\Projects\\Li_swelling_moose\\stripping2"
# workdir = "/home/srinath/Projects/moose_tests/supercomp_latest/plating"
workdir = "/media/srinath/e5447af9-9fc4-4429-a0ef-31a12ec2f776/home/srinath/ubuntu_20_backup/moose_tests/stripping/"
datadir = workdir
cd(datadir)
dirs = ["0.2mA","0.4mA","0.8mA"]
currents = [0.2, 0.4, 0.8] # ma/cm^2
Faraday = 96485.3329 # C/mol
rho = 1e-5 # m^3/mol
# dirs = ["0.4mA"]
pyplot()
p = plot(xlabel = L"        Plated Li $(\mu m)$",
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
        ylim = (0,5))

for (i,d) in enumerate(dirs)
    cd(d)
    if ( i == 1)
        cd("redo")
    end
    filename = "Rct_low_" * d * "_csv.csv"
    println(filename)
    df = DataFrame(CSV.File(filename))
    # if (i == 1)
    df11 = @where(df, :time .== 150.0)

    # else
    #     df11 = @where(df, :time .== 150.025)
    # end
    # E0 = df11[1,:over_potential]
    df2 = @where(df,:time .> 150)
    E0 = df2[1,:over_potential]
    cd("../")
    if (i == 1)
        cd("../")
    end
    labstr = string(currents[i])
    y = (df2[!,:over_potential] .- E0) ./ E0 * 5.0
    t = df2[!,:time]
    # if (i == 1)
    #     t = t
    # elseif ( i == 2)
    #     t = t./1.1
    # elseif (i == 3)
    #     t = t./1.3
    # else
    #     t = t./1.35
    # end

    plated_Li = t * 10.0 * 10.0 * currents[i]/ Faraday * rho * 1e6 # um
    # println(t, y)
    println(E0)
    plot!(plated_Li ,y, label=labstr, ls=:dash)

    # if ( i == 1)
        plot!(twinx(), plated_Li, -df2[!,:ext_pressure]*1e6,
                color=:black, guidefontsize=16, tickfontsize = 14,
                label="", ylabel=L"$\sigma (MPa)$", ylim = (0.0, 3.2))
    # end
    println(last(df2,10))
end
# annotate!()
# p.o[:legend](bbox_to_anchor = (1.05,1), loc=2, borderaxespad = 1.0)
cd(workdir)
png(p,"test_over_potential_plated_li.png")
