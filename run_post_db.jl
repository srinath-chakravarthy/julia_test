using JSON
using DataFrames, Plots, Statistics, DataFramesMeta
using Colors, ColorSchemes, ColorSchemeTools
using LsqFit, LaTeXStrings
using CSV, PyCall
sf = pyimport("scipy.ndimage.filters")
##
function replacenan(x,y)
    if length(x) != length(y)
        error("Array lengths do not match")
    end
    for i in eachindex(x)
        @inbounds x[i] = ifelse(isnan(x[i]), y[i], x[i])
    end
end


function smooth_util(x,y, ret::Bool=false)
    @. model(x,p) = 1.0/(1.0 + exp(-p[1]*(x-p[2])))
    p0 = [1.0, 0.25]
    fit = curve_fit(model, x, y, p0)
    newyy = model(x, coef(fit))
    ymax = [0.95, 0.98]
    @. intersect_y(ymax, coeffs) = coeffs[2] - (1.0/coeffs[1]) .* log.((1.0-ymax) / (ymax))
    xinter = nothing
    if coef(fit)[1] > 0.0
        xinter = intersect_y(ymax, coef(fit))
    end
    if !ret
        return newyy
    else
        return Dict("fit"=>fit, "model"=>model, "x"=>collect(0.25:0.05:3.0), "inter"=>xinter)
    end
end

function smooth_pack(x,y, ret::Bool = false)
    @. model(x,p) = 0.86/(1.0 + exp(-p[1]*(x-p[2])))
    p0 = [1.0, 0.25]
    fit = curve_fit(model, x, y, p0)
    newyy = model(x, coef(fit))
    if !ret
        return newyy
    else
        return Dict("fit"=>fit, "model"=>model, "x"=>collect(0.25:0.05:3.0))
    end
end
##
wt_ratio = collect(0.6:0.025:0.925)
lambda = [0.5, 0.75, 1.0, 1.25, 1.5, 1.667, 2.0, 2.5, 4.0, 5.0, 8.0]
workdir = "/home/srinath/repo/Projects/cathode_packing_data/new_ceder_data_rerun/"
stat = collect(1:5)
cd(workdir)
df = DataFrame([Float64, Float64, Float64, Float64, Float64, Float64],
    [:wt_ratio,:particle_size_ratio,:packing_density,:num_utilization, :utilization, :stat])
for w in wt_ratio
    wdir = "wt_ratio_" * string(w)
    if isdir(wdir)
        cd(wdir)
    else
        continue
    end
    cd("a_5.0")
    for l in lambda
        ldir = "ratio_" * string(l)
        pwd()
        if isdir(ldir)
            cd(ldir)
        else
            continue
        end
        ss::Int64 = 0
        for s in stat
            stat_dir = "stat_" * string(s)
            if isdir(stat_dir)
                cd(stat_dir)
            else
                continue
            end
            if isfile("results.json")
                res = JSON.parsefile("results.json")
                num_util = res["Results"]["Particle_fraction"]
                wt_util = res["Results"]["Active_weight_fraction"]
                pack = res["Results"]["Output_porosity"]
                ss += 1
                push!(df,[w,l,pack,num_util,wt_util, ss])
            end
            cd("../")
        end
        cd("../")
    end
    cd(workdir)
end
##
df_avg = by(df, [:wt_ratio, :particle_size_ratio],
num_utilization = :num_utilization=>mean, utilization = :utilization => mean, packing_density = :packing_density=>mean)
df_smooth= by(df_avg,[:wt_ratio],
    (:particle_size_ratio, :packing_density, :utilization, :num_utilization)
    => x->(particle_size_ratio = x.particle_size_ratio,
           packing_density = x.packing_density,
           utilization = x.utilization,
           pack2 = smooth_pack(x.particle_size_ratio, x.packing_density),
           util2 = smooth_util(x.particle_size_ratio, x.utilization),
           num_util2 = smooth_util(x.particle_size_ratio, x.num_utilization)))
## Now interpolate
wt_ratio1 = collect(0.65:0.005:0.925)
lambda1 = collect(0.5:0.025:5.0)
df_interp = DataFrame([Float64, Float64, Float64, Float64, Float64],
    [:wt_ratio,:particle_size_ratio,:packing_density,:utilization, :num_utilization])
for i in wt_ratio1
    for j in lambda1
        push!(df_interp, [i, j, 0.0, 0.0, 0.0])
    end
end
##
si = pyimport("scipy.interpolate")
points=convert(Matrix,df_smooth[:,[:wt_ratio, :particle_size_ratio]])
int_points=convert(Matrix,df_interp[:,[:wt_ratio, :particle_size_ratio]])
util = df_smooth[!,:util2]
num_util = df_smooth[!,:num_util2]
pack = df_smooth[!,:pack2]
int_util = si.griddata(points, util, int_points, method="linear")
int_num_util = si.griddata(points, num_util, int_points, method="linear")
int_pack = si.griddata(points, pack, int_points, method="linear")
df_interp[!,:utilization] = int_util
df_interp[!,:packing_density] = int_pack
df_interp[!,:num_utilization] = int_num_util
df_smooth_interp=  by(df_interp,[:wt_ratio],
    (:particle_size_ratio, :packing_density, :utilization, :num_utilization)
    => x->(particle_size_ratio = x.particle_size_ratio,
           packing_density = x.packing_density,
           utilization = x.utilization,
           pack2 = smooth_pack(x.particle_size_ratio, x.packing_density),
           util2 = smooth_util(x.particle_size_ratio, x.utilization),
           num_util2 = smooth_util(x.particle_size_ratio, x.num_utilization)))

##
#Contour Plots
cdict3 = Dict(:red =>  ((0.0, 1.0, 1.0),
                   (0.2, 1.0, 1.0),
                   (0.4, 1.0, 1.0),
                   (0.5, 1.0, 1.0),
                   (0.9, 0.9, 0.9),
                   (1.0, 0.0, 0.0)),

         :green => ((0.0, 0.0, 0.0),
                   (0.4, 0.0, 0.0),
                   (0.95, 0.9, 0.9),
                   (1.0, 1.0, 1.0)),

         :blue =>  ((0.0, 1.0, 1.0),
                   (0.2, 1.0, 1.0),
                   (0.4, 1.0, 1.0),
                   (0.5, 0.0, 0.0),
                   (1.0, 0.0, 0.0)))
scheme = make_colorscheme(cdict3)
## Ceder paper plot 2a
df2 = @where(df_smooth_interp, :particle_size_ratio .== 2.0)
pyplot()
# ylabel = L"$\lambda = \frac{D_{CAM}^{SM}}{D_{SE}}$"
xlabel = L"$f_{CAM} (wt \%)$"
ylabel = L"$\theta_{CAM} (\%)$"
# ylabel2 = L"$\rho (\%)$"
p = plot(xlabel = xlabel,
        ylabel = ylabel,
        xmirror = false,
        framestyle = :box,
        legend = :outertopright,
        legendfontsize = 14,
        legendtitlefontsize = 18,
        tickfontsize = 16,
        guidefontsize = 18,
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        grid = false,
        xlim = (65,93))
plot!(plot(df2[!,:wt_ratio]*100.0, sf.gaussian_filter1d(df2[!,:pack2], sigma=2.0)*100.0))
cd(workdir)
# png(p, "ceder_paper_plot2a.png")
