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
wt_ratio = collect(0.6:0.025:0.95)
lambda = [0.5, 0.75, 1.0, 1.25, 1.5, 1.667, 2.0, 2.5, 4.0, 5.0, 6.0, 8.0, 10.0]
workdir = "/home/srinath/Projects/cathode_packing_data/new_ceder_data_rerun/"
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
                if pack > 0.4
                    ss += 1
                    push!(df,[w,l,pack,num_util,wt_util, ss])
                else
                    cd("../")
                    continue
                end
            end
            cd("../")
        end
        cd("../")
    end
    cd(workdir)
end
##
df_avg = by(df, [:wt_ratio, :particle_size_ratio],
                num_utilization = :num_utilization=>mean,
                utilization = :utilization => mean,
                packing_density = :packing_density=>mean,
                stat = :stat=>maximum)
df_smooth= by(df_avg,[:wt_ratio],
    (:particle_size_ratio, :packing_density, :utilization, :num_utilization)
    => x->(particle_size_ratio = x.particle_size_ratio,
           packing_density = x.packing_density,
           utilization = x.utilization,
           pack2 = smooth_pack(x.particle_size_ratio, x.packing_density),
           util2 = smooth_util(x.particle_size_ratio, x.utilization),
           num_util2 = smooth_util(x.particle_size_ratio, x.num_utilization)))
cd(workdir)
CSV.write("orig_data_avg.csv", df_avg)
CSV.write("orig_data_avg_smooth.csv", df_smooth)
## Now interpolate
wt_ratio1 = collect(0.65:0.005:0.95)
lambda1 = collect(0.5:0.025:9.0)
df_interp = DataFrame([Float64, Float64, Float64, Float64, Float64],
    [:wt_ratio,:particle_size_ratio,:packing_density,:utilization, :num_utilization])
for i in wt_ratio1
    for j in lambda1
        push!(df_interp, [i, j, 0.0, 0.0, 0.0])
    end
end
##
si = pyimport("scipy.interpolate")
sf = pyimport("scipy.ndimage.filters")
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
CSV.write("data_smooth_interpolated.csv", df_smooth_interp)
##
#Contour Plots
pyplot()

cdict3 = Dict(:red =>  ((0.0, 1.0, 1.0),
                   (0.2, 1.0, 1.0),
                   (0.4, 1.0, 1.0),
                   (0.5, 1.0, 1.0),
                   (0.9, 1.0, 1.0),
                   (1.0, 0.0, 0.0)),

         :green => ((0.0, 0.0, 0.0),
                   (0.4, 0.0, 0.0),
                   (0.95, 1.0, 1.0),
                   (1.0, 1.0, 1.0)),

         :blue =>  ((0.0, 1.0, 1.0),
                   (0.2, 1.0, 1.0),
                   (0.4, 1.0, 1.0),
                   (0.5, 0.0, 0.0),
                   (1.0, 0.0, 0.0)))
xlabel = L"$\lambda = \frac{D_{CAM}^{SM}}{D_{SE}}$"
ylabel = L"$f_{CAM} (wt \%)$"
ylabel2 = L"$\rho (\%)$"
p = plot(xlabel = xlabel,
       ylabel = ylabel,
       xmirror = false,
       framestyle = :box,
       legend = :outertopright,
       legendfontsize = 14,
       legendtitlefontsize = 18,
       tickfontsize = 16,
       guidefontsize = 18,
       legendtitle = L"$\theta_{CAM}$",
       foreground_color_legend=nothing,
       background_color_legend=nothing,
       grid = false,
       ylim = (65,90))

scheme = make_colorscheme(cdict3)
x = unique(df_smooth_interp[!,:particle_size_ratio])
y = unique(df_smooth_interp[!,:wt_ratio])*100.0
y1 = 95.0 .* y ./(y .+ (99 .-y))
y1 = y./99 * 100.0
# y1 = y ./1.05
z = df_smooth_interp[!,:util2]
z2 = reshape(z, length(x), length(y))
z1 = sf.gaussian_filter(reshape(z, length(x), length(y)), sigma=5.0, mode=["nearest", "nearest"])
plot!(x,y1, z1'*100.0, st=:contourf, seriescolor=cgrad(scheme.colors), levels=[20,30, 40,50, 60,70, 80,90, 100])
plot!(x,y1, z1'*100.0, st=:contourf, seriescolor=cgrad(scheme.colors), levels=256)
plot!(x,y1, z1'*100.0, st=:contour, seriescolor=:black, levels=[95.0,98.0])
cd(workdir)
# plot!(x,y1, z1', st=:contourf, levels=256)
png(p,"test_contour_psi_0.285_phi_0.25_new_fcam.png")

## Ceder paper plot 2a
df2 = @where(df_smooth_interp, :particle_size_ratio .== 1.675)
# pyplot()
gr()
# ylabel = L"$\lambda = \frac{D_{CAM}^{SM}}{D_{SE}}$"
xlabel = L"$f_{CAM} (\mathrm{wt} \%)$"
ylabel = L"$\theta_{\mathrm{CAM}} (\%)$"
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
        xlim = (65,92))
x = df2[!,:wt_ratio]*100.0
x1 = 95.0.* x ./(x .+ (99 .-x))
# x1 = x./1.05
plot!(x1, sf.gaussian_filter1d(df2[!,:num_util2],sigma=3.0),label="", lw=4)
png(p, "ceder_paper_plot2a.png")
## Flatten z1 and put values into df_smooth_interp
function flatten_new(Mat)
    nrows,mcols = size(Mat)
    flattened = zeros(nrows*mcols)
    kk = 1
    for i = 1:mcols
        for j = 1:nrows
            flattened[kk] = Mat[j,i]
            kk += 1
        end
    end
    return flattened
end
## RBF trials
# rbfi = si.Rbf(df_avg[!,:wt_ratio], df_avg[!,:particle_size_ratio], df_avg[!,:utilization],"multiquadric", epsilon = 1e-4, smooth=5.0)
# # XI, YI = np.meshgrid(df_interp[!,:wt_ratio], df_interp[!,:particle_size_ratio])
# z = rbfi(df_smooth_interp[!,:wt_ratio], df_smooth_interp[!,:particle_size_ratio])
# # ZI = rbfi(XI,YI)
# # df_smooth_interp[!,:util3] = rbfi(df_interp[!,:wt_ratio], df_interp[!,:particle_size_ratio])
# p = plot(xlabel = xlabel,
#        ylabel = ylabel,
#        xmirror = false,
#        framestyle = :box,
#        legend = :outertopright,
#        legendfontsize = 14,
#        legendtitlefontsize = 18,
#        tickfontsize = 16,
#        guidefontsize = 18,
#        legendtitle = L"$\theta_{CAM}$",
#        foreground_color_legend=nothing,
#        background_color_legend=nothing,
#        grid = false,
#        ylim = (65,87))
#
# scheme = make_colorscheme(cdict3)
# x = unique(df_smooth_interp[!,:particle_size_ratio])
# y = unique(df_smooth_interp[!,:wt_ratio])*100.0
# # y1 = 95.0 .* y ./(y .+ (99 .-y))
# # y1 = y ./1.05
# # z = df_smooth_interp[!,:util3]*100.0
# z1 = sf.gaussian_filter(reshape(z, length(x), length(y)), sigma=0.1)*100.0
#
# # plot!(x,y1, z1', st=:contourf, seriescolor=cgrad(scheme.colors), levels=[10, 20,30, 40,50, 60,70, 80,90, 100])
# # plot!(x,y1, z1', st=:contourf, seriescolor=cgrad(scheme.colors), levels=256)
# # plot!(x,y1, z1', st=:contour, seriescolor=:black, levels=[95.0,98.0])
# cd(workdir)
# plot!(x,y, z1', st=:contourf,seriescolor=cgrad(scheme.colors), levels=256)
# png(p,"test_contour_psi_0.285_phi_0.25_new.png")
## Ceder f_cam critical lambda
df_smooth_interp[!,:fcam] .= df_smooth_interp[!,:wt_ratio] * 0.95/0.99
gd_wr = groupby(df_smooth_interp, :wt_ratio)
l_crit = Float64[]
f_cam_ceder = Float64[]
for (i,wrr) in enumerate(keys(gd_wr))
    global l_crit, f_cam_ceder
    dfx = gd_wr[wrr]
    aa = @where(dfx, :util2 .> 0.998)
    if nrow(aa) > 0
        # l_crit[i] = aa[1,:particle_size_ratio]
        # f_cam_ceder[i] aa[1,:fcam]
        append!(l_crit, 15/aa[1,:particle_size_ratio])
        append!(f_cam_ceder, aa[1,:wt_ratio])
    end
    # append!(f_cam_ceder, aa[1,:fcam])
end
## Simple plot of utilization vs particle size ratio for different f_cam_ceder
xlabel = L"$\lambda = \frac{D_{CAM}}{D_{SE}}$"
ylabel = L"$\theta_{\mathrm{CAM}} (\%)$"
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
        legendtitle = L"f_{CAM}")

wrr = [0.825, 0.85, 0.875, 0.9, 0.925]

colors = [:black,:blue,:red, :orange, :green]
for (i,wr) in enumerate(wrr)
    dfx = @where(df_smooth_interp, :wt_ratio .== wr)
    labstr = string(wr*100.0)
    plot!(15.0/dfx[!,:particle_size_ratio], dfx[!,:util2]*100.0, color = colors[i], label = labstr)
end
png(p, "simple_plot.png")
