using DataFrames, DataFramesMeta, CSV, Statistics
using Colors, ColorSchemes, ColorSchemeTools
using JSON
using Plots
using PyCall, Conda, LaTeXStrings
using LsqFit
import GR.meshgrid
using Contour

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

function run_post(filename)
    try
        # println(pwd())
        data = JSON.parsefile(filename; dicttype=Dict, inttype=Int64, use_mmap=true)
        if !haskey(data, "Results")
            error("Results file not found")
        else
            wt = 0.0
            wt1 = 0.0
            wt2 = 0.0
            wt2_py = 0.0
            a = data["Results"]
            wt = a["Output_porosity"]
            wt1 = a["Particle_fraction"]
            wt2 = a["Active_weight_fraction"]
            if isfile("cathode.json")
                data2 = JSON.parsefile("cathode.json"; dicttype=Dict, inttype=Int64, use_mmap=true)
                if haskey(data2, "Results")
                    for b in data2["Results"]
                        wt2_py = parse(Float64,b["Active_Weight_fraction"])
                    end
                end
            end
            return wt, wt1, wt2, wt2_py
        end
    catch
        error("File Not found")
    end

end

function replacenan(x,y)
    if length(x) != length(y)
        error("Array lengths do not match")
    end
    for i in eachindex(x)
        @inbounds x[i] = ifelse(isnan(x[i]), y[i], x[i])
    end
end


function smooth2(x, y; ret::Bool = false, ymax::Float64 = 0.9)
    @. model(x,p) = p[1]/(p[2] + exp(-p[3]*(x-p[4])))
    p0 = [0.5, 0.5, 1.0, 0.25]
    fit = curve_fit(model, x, y, p0)
    newyy = model(x, coef(fit))
    println(ymax)
    xinter = 100.0
    if !ret
        @. intersect_y(ymax, coeffs) = coeffs[4] - (1.0/coeffs[3]) .* log.((coeffs[1]-coeffs[2]*ymax) / (ymax))
        if coef(fit)[1] > 0.0
            xinter = intersect_y(ymax, coef(fit))
        end
    end
    if !ret
        return xinter
    else
        return newyy
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
        return Dict("fit"=>fit, "model"=>model, "x"=>collect(0.25:0.05:3.0), "inter"=>xinter)
    else
        return newyy
    end
end

function smooth_pack(x,y, ret::Bool=false)
    @. model(x,p) = 0.86/(1.0 + exp(-p[1]*(x-p[2])))
    p0 = [1.0, 0.25]
    fit = curve_fit(model, x, y, p0)
    newyy = model(x, coef(fit))
    if !ret
        return Dict("fit"=>fit, "model"=>model, "x"=>collect(0.25:0.05:3.0))
    else
        return newyy
    end
end


function replacenan(x,y)
    if length(x) != length(y)
        error("Arrays have to be of equal length")
    end
    for i in eachindex(x)
        x[i] = ifelse(isnan(x[i]), y[i], x[i])
    end
end

function smooth(x,y)
    si = pyimport("scipy.interpolate")
    ss = pyimport("scipy.signal")
    np = pyimport("numpy")
    z = []
    try
        itp = si.interp1d(x,y,kind="linear")
        window_size, poly_order = 171, 3
        xx = np.concatenate((np.arange(0.25,2.0, 0.005), np.array([2.0])))
        ysmooth = ss.savgol_filter(itp(xx), window_size, poly_order)
        itp2 = si.interp1d(xx,ysmooth,kind="linear")
        z = itp2(x)
    catch
        z = y
    end
    return z
end
## Get data
work_dir = "/home/srinath/Projects/cathode_packing_data/bimodal2"
plot_dir = "/home/srinath/Projects/cathode_packing_data/bimodal2/plots"
bimodal_dir="/home/srinath/Projects/cathode_packing_data/bimodal2/bimodal_data/"
unimodal_dir= "/home/srinath/Projects/cathode_packing_data/bimodal2/unimodal_data/"
cd(bimodal_dir)
df_bimodal = DataFrame(CSV.File("all_data_smooth_interp.csv"))
cd(unimodal_dir)
df_unimodal = DataFrame(CSV.File("all_data_smooth_interp.csv"))
cd(work_dir)

## Smooth Unimodal data contour
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
x = unique(df_unimodal[!,:particle_size_ratio])
y = unique(df_unimodal[!,:wt_ratio])*100.0
x1 = x
y1 = y
# y1 = 95.0 .* y ./(y .+ (99 .-y))
# y1 = y./99 * 100.0
# y1 = y ./1.05
z = df_unimodal[!,:util2]
z2 = reshape(z, length(x), length(y))
z1 = sf.gaussian_filter(reshape(z, length(x), length(y)), sigma=5.0, mode=["nearest", "nearest"])
z3 = flatten_new(z1)
df_unimodal[!,:util_smooth] = z3
plot!(x,y1, z1'*100.0, st=:contourf, seriescolor=cgrad(scheme.colors), levels=[20,30, 40,50, 60,70, 80,90, 100])
plot!(x,y1, z1'*100.0, st=:contourf, seriescolor=cgrad(scheme.colors), levels=256)
plot!(x,y1, z1'*100.0, st=:contour, seriescolor=:black, levels=[95.0,98.0])
cd(plot_dir)
# plot!(x,y1, z1', st=:contourf, levels=256)
png(p,"test_contour_unimodeal.png")
cd(work_dir)
##Plotting
pyplot()

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


## Compare unimodal to bimodal
xlabel = L"$\lambda = \frac{D_{CAM}^{SM}}{D_{SE}}$"
# ylabel = L"$\theta_{CAM} (wt \%)$"
ylabel = L"Specific Capacity $(\mathrm{mAh \;{g}^{-1}})$"
ylabel2 = L"$\rho (\%)$"
xlabel_real = L"$D_{SE} (\mu m)$"
ylabel_real = ylabel
title = L"$\psi =  0.285, \phi = 0.25$"
Dcambig = 16.2
brr = 0.3
bmr = 0.25
Dcamsmall = Dcambig * brr *1.06
dcams = round(Dcamsmall,digits=1)
colors = [:blue, :black, :red]
title_real = latexstring("\$D_{CAM}^{BIG} = $(Dcambig) \\mu m, D_{CAM}^{SM} = $(dcams) \\mu m, \\phi = $(bmr)\$")

# title_real = L"D_{CAM}^{BIG} = 14 \mu m, D_{CAM}^{SM} = 4.0 \mu m, \phi = 0.25"
p = plot(xlabel = xlabel_real,
        ylabel = ylabel_real,
        xmirror = false,
        framestyle = :box,
        legend = :outertopright,
        legendfontsize = 14,
        legendtitlefontsize = 16,
        tickfontsize = 14,
        guidefontsize = 16,
        legendtitle = L"$f_{CAM} (\mathrm{wt} \%)$",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        title = title_real,
        titlefontsize = 16,
        grid = false,
        ylim = (100,193.5),
        xlim = (0.95,6))
expt_data = DataFrame()
expt_data.dse = [1.0, 2.0, 3.0]
expt_data.cap1 = [185, 192, 190]
expt_data.cap2 = [169, 177, 176]
expt_data.cap3 = [153, 162, 163]
expt_data.util1 = expt_data.cap1/expt_data[2,2]
plot!(expt_data[!,:dse], expt_data[!,:cap1], st=:scatter, color=:red, marker=:star, label="", markersize=15)

# Dcambig = 15.0
qwrr = [0.85, 0.9, 0.925]
for (i,wrr) in enumerate(qwrr)
    df_uni = @where(df_unimodal, :wt_ratio .== wrr)
    df_bi = @where(df_bimodal, :wt_ratio .== wrr, :bimodal_radius_ratio .== brr,
                    :bimodal_mix_ratio .== bmr)
    labstr = string(wrr*100.0)
    # plot!(Dcambig./df_uni[!,:particle_size_ratio], df_uni[!,:util_smooth]*100.0, ls=:dash, color = colors[i], label=labstr)
    # plot!(Dcamsmall./df_bi[!,:particle_size_ratio], df_bi[!,:util2]*100.0, ls=:solid, color = colors[i], label=labstr)
    x = df_bi[!,:particle_size_ratio]
    y = df_bi[!,:util2]
    dict = smooth_util(x,y, false)
    model = dict["model"]
    fit = dict["fit"]
    newxx = collect(0.5:0.05:6.0)
    newyy = model(newxx, coef(fit))
    # plot!(Dcamsmall./df_bi[!,:particle_size_ratio], df_bi[!,:util2]*100.0, st=:scatter, color = colors[i], label=labstr, marker=:x)
    plot!(Dcamsmall./newxx, newyy*expt_data[2,2], color= colors[i], label = labstr)
end
# plot!(collect(0:2.0:100.1), [99.5], seriestype=:hline, ls=:dash, color = :black, label="")
cd(plot_dir)
png(p, "compare_uni_bi_2.png")
cd(work_dir)
## Contour Plots
using Contour
gr()
sf = pyimport("scipy.ndimage.filters")
xlabel_real = L"$f_{CAM} (\mathrm{wt} \%)$"
ylabel_real = L"$D_{SE} (\mu m)$"
xlabel = L"$\lambda = \frac{D_{CAM}^{BIG}}{D_{SE}}$"
Dcambig = 16.2
brr = 0.32
bmr = 0.25
Dcamsmall = Dcambig * brr
dcams = round(Dcamsmall,digits=1)
title_real = latexstring("\$D_{CAM}^{BIG} = $(Dcambig) \\mu m, D_{CAM}^{SM} = $(dcams) \\mu m, \\phi = $(bmr)\$")
p = plot(xlabel = xlabel,
        ylabel = xlabel_real,
        zlabel = L"$\theta_{CAM} (\mathrm{wt} \%)$",
        xmirror = false,
        frame_style = :box,
        legend = :none,
        legendfontsize = 14,
        legendtitlefontsize = 16,
        tickfontsize = 14,
        guidefontsize = 16,
        # legendtitle = "Distribution",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        contour_labels = true,
        title = "",
        titlefontsize = 16,
        grid = false)

x = unique(df_unimodal[!,:particle_size_ratio])
y = unique(df_unimodal[!,:wt_ratio])*100.0
z = df_unimodal[!,:util2]*100.0
z2 = reshape(z, length(x), length(y))
zsmooth = sf.gaussian_filter(z2,sigma=4.0)
contours_uni = contours(x,y,zsmooth,[90.0])
cl = first(Contour.levels(contours_uni))
line = first(Contour.lines(cl))

xs, ys = Contour.coordinates(line)
newys = smooth2(xs,ys/100.0, ret= true, ymax = 0.9) * 100.0
# plot!(x,y, zsmooth', st=:contour, lt= :dash, seriescolor=:black, levels=[98.99,99.0], xlim = (3,10))
# plot!(newys, Dcambig./xs, ls = :dash, label="Unimodal", color = :black)
xs_save = xs
ys_save = newys

# plot!(ys, Dcambig./xs, seriestype = :scatter, marker=:x)

bmr = 0.25
brr = 0.285
dcams = Dcambig * brr
df_bi = @where(df_bimodal, :bimodal_radius_ratio .== brr,
                :bimodal_mix_ratio .== bmr)
x = unique(df_bi[!,:particle_size_ratio])
y = unique(df_bi[!,:wt_ratio])*100.0
z = df_bi[!,:util2]*100.0
z2 = reshape(z, length(x), length(y))
zsmooth = sf.gaussian_filter(z2,sigma=1.0)
zsmooth2 = flatten_new(zsmooth)
contours_uni = contours(x,y,zsmooth,[90.0])
cl = first(Contour.levels(contours_uni))
line = first(Contour.lines(cl))

xs, ys = Contour.coordinates(line)
newys = smooth2(xs,ys/100.0, ret= true, ymax = 0.95) * 100.0

plot!(x/brr,y, zsmooth', st=:contour, levels=[50.0, 75.0, 90.0, 95.0, 99.0], label="", seriescolor = :black, clabels=true, margin=10mm)
# plot!(x/brr, y, zsmooth2, st=:wireframe, camera=(70,30), margin=10mm)
# plot!(x/brr,y, zsmooth', st=:contour, seriescolor=cgrad(scheme.colors), levels=256)
# plot!(x/brr,y, zsmooth', st=:contour, seriescolor=:black, levels=[98.99,99.0])
# plot!(xs_save, ys_save, color = :black, ls=:dash, xlim=(2.5,8), ylim=(70,92), label="")
# plot!(xs/brr, newys, ls = :solid, color = :black,xlim=(2.5,8), ylim=(70,92), label="")

cd(plot_dir)
png(p, "compare_uni_bi_contour_90util.png")
cd(work_dir)
p = plot(xlabel = ylabel_real,
        ylabel = xlabel_real,
        xmirror = false,
        framestyle = :box,
        legend = :outertopright,
        legendfontsize = 14,
        legendtitlefontsize = 16,
        tickfontsize = 14,
        guidefontsize = 16,
        legendtitle = "Distribution",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        title = title_real,
        titlefontsize = 16,
        grid = false,
        xlim = (1.6,5),
        ylim = (80,92))

plot!(dcams./xs, newys, label="Bimodal")
plot!(Dcambig./xs_save, ys_save, label="Unimodal")


cd(plot_dir)
png(p, "compare_uni_bi_90_util.png")
cd(work_dir)


## Variation with bimodal radius ratio_
using Contour
xlabel_real = L"$D_{CAM} (\mu m)$"
ylabel_real = L"$\phi$"
zlabel_real = L"$D_{SE} (\mu m)$"
p = plot(xlabel = xlabel_real,
        ylabel = ylabel_real,
        zlabel = zlabel_real,
        xmirror = false,
        framestyle = :box,
        legend = :outertopright,
        legendfontsize = 14,
        legendtitlefontsize = 16,
        tickfontsize = 14,
        guidefontsize = 16,
        legendtitle = "Distribution",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        title = title_real,
        titlefontsize = 16,
        grid = false)
df2 = @where(df_bimodal, :wt_ratio .== 0.925, :bimodal_mix_ratio .== 0.25, :bimodal_radius_ratio .> 0.24)
x = unique(df2[!,:particle_size_ratio])
y = unique(df2[!,:bimodal_radius_ratio])
z = df2[!,:util2]
z2 = reshape(z, length(x), length(y))
zsmooth = sf.gaussian_filter(z2,sigma=5.0)
plot!(x/brr,y, zsmooth', st=:contourf, seriescolor=cgrad(scheme.colors))
plot!(x/brr,y, zsmooth', st=:contourf, seriescolor=cgrad(scheme.colors), levels=256)
plot!(x/brr,y, zsmooth', st=:contour, seriescolor=:black, levels=[90.0,95.0])
# contours_uni = Contour.contours(x,y,zsmooth,levels=10)
# line = first(Contour.lines(cl))
# cl = first(Contour.levels(contours_uni))

# xs, ys = Contour.coordinates(line)

cd(plot_dir)
png(p, "test.png")
cd(work_dir)
## create db for fcam = 90% and lambda_critical = 99 bmr brr_dir, brr > 0.245 and brr <
using Plots.PlotMeasures
pyplot()
cdict4 = Dict(:red =>  ((0.0, 1.0, 1.0),
                   (0.2, 0.8, 0.8),
                   (0.4, 0.6, 0.6),
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

scheme2 = make_colorscheme(cdict4)

df2 = @where(df_bimodal, :wt_ratio .== 0.85,
            :bimodal_radius_ratio .> 0.245, :bimodal_radius_ratio .<0.405,
            :bimodal_mix_ratio .>0.195, :bimodal_mix_ratio .<0.505,
            :particle_size_ratio .> 0.75)
df_new = DataFrame([Float64, Float64, Float64, Float64, Float64],
    [:bimodal_radius_ratio,:bimodal_mix_ratio, :lambda_crit, :dcamsmall, :dse])
gd = groupby(df2, :bimodal_radius_ratio)
Dcambig = 16.2

xlabel_real = L"$D_{CAM}^{SM} (\mu m)$"
ylabel_real = L"$\phi$"
zlabel_real = L"$D_{SE} (\mu m)$"
zlabel = L"$\lambda$"
title = L"$f_{CAM} = 90 (\mathrm{wt} \%)$"
xlabel = L"$\psi$"
p = plot(xlabel = ylabel_real,
        ylabel = xlabel,
        zlabel = zlabel,
        xmirror = false,
        framestyle = :box,
        legend = :outertopright,
        legendfontsize = 14,
        legendtitlefontsize = 16,
        tickfontsize = 12,
        guidefontsize = 16,
        legendtitle = "Distribution",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        title = title,
        titlefontsize = 16,
        grid = false)

for (i,brr) in enumerate(keys(gd))
    df_brr = gd[brr]
    println(brr.bimodal_radius_ratio)
    gd_bmr = groupby(df_brr,:bimodal_mix_ratio)
    for (j, bmr) in enumerate(keys(gd_bmr))
        df_bmr = gd_bmr[bmr]
        xinter = smooth2(df_bmr[!,:particle_size_ratio], df_bmr[!,:util2], ret = false, ymax = 0.95)
        # dfxx = @where(df, :y .> 0.98)
        # xinter = 100.0
        # if nrow(dfxx) > 0
            # println(dfxx[1,:])
            # xinter = dfxx[1,1]
        # end
        println(xinter)
        # println(brr.bimodal_radius_ratio, bmr.bimodal_mix_ratio, xinter)
        push!(df_new, [brr.bimodal_radius_ratio, bmr.bimodal_mix_ratio, xinter,
            Dcambig*brr.bimodal_radius_ratio,
            Dcambig*brr.bimodal_radius_ratio/xinter])

        #
        # plot!(df_bmr[!,:particle_size_ratio], newy)
    end
end
sort!(df_new)

df_new2 = by(df_new, [:bimodal_mix_ratio],
            (:bimodal_radius_ratio, :lambda_crit, :dcamsmall, :dse)
            => x->(bimodal_radius_ratio = x.bimodal_radius_ratio,
                   lambda_crit = sf.gaussian_filter1d(x.lambda_crit, sigma = 1.0)./x.bimodal_radius_ratio,
                   dcamsmall = x.dcamsmall,
                   dse = sf.gaussian_filter1d(x.dse, sigma = 1.0)))
sort!(df_new2)
sf = pyimport("scipy.ndimage.filters")
x = unique(df_new2[!,:bimodal_radius_ratio])
y = unique(df_new2[!,:bimodal_mix_ratio])
z = df_new2[!,:lambda_crit]
z1 = reshape(z, length(x), length(y))
zsmooth = sf.gaussian_filter(z1, sigma=4.0)
df_new2[!,:l2] = flatten_new(zsmooth)
# plot!(y,x, zsmooth, st=:contourf, seriescolor=cgrad(scheme.colors), levels=10)
# plot!(y,x, zsmooth, st=:surface, levels=256, seriescolor=cgrad(scheme.colors), camera = (50,30), margin = 8mm)
# plot!(x,y, zsmooth', st=:contour, seriescolor=:black, levels=[2.8,2.6])
plot!(y,x, zsmooth, st=:wireframe, seriescolor = :black, camera = (10,30), margin = 8mm)
cd(plot_dir)
png(p,"3d_contour_non_dim_util_85.png")
cd(work_dir)
##
xlabel_real = L"$f_{CAM} (\mathrm{wt} \%)$"
ylabel_real = L"$\theta_{CAM} (\mathrm{wt} \%)$"
zlabel_real = L"$D_{SE} (\mu m)$"
title = L"$D_{CAM}^{BIG} = 15.0 \mu m, D_{CAM}^{SM} = 4.5 \mu m, D_{SE} = 2.5 \mu m$"
p = plot(xlabel = xlabel_real,
        ylabel = ylabel_real,
        zlabel = zlabel_real,
        xmirror = false,
        framestyle = :box,
        legend = :outertopright,
        legendfontsize = 14,
        legendtitlefontsize = 16,
        tickfontsize = 14,
        guidefontsize = 16,
        legendtitle = "Distribution",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        title = title,
        titlefontsize = 16,
        grid = false,
        xlim = (70,95))
brr = 0.3
bmr = 0.25
df2 = @where(df_bimodal, :particle_size_ratio .== 1.5, :bimodal_mix_ratio .== 0.25, :bimodal_radius_ratio .== 0.3)
df_uni = @where(df_unimodal, :particle_size_ratio .== 6.0)
plot!(df2[!,:wt_ratio]*100.0, df2[!,:util2]*100.0, label="Bimodal")
plot!(df_uni[!,:wt_ratio]*100.0, df_uni[!,:util2]*100.0, label="Unimodal")
cd(plot_dir)
png(p,"test3.png")
cd(work_dir)
xlabel_real = L"$f_{CAM} (\mathrm{wt} \%)$"
ylabel_real = L"$\theta_{CAM} (\mathrm{wt} \%)$"
zlabel_real = L"$D_{SE} (\mu m)$"
title = L"$D_{CAM}^{BIG} = 15.0 \mu m, D_{CAM}^{SM} = 4.5 \mu m, D_{SE} = 2.5 \mu m$"
p = plot(xlabel = xlabel_real,
        ylabel = ylabel_real,
        zlabel = zlabel_real,
        xmirror = false,
        framestyle = :box,
        legend = :outertopright,
        legendfontsize = 14,
        legendtitlefontsize = 16,
        tickfontsize = 14,
        guidefontsize = 16,
        legendtitle = "Distribution",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        title = "",
        titlefontsize = 16,
        grid = false,
        xlim = (70,95))
brr = 0.3
bmr = 0.25
df2 = @where(df_bimodal, :particle_size_ratio .== 1.5, :bimodal_mix_ratio .== 0.25, :bimodal_radius_ratio .== 0.3)
df_uni = @where(df_unimodal, :particle_size_ratio .== 6.0)
plot!(df2[!,:wt_ratio]*100.0, df2[!,:util2]*100.0, label="Bimodal")
plot!(df_uni[!,:wt_ratio]*100.0, df_uni[!,:util2]*100.0, label="Unimodal")
cd(plot_dir)
png(p,"test4.png")
cd(work_dir)
## New compare unimodal to bimodal2

Dcambig = 16.2
brr = 0.32
bmr = 0.25
Dcamsmall = Dcambig*brr
# unimodal util > 99
#bimodal util > 99
dfbi = @where(df_bimodal, :bimodal_radius_ratio .== 0.32, :bimodal_mix_ratio .== 0.25)
gd_wr = groupby(dfbi, :wt_ratio)
df_new = DataFrame([Float64, Float64, Float64, Float64],
    [:wt_ratio, :lambda_crit, :dcamsmall, :dse])

for (i,wrr) in enumerate(keys(gd_wr))
    df_wr = gd_wr[wrr]
    # println(brr.bimodal_radius_ratio)
    xinter = smooth2(df_wr[!,:particle_size_ratio], df_wr[!,:util2])
    println(xinter)
    push!(df_new, [wrr.wt_ratio, xinter,
            Dcamsmall,Dcamsmall/xinter])
end
sort!(df_new)
## Non dimensional 4d plots
pyplot()
xlabel_real = L"$f_{CAM} (\mathrm{wt} \%)$"
xlabel = L"$\lambda = \frac{D_{CAM}^{BIG}}{D_{SE}}$"
ylabel_real = L"$\theta_{CAM} (\mathrm{wt} \%)$"
zlabel_real = L"$D_{SE} (\mu m)$"
title = L"$f_{CAM} = 90 (\mathrm{wt} \%)$"
p = plot(xlabel = xlabel,
        ylabel = ylabel_real,
        zlabel = zlabel_real,
        xmirror = false,
        framestyle = :box,
        legend = :outertopright,
        legendfontsize = 14,
        legendtitlefontsize = 16,
        tickfontsize = 14,
        guidefontsize = 16,
        legendtitle = "",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        title = title,
        titlefontsize = 16,
        grid = false,
        ylim = (80,100),
        xlim = (2.5,5))
qbrr = [0.25,0.275, 0.3, 0.4]
qbmr = [0.25]

df2 = @where(df_bimodal, :wt_ratio .== 0.9)
gd = groupby(df2,:bimodal_mix_ratio)
alphas = [0.1, 0.2, 0.3, 0.4, 0.5]
k = 1
for (j,bmr) in enumerate(keys(gd))
    global k
    dfx = gd[bmr]
    if bmr.bimodal_mix_ratio in qbmr
        k += 1
        newyymin = []
        newyymax = []
        gd2 = groupby(dfx,:bimodal_radius_ratio)
        newxx1 = []

    end
end
plot!(collect(0:2.0:100.1), [98], seriestype=:hline, ls=:dash, color = :black, label="")
plot!(collect(0:2.0:100.1), [90], seriestype=:hline, ls=:dash, color = :black, label="")
# plot!(newxx, newyymin-yl)
# plot!(newxx, newyymin)
png(p, "testxx.png")
##
pyplot()
xlabel_real = L"$f_{CAM} (\mathrm{wt} \%)$"
xlabel = L"$\lambda = \frac{D_{CAM}^{BIG}}{D_{SE}}$"
ylabel_real = L"$\theta_{CAM} (\mathrm{wt} \%)$"
ylabel_real = L"$\rho (\%)$"
zlabel_real = L"$D_{SE} (\mu m)$"
title = L"$f_{CAM} = 90 (\mathrm{wt} \%)$"
p = plot(xlabel = xlabel,
        ylabel = ylabel_real,
        zlabel = zlabel_real,
        xmirror = false,
        framestyle = :box,
        legend = :outertopright,
        legendfontsize = 14,
        legendtitlefontsize = 16,
        tickfontsize = 14,
        guidefontsize = 16,
        legendtitle = "",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        title = title,
        titlefontsize = 16,
        grid = false,
        ylim = (60,90),
        xlim = (2.5,10))
qbrr = [0.25, 0.4]
qbmr = [0.2 0.4]

df2 = @where(df_bimodal, :wt_ratio .== 0.9, :bimodal_mix_ratio .== 0.25)
gd2 = groupby(df2,:bimodal_radius_ratio)
newxx = collect(1.0:0.01:6.0)
newxx2 = collect(0.5:0.01:6.0)
newyymax = []
newyymax = []
newyymin = []
newxxmin = []
df2 = @where(df_bimodal, :wt_ratio .== 0.9)
gd = groupby(df2,:bimodal_mix_ratio)
alphas = [0.1, 0.2, 0.3, 0.4, 0.5]
k = 1
for (j,bmr) in enumerate(keys(gd))
    global k
    dfx = gd[bmr]
    gd2 = groupby(dfx, :bimodal_radius_ratio)
    if bmr.bimodal_mix_ratio in qbmr
        for (i,brr) in enumerate(keys(gd2))
            global newxx, newyymin, newyymax, newxx2,newxxmax, newxxmin
            dfxx = gd2[brr]
            if brr.bimodal_radius_ratio in qbrr
                x = dfxx.particle_size_ratio
                y = dfxx.pack2
                dict = smooth_util(x,y)
                model = dict["model"]
                fit = dict["fit"]
                newyy = model(newxx2, coef(fit))
                # newxx = newxx/brr.bimodal_radius_ratio
                labstr = string(brr.bimodal_radius_ratio)
                if brr.bimodal_radius_ratio == 0.25
                    newyymin = newyy*100.0
                    newxxmin = newxx2/0.25
                    # plot!(dfxx.particle_size_ratio/brr.bimodal_radius_ratio, dfxx.util2*100.0)
                    plot!(newxx2./brr.bimodal_radius_ratio, newyy*100.0, color = :black, label="")
                end
                if brr.bimodal_radius_ratio == 0.4
                    newyymax = newyy*100.0
                    newxxmax = newxx2/0.4
                    plot!(newxx2./brr.bimodal_radius_ratio, newyy*100.0, color=:black, label="")
                    # plot!(dfxx.particle_size_ratio/brr.bimodal_radius_ratio, dfxx.util2*100.0)
                end

                # plot!(newxx, newyy*100.0, label=labstr)
            end
        end
    end
end
# plot!(newxx2, newyymin ewyymin], fillrange=[newyymax,90], alpha = 0.6, label="")
# plot!(newxxmin, newyymin, label="test", fillrange = newyymin)
# plot!(newxxmax, newyymax, label="test2")

# plot!(collect(0:2.0:100.1), [98], seriestype=:hline, ls=:dash, color = :black, label="")
plot!(collect(0:2.0:100.1), [90], seriestype=:hline, ls=:dash, color = :black, label="")
# plot!(newxx, newyymin-yl)
# plot!(newxx, newyymin)
cd(plot_dir)
png(p, "testxx_pack.png")
cd(work_dir)
