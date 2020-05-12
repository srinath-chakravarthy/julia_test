using DataFrames, DataFramesMeta, CSV, Statistics
using Colors, ColorSchemes, ColorSchemeTools
using JSON
using Plots
using PyCall, Conda, LaTeXStrings
using LsqFit
import GR.meshgrid
using Contour
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

function smooth2(x, y)
    @. model(x,p) = p[1]/(p[2] + exp(-p[3]*(x-p[4])))
    p0 = [0.5, 0.5, 1.0, 0.25]
    fit = curve_fit(model, x, y, p0)
    newyy = model(x, coef(fit))
    return newyy
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
ylabel = L"$\theta_{CAM} (wt \%)$"
ylabel2 = L"$\rho (\%)$"
xlabel_real = L"$D_{SE} (\mu m)$"
ylabel_real = ylabel
title = L"$\psi =  0.285, \phi = 0.25$"

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
        ylim = (90,101),
        xlim = (1,5))
Dcambig = 15.0
qwrr = [0.85, 0.9]
brr = 0.285
bmr = 0.25
Dcamsmall = Dcambig * brr
dcams = round(Dcamsmall,digits=2)
colors = [:blue, :orange]
title_real = latexstring("\$D_{CAM}^{BIG} = $(Dcambig) \\mu m, D_{CAM}^{SM} = $(dcams) \\mu m, \\phi = $(bmr)\$")
for (i,wrr) in enumerate(qwrr)
    df_uni = @where(df_unimodal, :wt_ratio .== wrr)
    df_bi = @where(df_bimodal, :wt_ratio .== wrr, :bimodal_radius_ratio .== brr,
                    :bimodal_mix_ratio .== bmr)
    labstr = string(wrr*100.0)
    plot!(Dcambig./df_uni[!,:particle_size_ratio], df_uni[!,:util2]*100.0, ls=:dash, color = colors[i], label=labstr)
    plot!(Dcamsmall./df_bi[!,:particle_size_ratio], df_bi[!,:util2]*100.0, ls=:solid, color = colors[i], label="")
end
plot!(collect(0:2.0:100), [99.5], seriestype=:hline,color = :black, label="")
cd(plot_dir)
png(p, "compare_uni_bi_1.png")
cd(work_dir)
## Contour Plots
using Contour
sf = pyimport("scipy.ndimage.filters")
xlabel_real = L"$f_{CAM} (\mathrm{wt} \%)$"
ylabel_real = L"$D_{SE} (\mu m)$"
brr = 0.333333
Dcamsmall = Dcambig * brr
dcams = round(Dcamsmall,digits=1)
title_real = latexstring("\$D_{CAM}^{BIG} = $(Dcambig) \\mu m, D_{CAM}^{SM} = $(dcams) \\mu m, \\phi = $(bmr)\$")
p = plot(xlabel = ylabel_real,
        ylabel = xlabel_real,
        xmirror = false,
        framestyle = :box,
        legend = :outertopright,
        legendfontsize = 14,
        legendtitlefontsize = 16,
        tickfontsize = 14,
        guidefontsize = 16,
        # legendtitle = "Distribution",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        title = title_real,
        titlefontsize = 16,
        grid = false)

x = unique(df_unimodal[!,:particle_size_ratio])
y = unique(df_unimodal[!,:wt_ratio])*100.0
z = df_unimodal[!,:util2]*100.0
z2 = reshape(z, length(x), length(y))
zsmooth = sf.gaussian_filter(z2,sigma=3.0)
contours_uni = contours(x,y,zsmooth,[99.5])
cl = first(Contour.levels(contours_uni))
line = first(Contour.lines(cl))

xs, ys = Contour.coordinates(line)
newys = smooth2(xs,ys/100.0) * 100.0
# plot!(x,y, zsmooth', st=:contour, lt= :dash, seriescolor=:black, levels=[98.99,99.0], xlim = (3,10))
# plot!(newys, Dcambig./xs, ls = :dash, label="Unimodal", color = :black)
xs_save = xs
ys_save = newys

# plot!(ys, Dcambig./xs, seriestype = :scatter, marker=:x)

bmr = 0.25
brr = 0.3
dcams = Dcambig * brr
df_bi = @where(df_bimodal, :bimodal_radius_ratio .== brr,
                :bimodal_mix_ratio .== bmr)
x = unique(df_bi[!,:particle_size_ratio])
y = unique(df_bi[!,:wt_ratio])*100.0
z = df_bi[!,:util2]*100.0
z2 = reshape(z, length(x), length(y))
zsmooth = sf.gaussian_filter(z2,sigma=0.1)
contours_uni = contours(x,y,zsmooth,[99.5])
cl = first(Contour.levels(contours_uni))
line = first(Contour.lines(cl))

xs, ys = Contour.coordinates(line)
newys = smooth2(xs,ys/100.0) * 100.0

plot!(dcams./x,y, zsmooth', st=:contourf, seriescolor=cgrad(scheme.colors))
plot!(dcams./x,y, zsmooth', st=:contourf, seriescolor=cgrad(scheme.colors), levels=256)
# plot!(x/brr,y, zsmooth', st=:contour, seriescolor=:black, levels=[98.99,99.0])
plot!(dcams./xs, newys, ls = :solid, color = :black,xlim=(1.5,5), ylim=(70,92), label="")
plot!(Dcambig./xs_save, ys_save, color = :black, ls=:dash, xlim=(1.5,5), ylim=(70,92), label="")
#
# bmr = 0.25
# brr = 0.3
# dcams = Dcambig * brr
# df_bi = @where(df_bimodal, :bimodal_radius_ratio .== brr,
#                 :bimodal_mix_ratio .== bmr)
# x = unique(df_bi[!,:particle_size_ratio])
# y = unique(df_bi[!,:wt_ratio])*100.0
# z = df_bi[!,:utilization]*100.0
# z2 = reshape(z, length(x), length(y))
# zsmooth_uni = sf.gaussian_filter(z2,sigma=5.0)
# contours_uni = contours(x,y,zsmooth_uni,[99.0])
# cl = first(Contour.levels(contours_uni))
# line = first(Contour.lines(cl))
#
# xs, ys = Contour.coordinates(line)
# plot!(ys, dcams./xs, ls = :solid, color = :red)
#
# bmr = 0.25
# brr = 0.4
# dcams = Dcambig * brr
# df_bi = @where(df_bimodal, :bimodal_radius_ratio .== brr,
#                 :bimodal_mix_ratio .== bmr)
# x = unique(df_bi[!,:particle_size_ratio])
# y = unique(df_bi[!,:wt_ratio])*100.0
# z = df_bi[!,:utilization]*100.0
# z2 = reshape(z, length(x), length(y))
# zsmooth_uni = sf.gaussian_filter(z2,sigma=5.0)
# contours_uni = contours(x,y,zsmooth_uni,[99.0])
# cl = first(Contour.levels(contours_uni))
# line = first(Contour.lines(cl))
#
# xs, ys = Contour.coordinates(line)
# plot!(ys, dcams./xs, ls = :solid, color = :green)


# bmr = 0.25
# df_bi = @where(df_bimodal, :bimodal_radius_ratio .== brr,
#                 :bimodal_mix_ratio .== bmr)
# x = unique(df_bi[!,:particle_size_ratio])
# y = unique(df_bi[!,:wt_ratio])*100.0
# z = df_bi[!,:utilization]*100.0
# z2 = reshape(z, length(x), length(y))
# zsmooth_uni = sf.gaussian_filter(z2,sigma=4.0)
# contours_uni = contours(x,y,zsmooth_uni,[99.0])
# cl = first(Contour.levels(contours_uni))
# line = first(Contour.lines(cl))
#
# xs, ys = Contour.coordinates(line)
# plot!(ys, dcams./xs, ls = :solid, label="Bimodal", color = :blue)
#
# bmr = 0.3
# df_bi = @where(df_bimodal, :bimodal_radius_ratio .== brr,
#                 :bimodal_mix_ratio .== bmr)
# x = unique(df_bi[!,:particle_size_ratio])
# y = unique(df_bi[!,:wt_ratio])*100.0
# z = df_bi[!,:utilization]*100.0
# z2 = reshape(z, length(x), length(y))
# zsmooth_uni = sf.gaussian_filter(z2,sigma=5.0)
# contours_uni = contours(x,y,zsmooth_uni,[99.0])
# cl = first(Contour.levels(contours_uni))
# line = first(Contour.lines(cl))
#
# xs, ys = Contour.coordinates(line)
# plot!(ys, dcams./xs, ls = :solid, label="Bimodal", color = :red)
#
# bmr = 0.4
# df_bi = @where(df_bimodal, :bimodal_radius_ratio .== brr,
#                 :bimodal_mix_ratio .== bmr)
# x = unique(df_bi[!,:particle_size_ratio])
# y = unique(df_bi[!,:wt_ratio])*100.0
# z = df_bi[!,:utilization]*100.0
# z2 = reshape(z, length(x), length(y))
# zsmooth_uni = sf.gaussian_filter(z2,sigma=5.0)
# contours_uni = contours(x,y,zsmooth_uni,[99.0])
# cl = first(Contour.levels(contours_uni))
# line = first(Contour.lines(cl))
#
# xs, ys = Contour.coordinates(line)
# plot!(ys, dcams./xs, ls = :solid, label="Bimodal", color = :green)

cd(plot_dir)
png(p, "compare_uni_bi_1.png")
cd(work_dir)
## Variation with bimodal radius ratio_
using Contour
p = plot(xlabel = xlabel_real,
        ylabel = ylabel_real,
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
