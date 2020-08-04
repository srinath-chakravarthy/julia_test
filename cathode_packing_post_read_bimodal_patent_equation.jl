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
            # if isfile("cathode.json")
            #     data2 = JSON.parsefile("cathode.json"; dicttype=Dict, inttype=Int64, use_mmap=true)
            #     if haskey(data2, "Results")
            #         for b in data2["Results"]
            #             wt2_py = parse(Float64,b["Active_Weight_fraction"])
            #         end
            #     end
            # end
            return wt, wt1, wt2
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
plot_dir = "/home/srinath/Projects/cathode_packing_data/bimodal2/paper_plots/"
bimodal_dir="/home/srinath/Projects/cathode_packing_data/bimodal2/new_data/"
cd(bimodal_dir)
df_bimodal = DataFrame(CSV.File("all_data_interpolated_smoothed.csv"))
cd(work_dir)

## Plotting
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

cdict4_reverse = Dict(:red =>  ((0.0, 0.0, 0.0),
                   (0.2, 0.0, 0.0),
                   (0.4, 0.6, 0.6),
                   (0.5, 0.0, 0.0),
                   (0.9, 0.9, 0.9),
                   (1.0, 1.0, 1.0)),

         :green => ((0.0, 1.0, 1.0),
                   (0.4, 0.5, 0.5),
                   (0.95, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         :blue =>  ((0.0, 1.0, 1.0),
                   (0.5, 0.0, 0.0),
                   (0.6, 1.0, 1.0),
                   (0.8, 1.0, 1.0),
                   (1.0, 1.0, 1.0)))

scheme2 = make_colorscheme(cdict4)
scheme_rev = make_colorschem(cdict4_rev)

df2 = @where(df_bimodal, :wt_ratio .== 0.85,
            :bimodal_radius_ratio .> 0.245, :bimodal_radius_ratio .<0.405,
            :bimodal_mix_ratio .>0.195, :bimodal_mix_ratio .<0.505,
            :particle_size_ratio .> 0.75)
df_new = DataFrame([Float64, Float64, Float64, Float64, Float64],
    [:bimodal_radius_ratio,:bimodal_mix_ratio, :lambda_crit, :dcamsmall, :dse])
gd = groupby(df2, :bimodal_radius_ratio)
Dcambig = 16.2
sf = pyimport("scipy.ndimage.filters")
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
