using DataFrames, DataFramesMeta, CSV, Statistics
using Colors, ColorSchemes, ColorSchemeTools
using JSON
using Plots
using PyCall, Conda, LaTeXStrings
using LsqFit
import GR.meshgrid
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
            a = data["Results"]
            for b in a
                wt = parse(Float64, b["Output_porosity"])
                wt1 = parse(Float64,b["Particle_fraction"])
                wt2 = parse(Float64,b["Active_Weight_fraction"])
            end
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
work_dir="/home/srinath/Projects/cathode_packing_data/bimodal_distribution_database"
cd(work_dir)
# cd("julia")

# df_avg = DataFrame(CSV.File("orig_data_julia_util_packing.csv"))
nstat = collect(1:5)
nmc_wt_ratio = collect(0.7:0.05:0.92)
# nmc_wt_ratio = [0.85, 0.9]
lambda_min = [0.25, 0.5, 0.75, 1.0, 1.1, 1.2, 1.25, 1.4, 1.5, 1.6666667, 1.8, 1.9, 2.0, 2.5,3.0]
bimodal_radius_ratio = [0.2, 0.25, 0.275, 0.285, 0.3, 0.4, 0.5]
bimodal_mix_ratio = [0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.3, 0.4, 0.5]
df = DataFrame([Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64],
    [:wt_ratio, :bimodal_radius_ratio,:bimodal_mix_ratio,
        :particle_size_ratio,:packing_density,:num_utilization, :utilization, :stat])
df_avg = DataFrame([Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64],
    [:wt_ratio, :bimodal_radius_ratio,:bimodal_mix_ratio,
        :particle_size_ratio,:packing_density,:num_utilization, :utilization, :stat])

for w in nmc_wt_ratio
    wt_directory = "wt_ratio_" * string(round(w, digits=3))
    if  !isdir(wt_directory)
        continue
    end
    cd(wt_directory)
    for brr in bimodal_radius_ratio
        brr_dir = "bimodal_radius_ratio_" * string(brr)
        if !isdir(brr_dir)
            continue
        end
        cd(brr_dir)
        for b1 in lambda_min
            bdir = "lambda_" * string(round(b1,digits=3))
            if !isdir(bdir)
                continue
            end
            cd(bdir)
            for bmr in bimodal_mix_ratio
                bmr_dir = "bimodal_mix_ratio_" * string(bmr)
                if !isdir(bmr_dir)
                    continue
                end
                cd(bmr_dir)
                ss = 0
                wt_avg = 0.0
                wt1_avg = 0.0
                wt2_avg = 0.0
                wt_list = Float64[]
                wt1_list = Float64[]
                wt2_list = Float64[]
                for nst in nstat
                    stat_dir = "stat_" * string(nst)
                    if !isdir(stat_dir)
                        continue
                    end
                    cd(stat_dir)
                    wt = 0.0
                    wt1 = 0.0
                    wt2 = 0.0
                    try
                        wt,wt1, wt2 = run_post("cathode.json")
                    catch
                        cd("../")
                        continue
                    end
                    println(w, brr, bmr, b1, ss)
                    ss += 1
                    wt_avg += wt
                    wt1_avg += wt1
                    wt2_avg += wt2
                    push!(wt_list, wt)
                    push!(wt1_list, wt1)
                    push!(wt2_list, wt2)
                    push!(df,[w,brr, bmr,b1,wt,wt1, wt2, ss])
                    cd("../")
                end #stat_dir
                if ss > 0
                    if ss > 1
                        wt_avg1 = convert(Float64,mean(sort(wt_list)[2:ss]))
                        wt1_avg1 = convert(Float64, mean(sort(wt1_list)[2:ss]))
                        wt2_avg1 = convert(Float64, mean(sort(wt2_list)[2:ss]))
                        push!(df_avg, [w,brr, bmr,b1,wt_avg1,
                                wt1_avg1,
                                wt2_avg1,
                                ss])
                    else
                        push!(df_avg, [w,brr, bmr,b1, wt_avg/ss,wt1_avg/ss, wt2_avg/ss, ss])
                    end
                end
                cd("../") #bmr_dir
            end

            cd("../") #lambda_dir
        end
        cd("../") #brr_dir
    end
    cd("../") #wt_ratio_dir
end
cd(work_dir)
## Now try exponential interpolation
df_smooth= by(df_avg,[:wt_ratio,:bimodal_mix_ratio, :bimodal_radius_ratio],
   (:particle_size_ratio, :packing_density, :utilization)
    => x->(particle_size_ratio = x.particle_size_ratio,
           packing_density = x.packing_density,
           utlization = x.utilization,
           util2 = smooth_util(x.particle_size_ratio, x.utilization, true),
           pack2 = smooth_pack(x.particle_size_ratio, x.packing_density, true)))
##
cd(work_dir)
cd("julia/SAIT_data/")
CSV.write("orig_data_avg_smooth_all.csv", df_smooth)
cd(work_dir)

##
nmc_wt_ratio1 = collect(0.7:0.01:0.9)
# lambda_min1 = [0.25,0.3,0.35,0.4,0.45, 0.5,0.55,0.6,
#     0.65,0.75,0.8,0.9, 1.0, 1.1, 1.2, 1.25,
#     1.4, 1.5, 1.6666667, 1.8, 1.9, 2.0, 2.5,3.0]
lambda_min1 = collect(0.75:0.025:3.0)
bimodal_radius_ratio2 = collect(0.2:0.01:0.4)
append!(bimodal_radius_ratio2, 0.275)
append!(bimodal_radius_ratio2, 0.285)
# bimodal_radius_ratio3 = collect(0.30:0.01:0.5)
# bimodal_radius_ratio1 = hcat(bimodal_radius_ratio2,bimodal_radius_ratio3)
bimodal_mix_ratio1 = sort(bimodal_radius_ratio2)
bimodal_radius_ratio1 = bimodal_radius_ratio2
bimodal_mix_ratio1 = collect(0.2:0.01:0.5)

df_interp = DataFrame([Float64, Float64, Float64, Float64, Float64, Float64],
    [:wt_ratio, :bimodal_radius_ratio,:bimodal_mix_ratio,
        :particle_size_ratio,:packing_density,:utilization])
for i in nmc_wt_ratio1
    for j in bimodal_radius_ratio1
        for k in bimodal_mix_ratio1
            for l in lambda_min1
                # println(i,j,k,l)
                push!(df_interp, [i, j, k, l, 0.0, 0.0])
            end
        end
    end
end
##
si = pyimport("scipy.interpolate")
points=convert(Matrix,df_smooth[:,[:wt_ratio, :bimodal_radius_ratio, :bimodal_mix_ratio, :particle_size_ratio]])
int_points=convert(Matrix,df_interp[:,[:wt_ratio, :bimodal_radius_ratio, :bimodal_mix_ratio, :particle_size_ratio]])
util = df_smooth[!,:util2]
pack = df_smooth[!,:pack2]
int_utilnear = si.griddata(points, util, int_points, method="nearest")
int_util = si.griddata(points, util, int_points, method="linear")
replacenan(int_util, int_utilnear)
int_packnear = si.griddata(points, pack, int_points, method="nearest")
int_pack = si.griddata(points, pack, int_points, method="linear")
replacenan(int_pack, int_packnear)
df_interp[!,:utilization] = int_util
df_interp[!,:packing_density] = int_pack
cd(work_dir)
cd("julia/SAIT_data")
wr = [0.85 0.9]
filebase = "all_data_smooth_interp_"
for wrr in wr
    dfx = @where(df_interp, :wt_ratio .== wrr)
    newf = filebase * string(wrr) * ".csv"
    CSV.write(newf, dfx)
end
cd(work_dir)
# df_interp = DataFrame(CSV.File("all_data_interp_julia_util_packing.csv"))
###
#Plotting
pyplot()
Dcambig = 14.0
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

df2 = @where(df_interp, :bimodal_mix_ratio .== 0.25, :bimodal_radius_ratio .== 0.285)
Dcamsmall = Dcambig*0.285

# x = get(df2[1,4], "x", nothing)
# Dse = Dcamsmall ./ x
# y = df2[!,1]*100.0
#
#
# function get_contour(df, x)
#     z = zeros(Float64, nrow(df2), length(x))
#     for r in range(1,stop=nrow(df2))
#         newx = get(df2[r,4],"x", nothing)
#         fit = get(df2[r,4],"fit", nothing)
#         model = get(df2[r,4], "model", nothing)
#         newy = model(newx,coef(fit))
#         z[r,:] = newy*100.0
#     end
#     return z
# end
# z = get_contour(df2, x)
# \frac{D_{CAM}^{SM}}{D_{CAM}^{BIG}} =
# \frac{W_{CAM}^{SM}{W_{CAM}}}
x = unique(df2[!,:particle_size_ratio])
y = unique(df2[!,:wt_ratio])*100.0
z = df2[!,:utilization]*100.0
z1 = reshape(z, length(x), length(y))
pyplot()
xlabel = L"$\lambda = \frac{D_{CAM}^{SM}}{D_{SE}}$"
ylabel = L"$f_{CAM} (wt \%)$"
ylabel2 = L"$\rho (\%)$"
xlabel_real = L"$D_{SE} (\mu m)$"
ylabel_real = ylabel
title = L"$\psi =  0.285, \phi = 0.25$"
title_real = L"D_{CAM}^{BIG} = 14 \mu m, D_{CAM}^{SM} = 4.0 \mu m, \phi = 0.25"
p = plot(xlabel = xlabel_real,
        ylabel = ylabel_real,
        xmirror = false,
        framestyle = :box,
        legend = :outertopright,
        legendfontsize = 14,
        legendtitlefontsize = 16,
        tickfontsize = 14,
        guidefontsize = 16,
        legendtitle = L"$\theta_{CAM}$",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        title = title_real,
        titlefontsize = 16,
        grid = false,
        ylim = (70,90))
# plot(x,y,z, st=:contour, seriescolor=:black, levels=[95.0])
plot!(Dcamsmall./x,y, z1', st=:contourf, seriescolor=cgrad(scheme.colors))
plot!(Dcamsmall./x,y, z1', st=:contourf, seriescolor=cgrad(scheme.colors), levels=256)
plot!(Dcamsmall./x,y, z1', st=:contour, seriescolor=:black, levels=[95.0,98.0])
cd(work_dir)
cd("julia/SAIT_data/")
png(p,"contour_psi_0.285_phi_0.25_new_real.png")
cd(work_dir)
#End Plotting
## Linear Plots
pyplot()
sf = pyimport("scipy.ndimage.filters")
# xlabel = L"$\lambda_{crit} = \frac{D_{CAM}^{SM}}{D_{SE}}_{crit}$"
ylabel = L"$D_{SE} (\mu m)$"
xlabel = L"$\phi (wt \%)$"
ylabel2 = L"$\rho (\%)$"
p = plot(xlabel = xlabel,
        ylabel = ylabel,
        xmirror = false,
        framestyle = :box,
        legend = :outertopright,
        legendfontsize = 14,
        legendtitlefontsize = 18,
        tickfontsize = 14,
        guidefontsize = 16,
        legendtitle = L"$D_{CAM}^{SM} (\mu m)$",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        title = L"$D_{CAM}^{BIG} = 14 (\mu m)$",
        titlefontsize = 18,
        grid = false,
        xlim = (0.2,0.4))
Dcambig = 14.0

wr = [0.85, 0.9]
ls = [:solid, :dash]
q_brr = [0.25, 0.285, 0.3, 0.35]
colors = [:red, :blue, :green, :black, :purple, :orange]
for (i,wrr) in enumerate(wr)
    for (j,brr) in enumerate(q_brr)
        dcamsmall = 14.0*brr
        dfxx = @where(df_interp, :wt_ratio .== wrr, :bimodal_radius_ratio .== brr)
        gd_bmr = groupby(dfxx, :bimodal_mix_ratio)
        xbmr = []
        lcrit = []
        for bmr in keys(gd_bmr)
            push!(xbmr, bmr.bimodal_mix_ratio)
            # --- Find intersection by fitting
            dict = smooth_util(gd_bmr[bmr][!,:particle_size_ratio], gd_bmr[bmr][!,:utilization])
            xinter = get(dict,"inter", nothing)[2]
            push!(lcrit,xinter)
        end
        lcrit_smooth = sf.gaussian_filter1d(lcrit,sigma=1.0)
        if i == 1
            labstr = string(round(dcamsmall, digits=1))
        else
            labstr = ""
        end
        plot!(xbmr, dcamsmall./lcrit_smooth, label=labstr, color=colors[j], ls = ls[i])
    end
end

cd(work_dir)
cd("julia/SAIT_data/")
png(p,"D_crit_variation_with_bmr_brr.png")
cd("..")
##

# plot!([df_new[!,:mean] df_new[!,:mean]], fillrange= [df_new[!,:mean]-df_new[!,:min] df_new[!,:mean]-df_new[!,:min]])
# plot!(ymean, ribbon=[yl yu], fillalpha = 0.3, labels=["y","",""])
# cd("julia")
# png(p,"test_phi.png")
# cd("..")

##
pyplot()

xlabel = L"$\lambda = \frac{D_{CAM}^{SM}}{D_{SE}}$"
xlabel_real = L"$D_{SE} (\mu m)$"
ylabel = L"$\theta_{CAM} (\%)$"
ylabel_real = ylabel
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
        legendtitle = L"$f_{CAM} (wt \%)$",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        grid = false,
        xlim = (0.25,3.0))
        # inset = (1, bbox(0.05,0.3, 0.4, 0.4, :bottom, :right)))
wr = [0.85, 0.9]
colors = [:blue, :orange, :green, :black, :red]
bmr = 0.25
brr = 0.4
for (ii, wrr) in enumerate(wr)
    global bmr, brr
    dfx = @where(df_interp, :wt_ratio .== wrr, :bimodal_mix_ratio .== bmr,
                :bimodal_radius_ratio .== brr)
    dfx2 = @where(df_avg, :wt_ratio .== wrr, :bimodal_mix_ratio .== bmr,
                :bimodal_radius_ratio .== brr)
    labstr = string(wrr*100.0)
    # plot!(dfx[!,:particle_size_ratio], dfx[!,:utilization]*100.0, label=labstr, color=colors[ii])
    plot!(dfx2[!,:particle_size_ratio], dfx2[!,:utilization]*100.0, st=:scatter, label="", color=colors[ii], marker=:x)
    f = smooth_util(dfx2[!,:particle_size_ratio], dfx2[!,:utilization])
    println(wrr)
    # println(f)
    fit = get(f,"fit",nothing)
    model = get(f,"model",nothing)
    newx = collect(0.25:0.05:3.0)
    newy = model(newx, coef(fit))
    plot!(newx, newy*100.0,label=labstr, color=colors[ii])
end
cd(work_dir)
cd("julia/SAIT_data/")
plot!(collect(0.25:2.0:100), [98], seriestype=:hline,color = :black, label="")
png(p, "util_psi_0.3_phi_0.25.png")
cd(work_dir)
##pyplot()

# xlabel = L"$\lambda = \frac{D_{CAM}^{SM}}{D_{SE}}$"
xlabel = L"$D_{SE} (\mu m)$"
ylabel = L"$\rho (\%)$"
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
        legendtitle = L"$f_{CAM}$",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        grid = false,
        ylim = (80,88),
        xlim = (2.0, 4.0))
wr = [0.75, 0.85, 0.9]
colors = [:blue, :orange, :green, :black, :red]
for (ii, wrr) in enumerate(wr)
    dfx = @where(df_interp, :wt_ratio .== wrr)
    gd_bmr = groupby(dfx, :bimodal_mix_ratio)
    df_bmr = get(gd_bmr, (bimodal_mix_ratio=0.25,),nothing)
    gd_brr = groupby(df_bmr,:bimodal_radius_ratio)
    # q_brr = [0.21, 0.275, 0.285, 0.3, 0.4]
    q_brr = [0.285]
    j = 1
    @. model(x,p) = 0.86/(1.0 + exp(-p[1]*(x-p[2])))
    newxx = collect(0.25:0.01:3.0)
    p0 = [1.0, 0.25]
    for brr in keys(gd_brr)
        # global j, colors
        if brr.bimodal_radius_ratio in q_brr
            dcamsmall = Dcambig * brr.bimodal_radius_ratio
            xdata = gd_brr[brr][!,:particle_size_ratio]
            ydata = gd_brr[brr][!,:packing_density]
            fit = nothing
            println("computing ...", brr.bimodal_radius_ratio)
            # println(brr)
            try
                fit = curve_fit(model, xdata, ydata, p0)
                println("computed complete", brr.bimodal_radius_ratio)
            catch
            end

            if !isnothing(fit)
                newyy = model(newxx, coef(fit))
                kstr = string(wrr)
                labstr = latexstring(kstr)
                plot!(dcamsmall./xdata, ydata .* 100.0, seriestype=:scatter, labels="", color = colors[ii], marker=:x)
                plot!(dcamsmall./newxx, newyy .* 100.0, label=labstr, color = colors[ii])
                # println(kstr)
                j += 1
            end
            println(".............................")
        end
    end
end
cd(work_dir)
cd("julia")
plot!(collect(0.25:2.0:100), [98], seriestype=:hline,color = :black, label="")
png(p, "test_real_packing_psi_0.285.png")
cd("../")
##
# xlabel = L"$\lambda = \frac{D_{CAM}^{SM}}{D_{SE}}$"
xlabel = L"$D_{SE} (\mu m)$"
ylabel = L"$\theta_{CAM} (\%)$"
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
        legendtitle = L"$\phi$",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        grid = false,
        xlim = (0.25,3.0),
        ylim = (40.0, 100.0))
        # inset = (1, bbox(0.05,0.3, 0.4, 0.4, :bottom, :right)))
wrr = 0.9
brr = 0.285
dcamsmall = Dcambig * brr
dfx = @where(df_interp, :wt_ratio .== wrr, :bimodal_radius_ratio .== 0.285)
# dfx2 = @where(df_smooth, :wt_ratio .== wrr, :bimodal_radius_ratio .== 0.285)
gd_bmr = groupby(dfx, :bimodal_mix_ratio)
# gd_bmr2 = groupby(dfx2, :bimodal_mix_ratio)
q_bmr = [0.2,  0.25,  0.3, 0.4, 0.5]
# q_bmr =[ 0.25]
colors = [:red, :green, :blue, :orange, :purple, :black, :yellow, :pink]
j = 1
for bmr in keys(gd_bmr)
    global j
    if bmr.bimodal_mix_ratio in q_bmr
        df_bmr = gd_bmr[bmr]
        newx = df_bmr[!,:particle_size_ratio]
        newy = df_bmr[!,:utilization]
        newy_smooth = sf.gaussian_filter1d(newy, sigma=20.0)
        labstr = string(bmr.bimodal_mix_ratio)
        plot!(newx, newy_smooth*100.0, color= colors[j], label=labstr)
        # plot!(newx2, newy2*100.0, color= colors[j],label="")
        j += 1
    end
end
cd(work_dir)
cd("julia")
png(p,"util_variation_real_psi_0.285_w_0.85_interp.png")
cd("..")
##
##
# xlabel = L"$\lambda = \frac{D_{CAM}^{SM}}{D_{SE}}$"
xlabel = L"$D_{SE} (\mu m)$"
ylabel = L"$\theta_{CAM} (\%)$"
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
        legendtitle = L"$\psi$",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        grid = false,
        xlim = (0.25,3.0),
        ylim = (40.0, 100.0))
        # inset = (1, bbox(0.05,0.3, 0.4, 0.4, :bottom, :right)))
wrr = 0.9
bmr = 0.25
dcamsmall = Dcambig * brr
dfx = @where(df_interp, :wt_ratio .== wrr, :bimodal_mix_ratio .== bmr)
# dfx2 = @where(df_smooth, :wt_ratio .== wrr, :bimodal_radius_ratio .== 0.285)
gd_brr = groupby(dfx, :bimodal_mix_ratio)
# gd_bmr2 = groupby(dfx2, :bimodal_mix_ratio)
q_brr = [0.2,  0.25,  0.285, 0.3, 0.4]
# q_bmr =[ 0.25]
colors = [:red, :green, :blue, :orange, :purple, :black, :yellow, :pink]
j = 1
for brr in keys(gd_brr)
    global j
    if bmr.bimodal_radius_ratio in q_brr
        df_brr = gd_brr[brr]
        newx = df_brr[!,:particle_size_ratio]
        newy = df_brr[!,:utilization]
        newy_smooth = sf.gaussian_filter1d(newy, sigma=20.0)
        labstr = string(brr.bimodal_radius_ratio)
        plot!(newx, newy_smooth*100.0, color= colors[j], label=labstr)
        # plot!(newx2, newy2*100.0, color= colors[j],label="")
        j += 1
    end
end
cd(work_dir)
cd("julia")
png(p,"util_variation_real_phi_0.25_w_0.9_interp.png")
cd("..")
###Plotting
pyplot()
Dcambig = 14.0
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

df2 = @where(df_interp, :bimodal_mix_ratio .== 0.25, :bimodal_radius_ratio .== 0.285)
Dcamsmall = Dcambig*0.285

# x = get(df2[1,4], "x", nothing)
# Dse = Dcamsmall ./ x
# y = df2[!,1]*100.0
#
#
# function get_contour(df, x)
#     z = zeros(Float64, nrow(df2), length(x))
#     for r in range(1,stop=nrow(df2))
#         newx = get(df2[r,4],"x", nothing)
#         fit = get(df2[r,4],"fit", nothing)
#         model = get(df2[r,4], "model", nothing)
#         newy = model(newx,coef(fit))
#         z[r,:] = newy*100.0
#     end
#     return z
# end
# z = get_contour(df2, x)
# \frac{D_{CAM}^{SM}}{D_{CAM}^{BIG}} =
# \frac{W_{CAM}^{SM}{W_{CAM}}}
x = unique(df2[!,:particle_size_ratio])
y = unique(df2[!,:wt_ratio])*100.0
z = df2[!,:packing_density]*100.0
z1 = reshape(z, length(x), length(y))
pyplot()
xlabel = L"$\lambda = \frac{D_{CAM}^{SM}}{D_{SE}}$"
ylabel = L"$f_{CAM} (wt \%)$"
ylabel2 = L"$\rho (\%)$"
xlabel_real = L"$D_{SE} (\mu m)$"
ylabel_real = ylabel
title = L"$\psi =  0.285, \phi = 0.25$"
title_real = L"D_{CAM}^{BIG} = 14 \mu m, D_{CAM}^{SM} = 4.0 \mu m, \phi = 0.25"
p = plot(xlabel = xlabel,
        ylabel = ylabel,
        xmirror = false,
        framestyle = :box,
        legend = :outertopright,
        legendfontsize = 14,
        legendtitlefontsize = 16,
        tickfontsize = 14,
        guidefontsize = 16,
        legendtitle = L"$\theta_{CAM}$",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        title = title,
        titlefontsize = 16,
        grid = false,
        ylim = (70,90))
# plot(x,y,z, st=:contour, seriescolor=:black, levels=[95.0])
plot!(x,y, z1', st=:contourf, seriescolor=cgrad(scheme.colors))
plot!(x,y, z1', st=:contourf, seriescolor=cgrad(scheme.colors), levels=256)
# plot!(Dcamsmall./x,y, z1', st=:contour, seriescolor=:black, levels=[95.0,98.0])
cd(work_dir)
cd("julia/SAIT_data/")
png(p,"contour_packing_psi_0.285_phi_0.25_new_ratio.png")
cd(work_dir)
##
pyplot()

xlabel = L"$\lambda = \frac{D_{CAM}^{SM}}{D_{SE}}$"
xlabel_real = L"$D_{SE} (\mu m)$"
ylabel = L"$\theta_{CAM} (\%)$"
ylabel_real = ylabel
ylabel2 = L"$\rho (\%)$"
title_real = L"D_{CAM}^{BIG} = 14 \mu m, D_{CAM}^{SM} = 4.0 \mu m"
p = plot(xlabel = xlabel_real,
        ylabel = ylabel_real,
        xmirror = false,
        framestyle = :box,
        legend = :outertopright,
        legendfontsize = 12,
        legendtitlefontsize = 14,
        tickfontsize = 14,
        guidefontsize = 16,
        legendtitle = L"$\phi (wt \%)$",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        title = title_real,
        titlefontsize = 16,
        grid = false,
        xlim = (1.3,3.0),
        ylim = (95,100))
        # inset = (1, bbox(0.05,0.3, 0.4, 0.4, :bottom, :right)))
wrr = 0.9
colors = [:blue, :orange, :green, :black, :red]
Dcambig = 14.0
brr = 0.285
dcamsmall = Dcambig*brr
dfx = @where(df_interp, :wt_ratio .== wrr,:bimodal_radius_ratio .== brr)
gd_bmr = groupby(dfx, :bimodal_mix_ratio)

for bmr in keys(gd_bmr)
    if bmr.bimodal_mix_ratio > 0.23 && bmr.bimodal_mix_ratio <0.40
        df_bmr = gd_bmr[bmr]
        println(bmr.bimodal_mix_ratio)
        labstr = string(bmr.bimodal_mix_ratio)
        plot!(dcamsmall./df_bmr[!,:particle_size_ratio], df_bmr[!,:utilization]*100.0,label=labstr)
    end
end
plot!(collect(0.25:2.0:100), [98], seriestype=:hline,color = :black, label="")
cd(work_dir)
cd("julia/SAIT_data/")
png(p,"test_f_0.90.png")
cd(work_dir)
##
pyplot()
sf = pyimport("scipy.ndimage.filters")

xlabel = L"$\lambda = \frac{D_{CAM}^{SM}}{D_{SE}}$"
ylabel = L"$\phi = \frac{W_{CAM}^{SM}}{W_{CAM}} (wt \%)$"
ylabel2 = L"$\rho (\%)$"
xlabel_real = L"$D_{SE} (\mu m)$"
ylabel_real = ylabel
title = L"$\psi =  0.285, \phi = 0.25$"
title_real = L"D_{CAM}^{BIG} = 14 \mu m, D_{CAM}^{SM} = 4.0 \mu m"
p = plot(xlabel = xlabel_real,
        ylabel = ylabel_real,
        xmirror = false,
        framestyle = :box,
        legend = :outertopright,
        legendfontsize = 14,
        legendtitlefontsize = 16,
        tickfontsize = 14,
        guidefontsize = 16,
        legendtitle = L"$\theta_{CAM}$",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        title = title_real,
        titlefontsize = 16,
        grid = false,
        ylim = (0.2,0.4))
        # inset = (1, bbox(0.05,0.3, 0.4, 0.4, :bottom, :right)))
wrr = 0.85
colors = [:blue, :orange, :green, :black, :red]
Dcambig = 14.0
brr = 0.285
dcamsmall = Dcambig*brr
dfx = @where(df_interp, :wt_ratio .== wrr,:bimodal_radius_ratio .== brr)
# gd_bmr = groupby(dfx, :bimodal_mix_ratio)
x = unique(dfx[!,:particle_size_ratio])
y = unique(dfx[!,:bimodal_mix_ratio])
z = dfx[!,:utilization]*100.0
z1 = reshape(z, length(x), length(y))
z1smooth = sf.gaussian_filter(z1,sigma=2.0)
plot!(Dcamsmall./x,y, z1smooth', st=:contourf, seriescolor=cgrad(scheme.colors))
plot!(Dcamsmall./x,y, z1smooth', st=:contourf, seriescolor=cgrad(scheme.colors), levels=256)
plot!(Dcamsmall./x,y, z1smooth', st=:contour, seriescolor=:black, levels=[98.0])

cd(work_dir)
cd("julia/SAIT_data/")
png(p, "test_new_contour_0.85.png")
cd(work_dir)
