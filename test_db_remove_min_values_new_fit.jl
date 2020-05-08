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


function smooth_util(x,y)
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
    return Dict("fit"=>fit, "model"=>model, "x"=>collect(0.25:0.05:3.0), "inter"=>xinter)
end

function smooth_pack(x,y)
    @. model(x,p) = 0.86/(1.0 + exp(-p[1]*(x-p[2])))
    p0 = [1.0, 0.25]
    fit = curve_fit(model, x, y, p0)
    newyy = model(x, coef(fit))
    return Dict("fit"=>fit, "model"=>model, "x"=>collect(0.25:0.05:3.0))
end


function replacenan(x,y)
    if length(x) != length(y)
        error("Arrays have to be of equal length")
    end
    for i in eachindex(x)
        x[i] = ifelse(isnan(x[i]), y[i], x[i])
    end
end
#

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
##
nmc_wt_ratio1 = collect(0.7:0.01:0.92)
# lambda_min1 = [0.25,0.3,0.35,0.4,0.45, 0.5,0.55,0.6,
#     0.65,0.75,0.8,0.9, 1.0, 1.1, 1.2, 1.25,
#     1.4, 1.5, 1.6666667, 1.8, 1.9, 2.0, 2.5,3.0]
lambda_min1 = collect(0.25:0.025:3.0)
bimodal_radius_ratio2 = collect(0.2:0.005:0.4)
# bimodal_radius_ratio3 = collect(0.30:0.01:0.5)
# bimodal_radius_ratio1 = hcat(bimodal_radius_ratio2,bimodal_radius_ratio3)
bimodal_mix_ratio1 = bimodal_radius_ratio2
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
points=convert(Matrix,df_avg[:,[:wt_ratio, :bimodal_radius_ratio, :bimodal_mix_ratio, :particle_size_ratio]])
int_points=convert(Matrix,df_interp[:,[:wt_ratio, :bimodal_radius_ratio, :bimodal_mix_ratio, :particle_size_ratio]])
util = df_avg[!,:utilization]
pack = df_avg[!,:packing_density]
int_utilnear = si.griddata(points, util, int_points, method="nearest")
int_util = si.griddata(points, util, int_points, method="linear")
replacenan(int_util, int_utilnear)
int_packnear = si.griddata(points, pack, int_points, method="nearest")
int_pack = si.griddata(points, pack, int_points, method="linear")
replacenan(int_pack, int_packnear)
df_interp[!,:utilization] = int_util
df_interp[!,:packing_density] = int_pack
df_smooth= by(df_interp,[:wt_ratio,:bimodal_mix_ratio, :bimodal_radius_ratio],
    (:particle_size_ratio, :packing_density, :utilization)
    => x->(particle_size_ratio = x.particle_size_ratio,
           packing_density = x.packing_density,
           utlization = x.utilization,
           pack2 = smooth(x.particle_size_ratio, x.packing_density),
           util2 = smooth(x.particle_size_ratio, x.utilization)))
df_smooth= by(df_interp,[:wt_ratio,:bimodal_mix_ratio, :bimodal_radius_ratio],
   (:particle_size_ratio, :packing_density, :utilization)
    => x->(util2 = smooth_util(x.particle_size_ratio, x.utilization),
           pack2 = smooth_pack(x.particle_size_ratio, x.packing_density)))
# display(p)
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

df2 = @where(df_smooth, :bimodal_mix_ratio .== 0.25, :bimodal_radius_ratio .== 0.25)
Dcamsmall = Dcambig*0.25

x = get(df2[1,4], "x", nothing)
Dse = Dcamsmall ./ x
y = df2[!,1]*100.0


function get_contour(df, x)
    z = zeros(Float64, nrow(df2), length(x))
    for r in range(1,stop=nrow(df2))
        newx = get(df2[r,4],"x", nothing)
        fit = get(df2[r,4],"fit", nothing)
        model = get(df2[r,4], "model", nothing)
        newy = model(newx,coef(fit))
        z[r,:] = newy*100.0
    end
    return z
end
z = get_contour(df2, x)
pyplot()
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
        ylim = (70,90))
# plot(x,y,z, st=:contour, seriescolor=:black, levels=[95.0])
plot!(Dse,y, z, st=:contourf, seriescolor=cgrad(scheme.colors))
plot!(Dse,y, z, st=:contourf, seriescolor=cgrad(scheme.colors), levels=256)
plot!(Dse,y, z, st=:contour, seriescolor=:black, levels=[95.0,98.0])
cd(work_dir)
cd("julia")
png(p,"test_contour_real_psi_0.25_phi_0.25.png")
cd("..")
#End Plotting
## Linear Plots
pyplot()
sf = pyimport("scipy.ndimage.filters")
xlabel = L"$\lambda_{crit} = \frac{D_{CAM}^{SM}}{D_{SE}}_{crit}$"
ylabel = L"$f_{CAM} (wt \%)$"
ylabel2 = L"$\rho (\%)$"
p = plot(xlabel = ylabel,
        ylabel = xlabel,
        xmirror = false,
        framestyle = :box,
        legend = :outertopright,
        legendfontsize = 14,
        legendtitlefontsize = 18,
        tickfontsize = 16,
        guidefontsize = 18,
        legendtitle = L"$\phi = \frac{W_{CAM}^{SM}}{W_{CAM}}$",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        grid = false,
        xlim = (70,90))

q_brr = [0.25, 0.3, 0.4]
q_bmr = [0.2, 0.5]
marker = [:solid, :dash, :dot]
colors = [:red, :green, :blue]


for (i,bmr) in enumerate(q_bmr)
    df_new = DataFrame([Float64],[:wt_ratio])
    for w in y
        push!(df_new,[w])
    end
    dfxx = @where(df_smooth, :bimodal_mix_ratio .== bmr)
    gd_brr = groupby(dfxx,:bimodal_radius_ratio)
    j = 1
    for brr in keys(gd_brr)
        if brr.bimodal_radius_ratio in q_brr
            l_crit_95 = Float64[]
            yy = Float64[]
            println(brr.bimodal_radius_ratio)
            for r in range(1,stop=nrow(gd_brr[brr]))
                xinter_95 = get(gd_brr[brr][r,4], "inter", nothing)[2]
                push!(l_crit_95,xinter_95)
                push!(yy, gd_brr[brr][r,1]*100)
            # println("$brr $xinter_95")
            end
            labstr = ""
            if i == 1
                labstr = latexstring(string(brr.bimodal_radius_ratio))
            end
            # labstr = ifelse(i == 1,latexstring(string(brr.bimodal_radius_ratio),"")
            # println(l_crit_95)
            l_crit_smooth = sf.gaussian_filter1d(l_crit_95, sigma=0.9)
            # println(l_crit_smooth)
            # plot!(yy, l_crit_smooth, label=labstr, ls=marker[i], color = colors[j])
            j += 1
            x = Symbol("brr"*string(brr.bimodal_radius_ratio))
            df_new[!,x] = l_crit_smooth
        end
    end
    df_new[!,:mean] = mean.(eachrow(df_new[!,2:4]))
    df_new[!,:max] = maximum.(eachrow(df_new[!,2:4]))
    df_new[!,:min] = minimum.(eachrow(df_new[!,2:4]))
    # plot!(df_new[!,:wt_ratio],df_new[!,Symbol("brr0.3")])
    ymean = df_new[!,:mean]
    yl = df_new[!,:min]-df_new[!,:mean]
    yu = abs.(df_new[!,:mean]-df_new[!,:max])
    labs = string(bmr)
    plot!(df_new[!,:wt_ratio],df_new[!,Symbol("brr0.3")],linecolor=:black,ls=marker[i], label=labs)
    plot!(df_new[!,:wt_ratio], ymean, ribbon=[yl yu], label="", linecolor=nothing, fillalpha=0.5, fillcolor = colors[i])
end
cd("julia")
png(p,"test_real_phi.png")
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
        legendtitle = L"$f_{CAM} (wt \%)$",
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        grid = false)
        # inset = (1, bbox(0.05,0.3, 0.4, 0.4, :bottom, :right)))
wr = [0.75, 0.85, 0.9]
colors = [:blue, :orange, :green, :black, :red]
for (ii, wrr) in enumerate(wr)
    dfx = @where(df_interp, :wt_ratio .== wrr)
    gd_bmr = groupby(dfx, :bimodal_mix_ratio)
    df_bmr = get(gd_bmr, (bimodal_mix_ratio=0.25,),nothing)
    gd_brr = groupby(df_bmr,:bimodal_radius_ratio)
    # q_brr = [0.21, 0.275, 0.285, 0.3, 0.4]
    q_brr = [0.3]
    j = 1
    @. model(x,p) = 1.0/(1.0 + exp(-p[1]*(x-p[2])))
    # @. model2(x,p) = 0.86/(1.0 + exp(-p[1]*(x-p[2])))
    newxx = collect(0.1:0.01:3.0)
    p0 = [1.0, 0.25]
    for brr in keys(gd_brr)
        # global j, colors
        if brr.bimodal_radius_ratio in q_brr
            xdata = gd_brr[brr][!,:particle_size_ratio]
            ydata = gd_brr[brr][!,:utilization]
            # ydata2 = gd_brr[brr][!,:packing_density]
            fit = nothing
            # fit2 = nothing
            println("computing ...", brr.bimodal_radius_ratio)
            # println(brr)
            try
                fit = curve_fit(model, xdata, ydata, p0)
                # fit2 = curve_fit(model2, xdata, ydata2, p0)
                println("computed complete", brr.bimodal_radius_ratio)
            catch
            end

            if !isnothing(fit)
                newyy = model(newxx, coef(fit))
                # newyy2 = model2(newxx, coef(fit2))
                kstr = string(wrr)
                labstr = latexstring(kstr)
                plot!(xdata, ydata .* 100.0, seriestype=:scatter, labels="", color = colors[ii], marker=:x)
                plot!(newxx, newyy .* 100.0, label=labstr, color = colors[ii])
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
png(p, "test_util.png")
cd("../")
##pyplot()

xlabel = L"$\lambda = \frac{D_{CAM}^{SM}}{D_{SE}}$"
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
        ylim = (65,88))
wr = [0.75, 0.85, 0.9]
colors = [:blue, :orange, :green, :black, :red]
for (ii, wrr) in enumerate(wr)
    dfx = @where(df_interp, :wt_ratio .== wrr)
    gd_bmr = groupby(dfx, :bimodal_mix_ratio)
    df_bmr = get(gd_bmr, (bimodal_mix_ratio=0.25,),nothing)
    gd_brr = groupby(df_bmr,:bimodal_radius_ratio)
    # q_brr = [0.21, 0.275, 0.285, 0.3, 0.4]
    q_brr = [0.3]
    j = 1
    @. model(x,p) = 0.86/(1.0 + exp(-p[1]*(x-p[2])))
    newxx = collect(0.1:0.01:3.0)
    p0 = [1.0, 0.25]
    for brr in keys(gd_brr)
        # global j, colors
        if brr.bimodal_radius_ratio in q_brr
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
                plot!(xdata, ydata .* 100.0, seriestype=:scatter, labels="", color = colors[ii], marker=:x)
                plot!(newxx, newyy .* 100.0, label=labstr, color = colors[ii])
                # println(kstr)
                j += 1
            end
            println(".............................")
        end
    end
end
cd("julia")
plot!(collect(0.25:2.0:100), [98], seriestype=:hline,color = :black, label="")
png(p, "test_packing.png")
cd("../")
##
