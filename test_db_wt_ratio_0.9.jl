using DataFrames, DataFramesMeta, CSV, Statistics
using JSON
using Plots
using PyCall, Conda, LaTeXStrings
using LsqFit
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
work_dir="/home/srinath/repo/Projects/cathode_packing_data/bimodal_distribution_database"
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
                    push!(df,[w,brr, bmr,b1,wt,wt1, wt2, ss])
                    cd("../")
                end #stat_dir
                if ss > 0
                    push!(df_avg, [w,brr, bmr,b1,wt_avg/ss,wt1_avg/ss, wt2_avg/ss, ss])
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
nmc_wt_ratio1 = collect(0.7:0.05:0.92)
lambda_min1 = [0.25,0.3,0.35,0.4,0.45, 0.5,0.55,0.6,
    0.65,0.75,0.8,0.9, 1.0, 1.1, 1.2, 1.25,
    1.4, 1.5, 1.6666667, 1.8, 1.9, 2.0, 2.5,3.0]
bimodal_radius_ratio2 = collect(0.2:0.005:0.4)
# bimodal_radius_ratio3 = collect(0.30:0.01:0.5)
# bimodal_radius_ratio1 = hcat(bimodal_radius_ratio2,bimodal_radius_ratio3)
bimodal_radisu_ratio1 = bimodal_radius_ratio2
bimodal_mix_ratio1 = collect(0.2:0.01:0.4)

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

si = pyimport("scipy.interpolate")
points=convert(Matrix,df_avg[:,[:wt_ratio, :bimodal_radius_ratio, :bimodal_mix_ratio, :particle_size_ratio]])
int_points=convert(Matrix,df_interp[:,[:wt_ratio, :bimodal_radius_ratio, :bimodal_mix_ratio, :particle_size_ratio]])
util = df_avg[!,:utilization]
pack = df_avg[!,:packing_density]
int_util = si.griddata(points, util, int_points, method="linear")
df_interp[!,:utilization] = int_util
int_pack = si.griddata(points, pack, int_points, method="linear")
df_interp[!,:packing_density] = int_pack
df_smooth= by(df_interp,[:wt_ratio,:bimodal_mix_ratio, :bimodal_radius_ratio],
    (:particle_size_ratio, :packing_density, :utilization)
    => x->(particle_size_ratio = x.particle_size_ratio,
           packing_density = x.packing_density,
           utlization = x.utilization,
           pack2 = smooth(x.particle_size_ratio, x.packing_density),
           util2 = smooth(x.particle_size_ratio, x.utilization)))


# CSV.write("all_data_smooth_interp.csv",df_smooth)
# dfs_wr = groupby(df_smooth, :wt_ratio)
# for (i,k) in enumerate(keys(dfs_wr))
#     fname = "all_data_wt_ratio_" * string(k.wt_ratio) * ".csv"
#     CSV.write(fname, dfs_wr[k])
# end
## Plotting

pyplot()

xlabel = L"$\lambda = \frac{D_{CAM}^{SM}}{D_{SE}}$"
ylabel = L"$\theta_{CAM} (\%)$"
ylabel2 = L"$\rho (\%)$"

p = plot(xlabel = ["" xlabel],
        ylabel = [ylabel ylabel2],
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
        grid = false)
dfx = @where(df_smooth, :wt_ratio .== 0.9)
gd_bmr = groupby(dfx, :bimodal_mix_ratio)
df_bmr = get(gd_bmr, (bimodal_mix_ratio=0.25,),nothing)
gd_brr = groupby(df_bmr,:bimodal_radius_ratio)
colors = [:blue, :orange, :green, :black, :red]
# q_brr = [0.21,0.25, 0.275, 0.285, 0.3, 0.4]
q_brr = [0.2, 0.25, 0.285, 0.3, 0.4]
j = 1
@. model(x,p) = 1.0/(1.0 + exp(-p[1]*(x-p[2])))
newxx = collect(0.1:0.01:3.0)
p0 = [1.0, 0.25]
for brr in keys(gd_brr)
    global j, colors
    if brr.bimodal_radius_ratio in q_brr
        xdata = gd_brr[brr][!,:particle_size_ratio]
        ydata = gd_brr[brr][!,:utlization]
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
            kstr = string(brr.bimodal_radius_ratio)
            labstr = latexstring(kstr)
            plot!(xdata, ydata .* 100.0, seriestype=:scatter, labels="", color = colors[j], marker = :x)
            plot!(newxx, newyy .* 100.0, label=labstr, color = colors[j])
            # println(kstr)
        end
        println(".............................")
        j += 1
    end
end
cd("julia")
plot!(collect(0.25:2.0:100), [97], seriestype=:hline,color = :black, label="")
png(p, "test_0.9.png")

cd("../")
# display(p)
###
