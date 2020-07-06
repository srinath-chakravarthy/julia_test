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

function l2(x,y)
    @. model(x, p) = p[1] * x[:,1] * p[2] * x[:,2] + p[3] * x[:,3] +
                        p[4] * x[:,2]^2 + p[5] * x[:,3]^2 + p[6] * x[:,2] * x[:,3] +
                        p[7] * x[:,1] * x[:,2] + p[8] * x[:,1] * x[:,3]
    # println(x)
    p0 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0,1.0,1.0]
    fit = curve_fit(model, x, y, p0)
    return model, fit
end


function l3(x,y)
    @. model(x, p) = p[1] * x[:,1]  + p[2] * x[:,2] +
                        p[3] * x[:,1]^2 + p[4] * x[:,2]^2 + p[5] * x[:,1] * x[:,2] +
                        p[6]
    # println(x)
    p0 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    fit = curve_fit(model, x, y, p0)
    return model, fit
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
plot_dir="/home/srinath/Projects/cathode_packing_data/bimodal_distribution_database/patent_plots/"
cd(work_dir)
cd("julia")
nmc_wt_ratio1 = collect(0.7:0.01:0.9)
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
df_interp = DataFrame(CSV.File("all_data_interp_julia_util_packing.csv"))
cd(work_dir)
cd("julia")
df_avg = DataFrame(CSV.File("orig_data_julia_util_packing.csv"))
cd(work_dir)
##

si = pyimport("scipy.interpolate")

df_smooth= by(df_interp,[:wt_ratio,:bimodal_mix_ratio, :bimodal_radius_ratio],
    (:particle_size_ratio, :packing_density, :utilization)
    => x->(particle_size_ratio = x.particle_size_ratio,
           packing_density = x.packing_density,
           utlization = x.utilization,
           pack2 = smooth(x.particle_size_ratio, x.packing_density),
           util2 = smooth(x.particle_size_ratio, x.utilization)))
# df_smooth= by(df_interp,[:wt_ratio,:bimodal_mix_ratio, :bimodal_radius_ratio],
#    (:particle_size_ratio, :packing_density, :utilization)
#     => x->(util2 = smooth_util(x.particle_size_ratio, x.utilization),
#            pack2 = smooth_pack(x.particle_size_ratio, x.packing_density)))
# display(p)
###
#minimum particle_size_Ratio for f_cam > 85% and utilizaton > 99%
df2 = @where(df_avg, :utilization .> 0.995, :wt_ratio .>= 0.85)
df2[!,:lambda] = df2[!,:particle_size_ratio] ./ df2[!,:bimodal_radius_ratio]

dfx= by(df2,[:bimodal_mix_ratio, :bimodal_radius_ratio],
          (:particle_size_ratio, :utilization)
          => x->(lambda = minimum(x.particle_size_ratio)))

dfx[!,:lambda] = dfx[!,3] ./ dfx[!,:bimodal_radius_ratio]
x = dfx[!,:bimodal_mix_ratio]
y = dfx[!,:bimodal_radius_ratio]
z = dfx[!,:lambda]
model, fit = l3(Matrix(dfx[!,1:2]), dfx[!,:lambda])
dfx[!,:lambda2] = model(dfx[!,1:2], coef(fit))
z1 = dfx[!,:lambda2]
cd(plot_dir)
CSV.write("orig_data_f_cam_gt85_utilmax_lambda_min.csv",df_smooth)
cd(work_dir)

### Maximum particle_size_Ratio for f_cam > 85% and utilizaton > 99%
df2 = @where(df_avg, :utilization .> 0.995, :wt_ratio .>= 0.85)
df2[!,:lambda] = df2[!,:particle_size_ratio] ./ df2[!,:bimodal_radius_ratio]

dfx= by(df2,[:bimodal_mix_ratio, :bimodal_radius_ratio],
          (:particle_size_ratio, :utilization)
          => x->(lambda = maximum(x.particle_size_ratio)))

dfx[!,:lambda] = dfx[!,3] ./ dfx[!,:bimodal_radius_ratio]
x = dfx[!,:bimodal_mix_ratio]
y = dfx[!,:bimodal_radius_ratio]
z = dfx[!,:lambda]
model, fit = l3(Matrix(dfx[!,1:2]), dfx[!,:lambda])
dfx[!,:lambda2] = model(dfx[!,1:2], coef(fit))
z1 = dfx[!,:lambda2]
cd(plot_dir)
CSV.write("orig_data_f_cam_gt85_utilmax_lambda_max.csv",df_smooth)
cd(work_dir)


##
df2 = @where(df_smooth, :bimodal_mix_ratio .== 0.25, :bimodal_radius_ratio .== 0.285)
Dcambig = 14.0
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
