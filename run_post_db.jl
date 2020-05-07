using JSON
using DataFrames, Plots
wt_ratio = collect(0.6:0.05:0.95)
lambda = [0.5, 0.75, 1.0, 1.25, 1.5, 1.667, 2.0, 2.5, 4.0, 5.0, 8.0]
workdir = "/home/srinath/Projects/cathode_packing_data/new_ceder_data_rerun/"
stat = collect(1:5)
cd(workdir)
##
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
                wt_util = res["Results"]["Active_wt_fraction"]
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

df_smooth= by(df_interp,[:wt_ratio],
    (:particle_size_ratio, :packing_density, :utilization)
    => x->(particle_size_ratio = x.particle_size_ratio,
           packing_density = x.packing_density,
           utilization = x.utilization,
           pack2 = smooth(x.particle_size_ratio, x.packing_density),
           util2 = smooth(x.particle_size_ratio, x.utilization)))
