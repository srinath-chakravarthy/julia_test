using DataFrames, CSV, Glob
using Plots
# workdir = ("/home/srinath/data/Projects/moose_tests/supercomp_latest/curved_surface/dw_0.2_dh_0.05/no_initial_pressure/pr_2.0")
workdir = ("/home/srinath/github/electro_chemo_mech2/problems/benchmark4_fullCoupled/csv")
cd(workdir)
## read files
files=glob("full_model_current_*.csv", workdir)
files2=glob("full_model_interface_*.csv", workdir)
files3=glob("full_model_metal_flux_*.csv", workdir)
# files4 = glob("full_model_contact_pressure_*.csv", workdir)
post_proc_file = "full_model.csv"
dfs = DataFrame.(CSV.File.(files))
dfs2 = DataFrame.(CSV.File.(files2))
dfs3 = DataFrame.(CSV.File.(files3))
# dfs4 = DataFrame.(CSV.File.(files4))
df_post = DataFrame(CSV.File(post_proc_file))
## add index col for dfs
for i in 1:length(dfs)
    dfs[i][!,:time] .= df_post[i,:time]
    dfs2[i][!,:time] .= df_post[i,:time]
    dfs3[i][!,:time] .= df_post[i,:time]
    # dfs4[i][!,:time] .= df_post[i,:time]
    # dfs3[i][!,:flux] .= sqrt.(dfs3[i][!,:li_metal_flux_x].^2 + dfs3[i][!,:li_metal_flux_y].^2)
    # dfs[i][!,:flux] .= sqrt.(dfs[i][!,:li_ion_flux_x].^2 + dfs[i][!,:li_ion_flux_y].^2)
end
df = reduce(vcat, dfs)
df2 = reduce(vcat, dfs2)
df3 = reduce(vcat, dfs3)
# df4 = reduce(vcat, dfs4)

dfg = groupby(df, :time)
dfg2 = groupby(df2, :time)
dfg3 = groupby(df3, :time)
# dfg4 = groupby(df4, :time)
gr();
p = plot()
p2 = plot()
# p3 = plot()
for (i,t) in enumerate(keys(dfg))
    println(t)
    if (i == 12)
        # plot!(p, dfg3[i][!,:x], dfg3[i][!,:bndliflux], label="metal_flux")
        plot!(p2, dfg2[i][!,:x], dfg2[i][!,:mech_contact_normal_lm])
        plot!(p, dfg2[i][!,:x], abs.(dfg2[i][!,:thermal_lm]), label="lm")
        plot!(p, dfg[i][!,:x], dfg[i][!,:bndliflux], label="ion_flux")
    end
end
png(p, "test_ceramic.png")
png(p2, "test_metal.png")
# png(p3, "test_bound.png")
# p2 = plot()
# plot!(df_post)
