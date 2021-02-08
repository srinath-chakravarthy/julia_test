using DataFrames, CSV, Glob
using Plots
workdir = ("/home/srinath/data/Projects/moose_tests/supercomp_latest/curved_surface/dw_0.2_dh_0.05/rst")
cd(workdir)
## read files
files=glob("curved_restart3_bound_flux_ceramic*.csv", workdir)
files2=glob("curved_restart3_bound_flux_interlayer*.csv", workdir)
post_proc_file = "curved_restart3.csv"
dfs = DataFrame.(CSV.File.(files))
df_post = DataFrame(CSV.File(post_proc_file))
## add index col for dfs
for i in 1:length(dfs)
    dfs[i][!,:time] .= df_post[i,:time]
end
df = reduce(vcat, dfs)

dfg = groupby(df, :time)
