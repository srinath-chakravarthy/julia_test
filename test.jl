using CSV, DataFrames, Statistics
workdir = "/home/srinath/github/electro_chemo_mech2/problems/mcmeeking_paper_results/data"
cd(workdir)
filebase = "carbon_black_equil_potential.csv"
filebase2 = "carbon_black_equil_potential_mV.csv"
df = DataFrame(CSV.File(filebase))
df[!,2] .*= 1000
CSV.write(filebase2, df)
