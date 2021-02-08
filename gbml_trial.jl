# Machine learning for elastic properties using Gradient Boost Machine learning
# Ref (Mark Asta and Kristen Persson paper )
## Modules
using DataFrames, DataFramesMeta, Statistics, LocallyWeightedRegression
using LossFunctions, CSV, JSON
using Plots, PlotUtils
using PyCall
using XGBoost
## Function to compute Holder means
function holder_mean(values::AbstractArray, power::Float64 = 0.0, weights = nothing, weights_norm = nothing)
    """
        holder_mean(values, power, weights, weights_norm)
    Calculate (possibly weighted) Holder (or power) mean
    :param: values: list or array of (typically elemental property) values
    :param: power: the power to which the Holder mean is to be calculated
    :param: weights: option list or array of weights
    :param: weights_norm: optional normalization method to be applied to weights
    :return: holder_mean
    """
    # println(power)
    alpha = nothing
    if isnothing(weights)
        alpha = 1 / length(values)
    else
        if length(values) != length(weights)
            throw("Length of weights has to be equal to length of values")
        else
            if weights_norm == "max" & max(weights) != 1.0
                weights .= weights ./ max(weights)
            elseif weights_norm == "sum" & sum(weights)
                weights .= weights ./ sum(weights)
            else
                throw("Either sum of weights or max weights has to be 1.0")
            end
            alpha = 1.0/sum(weights)
        end
    end
    if power == 0.0 # geometric mean
        if any(abs(value) < 1e-12 for value in values)
            return 0.0
        elseif any(value < 0 for value in values)
            return 0.0
        end

        if isnothing(weights)
            return exp(alpha * sum(log.(values)))
        else
            return exp(alpha * sum(weights .* log.(values)))
        end
    elseif power == 1.0 # arithmetic mean
        if isnothing(weights)
            return mean(values)
        else
            return mean(values, weights(weights))
        end

    end
    if any(value < 0  for value in values)
        if power % 2 == 0.0
            return 0.0
        end
    end

    if isnothing(weights)
        return (alpha * sum(values .^ power))^(1.0/power)
    else
        return (alpha * sum(weights .* values .^ power))^(1.0/power)
    end
end
##  Function to compute periodic_table data
composition = pyimport("pymatgen.core.composition")
element = pyimport("pymatgen.core.periodic_table")
function compute_atomic_fraction_weights(compound)
    """Compute the weights from atomic fractions using pymatgen interface
       return array of weights ... array size = size of compound"""

    # weights = DataFrame([Array{Float64,1}], [:weights])
    compo = composition.Composition(compound)
    w = []
    for (element_key,amount) in compo.get_el_amt_dict()
        elem = element.Element(element_key)
        push!(w, compo.get_atomic_fraction(elem))
    end
    # push!(weights, [w])
    # println(w)
    return w
end

# ## Data loading
# workdir = "/home/srinath/Projects/Machine_Learning"
# datadir = workdir * "/data"
# cd(datadir)
# df = DataFrame(CSV.File("all_mp_basic.csv"))
# composition = pyimport("pymatgen.core.composition")
# periodic_table = pyimport("pymatgen.core.periodic_table")
# df_elastic_hill = DataFrame()
# original_columns = names(df)
# for row in eachrow(df)
## adf
values = [1.0, 2.0, 3.0]
println(holder_mean(values, 0.0))
println(holder_mean(values, 1.0))
println(holder_mean(values, 2.0))
println(holder_mean(values, -1.0))
## Data loading
workdir = "/home/srinath/Projects/Machine_Learning"
datadir = workdir * "/data"
cd(datadir)
df1 = DataFrame(CSV.File("Elastic_mp.csv"))
df_sites = DataFrame(CSV.File("cell_nsite.csv"))
df = dropmissing!(innerjoin(df1, df_sites, on = :task_id), "pretty_formula")
original_columns = names(df)
## compute
actual_names = ["pretty_formula", "e_above_hull", "band_gap"
                , "density",
                "formation_energy_per_atom"]
# i = 1
transform!(df, ["elasticity"] => x-> replace.(x, "'"=>"\""))
transform!(df, ["volume","nsites"] => (x,y)-> log.(x./y))
df_elastic_hill = df[!,actual_names]
a = JSON.parse.(df[!,"elasticity_function"])
insertcols!(df_elastic_hill, "K_VRH" => get.(a, "K_VRH", nothing))
insertcols!(df_elastic_hill, "G_VRH" => get.(a, "G_VRH", nothing))
insertcols!(df_elastic_hill, "logV" => df[!,"volume_nsites_function"])

## --- Create subdataframe containing positive values of K, G
df_sub = @where(df_elastic_hill, :K_VRH .> 0.0, :G_VRH .> 0.0)
transform!(df_sub, "K_VRH" => ByRow(log) => "logK")
transform!(df_sub, "G_VRH" => ByRow(log) => "logG")
## --- Now add the compositional descriptors to the df_sub....
c_descriptors = ["atomic_mass", "atomic_number", "atomic_radius",
                "boiling_temperature", "electronegativity", "group_number",
                "melting_temperature", "row_number"]
holder_means = collect(-4.0:1.0:4.0)
new_col_names = []
for c in c_descriptors
    for h in holder_means
        push!(new_col_names, c * "_mu_" * string(h))
    end
end
df2 = DataFrame([String,Array{Float64,1}],[:pretty_formula, :weights])

for i in collect(1:nrow(df_sub))
    w = compute_atomic_fraction_weights(df_sub[i,"pretty_formula"])
    push!(df2, [df_sub[i,"pretty_formula"], w])
end
insertcols!(df_sub, "weights" => df2[!,:weights])
