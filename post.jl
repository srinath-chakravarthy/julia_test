using cathode_post_process
using JSON

if !isfile("results.json")
    poro, num_util, wt_util = post_process(assume_particle_type = true)
    println("$poro $num_util $wt_util")
    Results = Dict()
    Results["Results"] = Dict("Particle_fraction" => num_util, "Active_weight_fraction" => wt_util, "Output_porosity" => poro)
    open("results.json","w") do f
        JSON.print(f, Results, 4)
    end
end
