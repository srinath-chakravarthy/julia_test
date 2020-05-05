using DataFrames, DataFramesMeta, CSV, Statistics
using Colors, ColorSchemes, ColorSchemeTools
using JSON
using Plots
using PyCall, Conda, LaTeXStrings
using LsqFit
using LightGraphs, MetaGraphs, SimpleWeightedGraphs

cat = pyimport("packing.cathode")
ovitoio = pyimport("ovito.io")
ovitomod = pyimport("ovito.modifiers")
ovitodata = pyimport("ovito.data")
ovitopipe = pyimport("ovito.pipeline")
cd("/home/srinath/Projects/cathode_packing_data/ceder_data/wt_ratio_0.8/a_5.0/ratio_4.0")
try
    catfile = cat.Cathode.from_json("cathode.json")
catch
    error("file not found")
end
electrode_type = [1]
electrolyte_type = [2]
pipeline = ovitoio.import_file("dump_final.cfg")
voro = ovitomod.VoronoiAnalysisModifier(generate_bonds=true, use_radii=true)
pipeline.modifiers.append(voro)
pipeline.modifiers.append(
    ovitomod.ComputePropertyModifier(operate_on="bonds", output_property="Length", expressions=["BondLength"]))
##
data = pipeline.compute()
id = convert(Array{Int64},data.particles.identifiers)
pos = data.particles.position
radius = convert(Array{Float64},data.particles.get("Radius"))
ptype = convert(Array{Int64},data.particles.particle_types)
bond_topology = data.particles.bonds.topology
bond_lengths = data.particles.bonds.get("Length")
pipeline.compute()
pcount = data.particles.count
bonds_enum = ovitodata.BondsEnumerator(data.particles.bonds)
G = MetaGraph(pcount)
# give ids to all particles in the graph and the particle type
for pid in range(1,stop=pcount)
    set_props!(G,pid,Dict(:id=>id[pid], :ptype=>ptype[pid], :node_type=>nothing, :has_path=>false, :volume=>radius[pid]^3, :active=>false))
    if ptype[pid] in electrode_type
        set_prop!(G,pid, :active, true)
        if pos[pid,3] - radius[pid] < 1e-3
            set_prop!(G,pid, :has_path, true)
        end
    end
end

# --- Add graph edges and weights
for pid in range(1,stop = pcount)
    for bid in bonds_enum.bonds_of_particle(pid-1)
        a = convert(Int64, get(bond_topology,(bid, 0))) + 1
        b = convert(Int64, get(bond_topology,(bid, 1))) + 1
        blength = bond_lengths[bid + 1]
        # if pid == 1
        #     println("$bid $a $b $blength")
        # end
        if blength < (radius[a] + radius[b] + 1.0e-3)
            if a == pid
                if ptype[a] in electrode_type && ptype[b] in electrode_type
                    nothing
                elseif ptype[a] in electrolyte_type && ptype[b] in electrolyte_type
                    add_edge!(G, a, b)
                    set_prop!(G, a, b, :weight, 0.1)
                elseif ptype[a] in electrode_type && ptype[b] in electrolyte_type
                    add_edge!(G, a, b)
                    set_prop!(G, a, b, :weight, 1.0)
                else
                    add_edge!(G,a , b)
                    set_prop!(G, a, b, :weight, 50.0)
                end
                if ptype[a] in electrode_type
                    # println("Source $a")
                    if pos[a,3] > radius[a] > 1e-3
                        set_prop!(G,a, :node_type, "source")
                    end
                end
            else
                if ptype[a] in electrode_type && ptype[b] in electrode_type
                    nothing
                elseif ptype[a] in electrolyte_type && ptype[b] in electrolyte_type
                    add_edge!(G, b, a)
                    set_prop!(G, b, a, :weight, 0.1)
                elseif ptype[b] in electrode_type && ptype[a] in electrolyte_type
                    add_edge!(G, b, a)
                    set_prop!(G, b, a, :weight, 1.0)
                else
                    add_edge!(G, b, a)
                    set_prop!(G, b, a, :weight, 50.0)
                end
                if ptype[b] in electrode_type
                    if pos[b,3] - radius[b] > 1e-3
                        # println("Source $b")
                        set_prop!(G,b, :node_type, "source")
                    end
                end
            end
            if ptype[a] in electrolyte_type
                if pos[a,3] - radius[a] < 1e-3
                    set_prop!(G,a, :node_type, "target")
                end
            end
            if ptype[b] in electrolyte_type
                if pos[b,3] - radius[b] < 1e-3
                    set_prop!(G,b, :node_type, "target")
                end
            end
        end
    end
end
## Assign distances
dd = convert.(Float64,adjacency_matrix(G))
for e in edges(G)
    s = LightGraphs.src(e)
    d = LightGraphs.dst(e)
    dd[s,d] = get_prop(G, e, :weight)
end
source = filter_vertices(G,:node_type, "source")
target = filter_vertices(G,:node_type, "target")
for s in source
    for t in target
        xx = a_star(G,s,t,dd)
        if length(xx) > 2
            set_prop!(G,s, :has_path, true)
            # println("Path found between $s and $t")
            break
        end
    end
end
##
active_particles = filter_vertices(G,:active,true)
hp = filter_vertices(G,:has_path, true)

total_vol = 0.0
active_vol = 0.0
j_active = []
for a in active_particles
    global total_vol
    total_vol += get_prop(G,a, :volume)
    push!(j_active, a)
end

j_source_has_path = []
for h in hp
    global active_vol
    active_vol += get_prop(G, h, :volume)
    push!(j_source_has_path, h)
end
