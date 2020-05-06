__precompile__()
module cathode_post_process
    using JSON
    using PyCall, LaTeXStrings
    using LightGraphs, MetaGraphs, SimpleWeightedGraphs, SparseArrays
    export post_process, load_file, find_paths!, find_utilization

    const cat = PyNULL()
    const ovitoio = PyNULL()
    const ovitomod = PyNULL()
    const ovitodata = PyNULL()
    const ovitopipe = PyNULL()
    const np = PyNULL()

    function __init__()
        copy!(cat, pyimport("packing.cathode"))
        copy!(ovitoio, pyimport("ovito.io"))
        copy!(ovitomod, pyimport("ovito.modifiers"))
        copy!(ovitodata, pyimport("ovito.data"))
        copy!(ovitopipe, pyimport("ovito.pipeline"))
        copy!(np, pyimport("numpy"))
    end
    ##
    # ---- Load JSON file for results
    function load_file(filename)
        catfile = nothing
        catfile = cat.Cathode.from_json(filename)
        return catfile
    end

    # function find_astar_paths!(G::AbstractGraph)
    #     dd = convert.(Float64,adjacency_matrix(G))
    #     for e in edges(G)
    #         s = LightGraphs.src(e)
    #         d = LightGraphs.dst(e)
    #         dd[s,d] = get_prop(G, e, :weight)
    #     end
    #     source = filter_vertices(G,:node_type, "source")
    #     source1 = collect(source)
    #     target = filter_vertices(G,:node_type, "target")
    #     for s in source
    #         for t in target
    #             xx = a_star(G, s, t, dd,
    #                 n->ifelse(get_prop(G,n,:ptype) in electrode_type, Inf, 0)
    #             if length(xx) > 0
    #                 set_prop!(G, s, :has_path, true)
    #                 break
    #         end
    #     end
    # end

    function find_paths!(G::AbstractGraph)
        dd = convert.(Float64,adjacency_matrix(G))
        for e in edges(G)
            s = LightGraphs.src(e)
            d = LightGraphs.dst(e)
            dd[s,d] = get_prop(G, e, :weight)
        end
        source = filter_vertices(G,:node_type, "source")
        source1 = collect(source)
        target = filter_vertices(G,:node_type, "target")

        for s in source
            set_prop!(G, s, :has_path, false)
        end
        all_se = collect(filter_vertices(G,:ptype, 3))
        for s in source
            for t in target
                fsource = filter(x->x!=s, source1)
                if has_path(G,s, t, exclude_vertices = fsource)
                    set_prop!(G,s, :has_path, true)
                    break
                end
            end
        end
    end

    function find_utilization(G::AbstractGraph)
        active_particles = filter_vertices(G,:active,true)
        hp = filter_vertices(G,:has_path, true)
        ac = filter_vertices(G,:active_connected, true)

        total_vol = 0.0
        active_vol = 0.0
        j_active = []
        for a in active_particles
            total_vol += get_prop(G,a, :volume)
            push!(j_active, a)
        end

        j_source_has_path = []
        for h in hp
            active_vol += get_prop(G, h, :volume)
            push!(j_source_has_path, h)
        end
        for a in ac
            active_vol += get_prop(G,a,:volume)
            push!(j_source_has_path,a)
        end

        num_util = length(j_source_has_path)/length(j_active)
        wt_util = active_vol/total_vol
        println(length(j_active))
        println(length(collect(ac)))
        println(length(collect(hp)))
        return num_util, wt_util
    end
    function create_graph(electrode_type::Array{Int64}, electrolyte_type::Array{Int64},filename="dump_final.cfg")
        poro = 0.0
        total_particle_volume = 0.0
        opipeline = ovitoio.import_file("dump_final.cfg")
        try
            voro = ovitomod.VoronoiAnalysisModifier(generate_bonds=true, use_radii=true)
        catch
            error("could not compute graph")
        end
        opipeline.modifiers.append(voro)
        opipeline.modifiers.append(
            ovitomod.ComputePropertyModifier(operate_on="bonds", output_property="Length", expressions=["BondLength"]))
        ##
        data = opipeline.compute()
        id = convert(Array{Int64},data.particles.identifiers)
        pos = convert(Matrix{Float64}, np.array(data.particles.position))
        radius = convert(Array{Float64},data.particles.get("Radius"))
        ptype = convert(Array{Int64},data.particles.particle_types)
        bond_topology = data.particles.bonds.topology
        bond_lengths = data.particles.bonds.get("Length")
        opipeline.compute()
        pcount = data.particles.count
        bonds_enum = ovitodata.BondsEnumerator(data.particles.bonds)
        ##Compute packing density
        zmax = maximum(pos[:,3])
        izmax = argmax(pos[:,3])
        rzmax = radius[izmax]
        set!(data.cell_, (2,2), zmax + rzmax)
        opipeline.compute()
        cell = data.cell
        cell_volume = data.cell.volume
        for r in radius
            total_particle_volume += (4.0/3.0) * pi * r^3
        end
        poro = total_particle_volume/cell_volume
        println(cell_volume)
        println(total_particle_volume)
        println(poro)
        G = MetaGraph(pcount)
        # give ids to all particles in the graph and the particle type
        for pid in range(1,stop=pcount)
            set_props!(G,pid,Dict(:id=>id[pid], :ptype=>ptype[pid],
                                  :node_type=>nothing, :has_path=>false,
                                  :volume=>radius[pid]^3, :active=>false,
                                  :active_connected=>false))
            if ptype[pid] in electrode_type
                set_prop!(G,pid, :active, true)
                if pos[pid,3] - radius[pid] < 1e-3
                    set_prop!(G,pid, :active_connected, true)
                end
            end
        end
        # --- Add graph edges and weights
        for pid in range(1,stop = pcount)
            for bid in bonds_enum.bonds_of_particle(pid-1)
                a = convert(Int64, get(bond_topology,(bid, 0))) + 1
                b = convert(Int64, get(bond_topology,(bid, 1))) + 1
                # blength = bond_lengths[bid + 1]
                blength = get(bond_lengths,bid,)
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
                        if pos[a,3]- radius[a] < 1e-3
                            set_prop!(G,a, :node_type, "target")
                        end
                    end
                    if ptype[b] in electrolyte_type
                        if pos[b,3]- radius[b] < 1e-3
                            set_prop!(G,b, :node_type, "target")
                        end
                    end
                end
            end
        end
        return G, poro
    end


    function post_process(;directory::String="./",
                          filename::AbstractString="dump_final.cfg",
                          electrode_type::Array{Int64,1}=[1],
                          electrolyte_type::Array{Int64,1}=[1],
                          compute_paths::Bool=true,
                          assume_particle_type = false)
        poro = 0.0
        num_util = 0.0
        wt_util = 0.0
        catfile = load_file("cathode.json")
        if !assume_particle_type
            electrode_type = Int64[]
            electrolyte_type = Int64[]
            for a in catfile.cmateriallist
                if a.active
                    if a.active_mat
                        push!(electrode_type, a.sim_particle_type)
                    else
                        push!(electrolyte_type, a.sim_particle_type)
                    end
                end
            end
        end
        try
            G, poro = create_graph(electrode_type, electrolyte_type)
        catch
            println("No graph could be computed")
            return poro, num_util, wt_util
        end
        find_paths!(G)
        num_util, wt_util = find_utilization(G)

        return poro, num_util, wt_util
    end
end
