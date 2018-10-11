@with_kw mutable struct Material

    # Properties of the Material
    # Manditory are: nuclides, atomic_density and density

    name :: String                                                         # Name of material
    nuclides :: Array{Any,1} = ["",:]                                    # Array of Nuclides names contained in the material
    nucs :: Array{Nuclide,1}                                               # Array of Nuclides contained in the material
    atomic_density :: Array{Float64,1}                                     # The Atomic densities of the nuclides in atom/b-cm²
    density :: Float64                                                     # The Total Density of the material in g/cm³

    id :: Int = 0                                                          # For multi material problems
    n_nuclides :: Int64  = 0                                               # Number of nuclides contained

    last_macro :: Float64 = 0.0


    T_atomic :: Float64 = 0.0                                       # The Total Atomic denisty in atom/b-cm²
    weights :: Array{Float64,1} = [0.0,]                            # The relative density of the nuclides



    function Material(name :: String, nuclides :: Array{Any,1}, nucs :: Array{Nuclide,1}, atomic_density :: Array{Float64,1}, density :: Float64, id :: Int, n_nuclides :: Int64, last_macro :: Float64, T_atomic :: Float64, weights :: Array{Float64,1})

        T_atomic = sum(atomic_density)                                     # The Total Atomic denisty in atom/b-cm²
        weights = density*atomic_density/T_atomic

        # The above must be changed

        n_nuclides = length(nucs)

        for i = 1:n_nuclides
            nuclides[i] = nucs[i].Name
        end


        new(name, nuclides, nucs, atomic_density, density, id, n_nuclides, last_macro, T_atomic, weights)


    end

    #length(atomic_denisty) == n_nuclides || error("Number of Atomic Densities and Nuclides must be of equal length")
    #length(atomic_denisty) == n_nuclides ?:throw(ArguementError("asdadasd"))

    end

function (obj :: Material)(E :: Float64)

    T_macro_xs=0

    for i = 1:obj.n_nuclides
        T_macro_xs += obj.weights[i].*obj.nucs[i](E)
    end

    obj.last_macro = T_macro_xs

    return T_macro_xs

end

function (obj :: Material)(E :: Float64, Peturb)

    T_macro_xs=0

    for i = 1:obj.n_nuclides
        T_macro_xs += obj.weights[i].*obj.nucs[i](E, Peturb[i,:])
    end

    obj.last_macro = T_macro_xs

    return T_macro_xs

end
