####
#
#   Material Class, including two methods for passing interpolation request to nuclides
#
#   Julia Version: V1.0
####
@with_kw mutable struct Material

    # Properties of the Material
    # Manditory are: nuclides, atomic_density and density

    name :: String                                                         # Name of material
    nuclides :: Array{Any,1} = ["",:]                                      # Array of Nuclides names contained in the material
    nucs :: Array{Nuclide,1}                                               # Array of Nuclides contained in the material
    atomic_density :: Array{Float64,1}                                     # The Atomic densities of the nuclides in atom/b-cm²
    density :: Float64                                                     # The Total Density of the material in g/cm³

    id :: Int = 0                                                          # For multi material problems
    n_nuclides :: Int64  = 0                                               # Number of nuclides contained

    last_macro :: Float64 = 0.0                                            # Last interpolated MacroXS


    T_atomic :: Float64 = 0.0                                       # The Total Atomic denisty in atom/b-cm²
    weights :: Array{Float64,1} = [0.0,]                            # The relative density of the nuclides



    function Material(name :: String, nuclides :: Array{Any,1}, nucs :: Array{Nuclide,1}, atomic_density :: Array{Float64,1}, density :: Float64, id :: Int, n_nuclides :: Int64, last_macro :: Float64, T_atomic :: Float64, weights :: Array{Float64,1})

        T_atomic = sum(atomic_density)                                     # The Total Atomic denisty in atom/b-cm²
        weights = atomic_density/T_atomic

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

# Passing interpolation request nuclides class and storing result
function (obj :: Material)(E :: Float64)

    T_macro_xs=0

    for i = 1:obj.n_nuclides
        T_macro_xs += obj.weights[i].*obj.nucs[i](E)
    end

    T_macro_xs = T_macro_xs*obj.denisty

    obj.last_macro = T_macro_xs

    return T_macro_xs

end

function (obj :: Material)(indx :: Int64, E :: Float64)

    T_macro_xs=0

    for i = 1:obj.n_nuclides
        T_macro_xs += obj.weights[i].*obj.nucs[i](indx,E)
    end

    T_macro_xs = T_macro_xs*obj.denisty

    obj.last_macro = T_macro_xs

    return T_macro_xs

end

# Passing interpolation request nuclides class and storing result for ad hoc sampling
function (obj :: Material)(E :: Float64, Peturb)

    T_macro_xs=0

    for i = 1:obj.n_nuclides
        T_macro_xs += obj.weights[i].*obj.nucs[i](E, Peturb[i,:])
    end

    T_macro_xs = T_macro_xs*obj.denisty

    obj.last_macro = T_macro_xs

    return T_macro_xs

end

function (obj :: Material)(indx::Int64, E :: Float64, Peturb)

    T_macro_xs=0

    for i = 1:obj.n_nuclides
        T_macro_xs += obj.weights[i].*obj.nucs[i](indx, E, Peturb[i,:])
    end

    T_macro_xs = T_macro_xs*obj.denisty

    obj.last_macro = T_macro_xs

    return T_macro_xs

end


@with_kw mutable struct Material_Tendl

    # Properties of the Material
    # Manditory are: nuclides, atomic_density and density

    name :: String                                                         # Name of material
    nuclides :: Array{Any,1} = ["",:]                                      # Array of Nuclides names contained in the material
    nucs :: Array{Nuclide_Tendl,1}                                               # Array of Nuclides contained in the material
    atomic_density :: Array{Float64,1}                                     # The Atomic densities of the nuclides in atom/b-cm²
    density :: Float64                                                     # The Total Density of the material in g/cm³
    n_files :: Array{Int64,1} =[0,]

    id :: Int = 0                                                          # For multi material problems
    n_nuclides :: Int64  = 0                                               # Number of nuclides contained

    last_macro :: Float64 = 0.0                                            # Last interpolated MacroXS
    last_macroBounds :: Array{Float64,1} = [0.0,]

    T_atomic :: Float64 = 0.0                                       # The Total Atomic denisty in atom/b-cm²
    weights :: Array{Float64,1} = [0.0,]                            # The relative density of the nuclides



    function Material_Tendl(name :: String, nuclides :: Array{Any,1}, nucs :: Array{Nuclide_Tendl,1}, atomic_density :: Array{Float64,1}, density :: Float64, n_files :: Array{Int64,1}, id :: Int, n_nuclides :: Int64, last_macro :: Float64, last_macroBounds :: Array{Float64,1}, T_atomic :: Float64, weights :: Array{Float64,1})

        T_atomic = sum(atomic_density)                                     # The Total Atomic denisty in atom/b-cm²
        weights = atomic_density/T_atomic

        # The above must be changed

        n_nuclides = length(nucs)
        n_files = zeros(n_nuclides)
        for i = 1:n_nuclides
            nuclides[i] = nucs[i].Name
            n_files[i]=size(nucs[i].XS[1].xs)[1]
        end


        new(name, nuclides, nucs, atomic_density, density,n_files, id, n_nuclides, last_macro,last_macroBounds, T_atomic, weights)


    end

    #length(atomic_denisty) == n_nuclides || error("Number of Atomic Densities and Nuclides must be of equal length")
    #length(atomic_denisty) == n_nuclides ?:throw(ArguementError("asdadasd"))

    end

# Passing interpolation request nuclides class and storing result

# Passing interpolation request nuclides class and storing result for ad hoc sampling
function (obj :: Material_Tendl)(E :: Float64, Sample :: Array{Int64,1})

    T_macro_xs=0
    bounds = [0.0,0.0]

    for i = 1:obj.n_nuclides
        a = obj.weights[i].*obj.nucs[i](E, Sample[i])
        T_macro_xs += a[1]
        bounds +=a[2]
    end

    T_macro_xs = T_macro_xs * obj.density
    bounds = bounds*obj.density

    obj.last_macro = T_macro_xs
    obj.last_macroBounds = bounds

    return T_macro_xs, bounds

end

function (obj :: Material_Tendl)(indx::Int64, E :: Float64, Sample :: Array{Int64,1})

    T_macro_xs=0
    bounds = [0.0,0.0]

    for i = 1:obj.n_nuclides
        a = obj.weights[i].*obj.nucs[i](indx, E, Sample[i])
        T_macro_xs += a[1]
        bounds +=a[2]
    end

    T_macro_xs = T_macro_xs * obj.density
    bounds = bounds*obj.density

    obj.last_macro = T_macro_xs
    obj.last_macroBounds = bounds

    return T_macro_xs, bounds

end
