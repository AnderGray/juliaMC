####
#
#   Crossection and Nuclide classes including interpolators
#
#   Julia Version: V1.0
####
@with_kw mutable struct Cross_section
    mt :: Int32
    reaction :: String
    energy_grid :: Array{Float64,1}
    xs :: Array{Float64,1}
    last_xs_value :: Float64 = 0
end


# Linear interpolation of XS
function (obj::Cross_section)(E :: Float64)

    inx = findInter(E, obj.energy_grid)

    y0=obj.xs[inx];
    y1=obj.xs[inx+1];
    x0=obj.energy_grid[inx];
    x1=obj.energy_grid[inx+1];

    XS_e = y0+(E - x0)*(y1-y0)/(x1-x0);

    obj.last_xs_value = XS_e

    return XS_e

end

function (obj::Cross_section)(E :: Float64, Peturb :: Float64)

    inx = findInter(E, obj.energy_grid)

    y0=obj.xs[inx];
    y1=obj.xs[inx+1];
    x0=obj.energy_grid[inx];
    x1=obj.energy_grid[inx+1];

    XS_e = y0+(E - x0)*(y1-y0)/(x1-x0);

    XS_e = XS_e*Peturb          # purtubation of xs after interpolation

    obj.last_xs_value = XS_e

    return XS_e

end



@with_kw mutable struct Nuclide

    Name :: String
    XS :: Array{Cross_section,1}
    last_T_xs_value :: Float64 = 0           # The last calculated XS value, so we don't have to interpolate again

end

function (obj :: Nuclide)(E :: Float64)

    T_micro_xs = 0
    for i = 1:length(obj.XS)
        T_micro_xs += obj.XS[i](E)
    end
    obj.last_T_xs_value = T_micro_xs
    return T_micro_xs
end


function (obj :: Nuclide)(E :: Float64, Peturb)

    T_micro_xs = 0
    for i = 1:length(obj.XS)
        T_micro_xs += obj.XS[i](E, Peturb[i])
    end
    obj.last_T_xs_value = T_micro_xs
    return T_micro_xs
end




# NOT USING THE BELOW

# This is the constructor for a summed crossection class, this can include the total microscopic XS for nuclides and total macroscopic XS for materials

@with_kw mutable struct Summed_xs

    mt :: Int32
    energy_grid :: Array{Float64,1}
    xs :: Array{Float64,1} = ones(2,1)

    comp :: Array{Cross_section,1}
    wgt :: Array{Float64,1} = ones(length(comp,1))


    function Summed_xs(mt :: Int32, energy_grid , xs, comp, wgt)

        if !check_grids(obj)
            for i = 1:length(comp)
                comp[i] = reshapeXS(comp[i],energy_grid)
            end
        end

        xs = sum(wgt.*comp.xs) # add wieghting here

        new(mt, energy_grid, xs, comp, wgt)
    end

end

function (obj::Summed_xs)(E :: Float64)

    inx = findInter(E, obj.energy_grid)

    y0=obj.xs[inx];
    y1=obj.xs[inx+1];
    x0=obj.energy_grid[inx];
    x1=obj.energy_grid[inx+1];

    XS_e = y0+(E - x0)*(y1-y0)/(x1-x0);

    return XS_e

end


function check_grids(SumXs :: Summed_xs)

    test = true

    for i =2:1:length(SumXs.composites)
        if SumXs.energy_grid != SumXs.compostes[i].energy_grid
            test=false
        end
    end

    return test
end

function reshapeXS(xs :: Cross_section, energy :: Array{Float64, 1})

    newXs :: Array{Float,1}(length(xs.energy_grid))

    for i=1:length(energy)
        newXs[i] = xs(energy[i])
    end
    xs.xs=newXs
    xs.energy_grid=energy;

    return xs


end
