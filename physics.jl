####
#
#   Conatains transport() and transportUQ() functions for transport of particles
#   using vanila monte-carlo and ad hoc sampling
#
#   Also contained are functions for selecting nulicides and reactions, and an
#   elastic scattering function
#
#   Julia Version: V1.0
####
function transport(p :: Particle)

    E = p.E
    ind = p.energyIndex
    #=
    println("before")
    println(ind)
    println(E)
    =#
    #Total_Macro = p.mat(E)                      # interpolation of XS and construction on MacroXS
    #Total_Macro = p.mat(ind,E)
    Total_Macro = p.mat(ind,E)
    d = -log(rand())/Total_Macro                # distance to next collision

    # pre-collision details stored
    p.last_xyz = p.xyz
    p.last_E = E
    p.last_wgt = p.wgt
    p.last_uvw = p.uvw
    p.last_d = d
    p.last_index = ind

    # new particle position. p.uvw always of unit length
    p.xyz=p.uvw*d


    #println("Particle has traveld $d and is now at $(p.xyz)")

    s = select(p.mat);                       # first select nuclide
    r = select(p.mat.nucs[s]);               # then reaction. Select() has been overloaded.
    react = p.mat.nucs[s].XS[r].reaction

    if react == "elastic_scatter"
        p = scatter(p)
        #println("scatter at $(p.xyz)")
#=
        if p.E < 1e6
            p.alive = false
            p.wgt = 0;
            #p.last_E = 1.42e8;          # cannot remember why I did this...
            #println("Particle $(p.id) has lost too much energy at $(p.xyz)")
        end
        =#

    elseif react == "absorption"
        p.alive = false
        p.wgt=0
        p.last_reaction="absorption"
        #println("absorption at $(p.xyz)")
    else
        error("Something has gone terribly wrong")
    end

    return p

end

function transport(p :: Particle, Sample :: Array{Int64,1})

    E = p.E
    ind = p.energyIndex
    #=
    println("before")
    println(ind)
    println(E)
    =#
    #Total_Macro = p.mat(E)                      # interpolation of XS and construction on MacroXS
    #Total_Macro = p.mat(ind,E)
    temp = p.mat(ind, E, Sample)
    #temp = p.mat(E, Sample)
    Total_Macro = temp[1]
    Total_bounds = temp[2]
    ran = rand()
    d = -log(ran)/Total_Macro
    d_bounds = -log(ran)./Total_bounds
    d_bounds = sort(d_bounds)

    # pre-collision details stored
    p.last_xyz = p.xyz
    p.last_E = E
    p.last_wgt = p.wgt
    p.last_uvw = p.uvw
    p.last_d = d
    p.last_d_bounds = d_bounds
    p.last_index = ind

    # new particle position. p.uvw always of unit length
    p.xyz=p.uvw*d


    #println("Particle has traveld $d and is now at $(p.xyz)")

    s = select(p.mat);                       # first select nuclide
    r = select(p.mat.nucs[s]);               # then reaction. Select() has been overloaded.
    p.last_Atomic_Weight = p.mat.nucs[s].atomicWeight;
    react = p.mat.nucs[s].XS[r].reaction;

    if react == "elastic_scatter"
        p = scatterOpenMC(p)
        #println("scatter at $(p.xyz)")
#=
        if p.E < 1e6
            p.alive = false
            p.wgt = 0;
            #p.last_E = 1.42e8;          # cannot remember why I did this...
            #println("Particle $(p.id) has lost too much energy at $(p.xyz)")
        end
        =#

    elseif react == "absorption"
        p.alive = false
        p.wgt=0
        p.last_reaction="absorption"
        #println("absorption at $(p.xyz)")
    else
        error("Something has gone terribly wrong")
    end

    return p

end

function transportUQ(p :: Particle, Sample)

    E = p.E
    ind = p.energyIndex
    # only difference to transport(). Perturbation passed right down to XS interpolator
    temp = p.mat(ind, E, Sample)
    Total_Macro = temp[1]
    Total_bounds = temp[2]
    ran = rand()
    d = -log(ran)/Total_Macro
    d_bounds = -log(ran)./Total_bounds
    d_bounds = sort(d_bounds)

    p.last_xyz = p.xyz
    p.last_E = E
    p.last_index = ind
    p.last_wgt = p.wgt
    p.last_uvw = p.uvw
    p.last_d = d
    p.last_d_bounds = d_bounds

    p.xyz=p.uvw*d


    #println("Particle has traveld $d and is now at $(p.xyz)")

    s = select(p.mat);
    r = select(p.mat.nucs[s]);
    p.last_Atomic_Weight = p.mat.nucs[s].atomicWeight;

    react = p.mat.nucs[s].XS[r].reaction

    if react == "elastic_scatter"
        p = scatter(p)
        #println("scatter at $(p.xyz)")
        #=
        if p.E < 1e6
            p.alive = false
            p.wgt = 0;
            p.last_E = 1.42e8;
            #println("Particle $(p.id) has lost too much energy at $(p.xyz)")
        end
        =#
    elseif react == "absorption"
        p.alive = false
        p.wgt=0
        p.last_reaction="absorption"

        #println("absorption at $(p.xyz)")

    else
        error("Something has gone terribly wrong")
    end

    return p

end


# Function for selecting a nuclide in the material for reaction
# NOTE: energy does not need to be passed here. XS only interpolated once and
# stored in last_T_xs_value
function select(mat :: Material_Tendl)

    n=length(mat.nucs)
    cdf_values = zeros(n+1)


    ##
    #   Problem here! Change descrete cdf class!
    ##
    # cdf values taken to be the macroscopic XS of each nuclide
    for i =1:n
        cdf_values[i+1] = mat.nucs[i].last_T_xs_value.*mat.weights[i].+cdf_values[i]
    end
    #cdf_values=filter(e->e!=0.0,cdf_values)
    println(cdf_values)
    a=sum(cdf_values)
    cdf_values=cdf_values./a
    println(cdf_values)
    # Descrete cdf created and sampled
    #println(cdf_values)
    a = Array{String, 1}(UndefInitializer(),n)          # intitalization of a array of strings
    for i =1:n
        a[i]=mat.nuclides[i]
    end
    dist = Descrete_CDF(a,cdf_values)
    nad, selection = dist()

    #println("Selected nuclide is: $(nad.Name) and $selection")

    return selection
end

# selection of XS in nuclide
function select(nuc :: Nuclide_Tendl)

    n=length(nuc.XS)
    cdf_values = zeros(n+1)

    #cdf constructed from microscopic XS of reactions
    for i =1:n
        cdf_values[i+1] = nuc.XS[i].last_xs_value+cdf_values[i]
    end
    #cdf_values=filter(e->e!=0.0,cdf_values)
    cdf_values=cdf_values/nuc.last_T_xs_value
    a = Array{String, 1}(UndefInitializer(),n)          # intitalization of a array of strings
    for i =1:n
        a[i]=nuc.XS[i].reaction
    end

    dist = Descrete_CDF(a,cdf_values)
    nad, selection = dist()

    #println("Selected reactions is $nad and $selection")

    return selection

end


function scatter(p :: Particle)

    direction_Inclination = Uniform(0,2*pi)          # Spherical coordinates for direction, the default is an
    #direction_Azamuthal = Uniform(0,pi)             # isotropic distribution. THIS COULD BE BIASED

    d_Inc = rand(direction_Inclination)           # also be sampled here but will just be a point source for the time being
    d_Az = acos(1-2*rand())

    directions = zeros(3)

    directions[1] = sin(d_Az) * cos(d_Inc)     # Conversion of the spherical coordinates that have been sampled to the
    directions[2] = sin(d_Az) * sin(d_Inc)     # carteesian used by particle. Here we have done this by vectorization
    directions[3] = cos(d_Az)                   # but is debatable if this is fast than looping in julia

    #println("The sampled direction is $directions")

    p.E = p.E - rand()*p.E;
    #p.E = p.E - rand(Uniform(0,0.2))*p.E;     # reduce energy by a random amount
    p.uvw = directions[:];
    p.last_reaction="scatter"

    return p

end

function scatterOpenMC(p :: Particle)

    E = p.E;                    #Lab Frame energy
    awr = p.last_Atomic_Weight;       #Atomic wieght ratio of traget
    uvw = p.uvw;                #Lab direction
    vel = sqrt(E);              #Lab Speed

    v_n = vel*uvw               #Lab frame velocity
    v_t = [0,0,0];                    #Target velocity

    #centre of mass velocity
    v_cm = (v_n + awr .* v_t)./(awr + 1.0);


    v_n = v_n - v_cm;   #Neutron vel in CoM
    vel = norm(v_n);    #Neutron speed in CoM
    mu_cm = 2.0 *rand() - 1.0;   #isotropic
    u_cm = v_n/vel;     #Direction in CoM

    phi = 2.0 *pi*rand();

    ub = rotate_angle(u_cm, mu_cm, phi);
    v_n =  vel * ub;

    v_n = v_n + v_cm;

    E = dot(v_n,v_n);
    vel = sqrt(E);

    u = v_n/vel;

    p.E = E;
    p.uvw = u;

    p.last_reaction="scatter"

    return p

end
