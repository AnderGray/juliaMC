####
#
#   Source class and generate function
#
#   Julia Version: V1.0
####
@with_kw struct Source

    E_distribution :: Sampleable = Truncated(Normal(14.2e7,6e6),1e6+1,3e8-1)   # The distribution of neutron energies, default set to normal
    position :: Array{Float64,1} = [0,0,0]                        # Position Where neutron is born, default is origin
    direction_Inclination :: Sampleable = Uniform(0,2*pi)          # Spherical coordinates for direction, the default is an
    direction_Azamuthal :: Sampleable = Uniform(0,1)             # isotropic distribution

end




# The generate function will generate a particle bank of length n from a source distribution s. This function could be altered later to also generate secondary neutrons to be added to an existing particel bank

@everywhere function generate(s :: Source, mat :: Material_Tendl, n :: Int64, grid :: Array{Float64,1})

    bank = Array{Particle,1}(UndefInitializer(),n)              # Initialization of particle bank array
    for i = 1:n
        bank[i]=Particle(id=i, mat = mat)           # Loop that fills the bank will default neutrons
    end

    # Delta dirac in energy
    Energys = 14.2e6*ones(n)#rand(s.E_distribution, m)             # Sampling of the Sources energy and direction distributions, position can
    d_Inc = rand(s.direction_Inclination,n)         # also be sampled here but will just be a point source for the time being
    d_Az = acos.(1 .- 2*rand(n))

    directions = zeros(n,3)

    directions[:,1] = sin.(d_Az) .* cos.(d_Inc)     # Conversion of the spherical coordinates that have been sampled to the
    directions[:,2] = sin.(d_Az) .* sin.(d_Inc)     # carteesian used by particle. Here we have done this by vectorization
    directions[:,3] = cos.(d_Az)                   # but is debatable if this is fast than looping in julia

    for i = 1:n
        bank[i].E = Energys[i]                      # Filling the particle bank with sampled information.
        bank[i].energyIndex = binarySearch(Energys[i], grid)
        bank[i].xyz = s.position
        bank[i].uvw = directions[i,:]
        #println(bank[i].uvw)
        #println(bank[i].energyIndex)
    end

    return bank
end
