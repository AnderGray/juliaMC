@with_kw mutable struct Particle
    
    #Particle information
    id :: Int32 = 0                             # So the particle can be tracked throughout the simulation
    Type :: Char = 'n'                          # i.e (n, p, e, ...)
    mat :: Material                             # The material the particle currenty is in
    
    #Particle Coordinates
    xyz :: Array{Float64,1}=[0,0,0]             # Carteesian Position
    uvw :: Array{Float64,1}=[0,0,0]             # Carteesian Direction
     
    #Physical Properties
    E :: Float64 = 0                            # Particles energy in eV
    wgt :: Float64 = 1                          # Weight of the particle for non-analogue 
    alive :: Bool = true                        # The state of the particles existance
    
    #precollision Properties
    last_xyz :: Array{Float64,1} = [0,0,0]
    last_uvw :: Array{Float64,1} = [0,0,0]
    last_E :: Float64 = 0
    last_wgt :: Float64 = 0
    last_reaction :: String = ""
    last_d :: Float64 = 0
    
    
    
    
    end