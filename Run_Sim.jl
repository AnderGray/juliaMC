####
#   This is the main script for running juliaMC
#   Here, simulation variables are declared, for example the number or nuclides, Crossections and numbers of particles
#   The various simulations (TMC, FLY, Vanila) are also launched from here
#   The julia prerequisit packages can be installed by declaring "include("setup.jl")"
#   This script can be launched from the julia REPL by "include("Run_Sim.jl")"
#
#   To launch this script in parallel, in the REPL "using Distributed" followed by "addprocs(N)"
#   where N is your number of cores. Serial otherwise.
#
#       Julia version: V1.0
#
#           By: Ander Gray,         University of Liverpool, Culham Centre for Fusion Energy
####
using Distributed
using SharedArrays
@everywhere using LinearAlgebra             #@everywhere for distributed computing
include("juliaMC.jl");



# Function for Crossections, it is a linearly decreasing function with a gausian + cosine mixture
# A defines the mean of the gausian, shifts resonance peaks
en_grid = [i for i=1e6:2e5:3e8];
function xs(x,A)
    if x<14.e7
        y=(2e8*cos.(3.0*x).*pdf.(Normal(A,1e7),x)-5e-8x +3)+20;
    else
        y= 1000         #To increase the XS at the source energy
    end
end


## Two Crossections and 2 Nuclides. Crossections are properties of nuclides,
## and nuclides properties of the material class
scat1 = Cross_section(mt=1,reaction="elastic_scatter",energy_grid=en_grid, xs= xs.(en_grid,1e8));
absp1 = Cross_section(mt=2,reaction="absorption",energy_grid=en_grid, xs= xs.(en_grid,5e7));

scat2 = Cross_section(mt=1,reaction="elastic_scatter",energy_grid=en_grid, xs= xs.(en_grid,1.5e8));
absp2 = Cross_section(mt=2,reaction="absorption",energy_grid=en_grid, xs= xs.(en_grid,8e7));

nuclide1 = Nuclide(Name="Fe56",XS=[scat1,absp1]);
nuclide2 = Nuclide(Name="Fe54",XS=[scat2,absp2]);

material1 = Material(name="IronMixture", nucs=[nuclide1,nuclide2], atomic_density = [0.6,0.6], density = 16.0,id=1);

## Energy grid for tally class
tally_grid = [i for i=1e6:2e5:3e8];
tally1 = Flux_tally(energy_bins = tally_grid)

## material and tally are properties of the juliaMC class

simulation1 = juliaMC(material=material1,n=100000, grid = en_grid, Tally_batch=tally1, n_batch=10)
#@time runMovie(simulation1)    ##option to make gif of sim, warning: takes for ever
println("Vanila MC")
a = @time runPar(simulation1);

simulation1 = juliaMC(material=material1,n=100000, grid = en_grid, Tally_batch=tally1,n_batch=10)
println("TMC")
b = @time runTotalMonteCarlo(simulation1,100);
#=
simulation1 = juliaMC(material=material1,n=1000, grid = en_grid, Tally_batch=tally1,n_batch=10)
println("FlySampling")
c = @time runFlySampling(simulation1);

=#
#For creating plot
#plotTally(b,c)
plotTally(b,a,a)

##For saving in csv format. HDF5 preferible but not currently compatible with risk cluster
#en = (tally_grid[2:end]+tally_grid[1:end-1])/2;
#df = DataFrame(energy=en[1:end], XSTMC = b.Tally[1:end], stdTMC= b.std[1:end],XSFLY = c.Tally[1:end], stdFLY = c.std[1:end], XSVIN = a.Tally[1:end],stdVIN = a.std[1:end])
#CSV.write("simulationSave.csv",df);

#plotTally(c)
#=
#a=[1,5,10,50,100,500,1000,5000]

a = [10000, 50000, 100000, 500000, 1000000,5000000,10000000,50000000,100000000,500000000]
tal = zeros(length(a))
index=1
c=deepcopy(tally1)
for i in a

    simulation1 = juliaMC(material=material1,n=i, Tally_batch=tally1,n_batch=10)
    c=runFlySampling(simulation1)
    p = zeros(length(c.Tally))


    #println(c.std)
    #println(c.Tally)

    for j = 1:length(c.Tally)-28
        p[j] = c.std[j]/c.Tally[j]
    end
    tal[index]=norm(p)
    index+=1
end

#println(tal)

answer = plot(a,tal);

savefig(answer, "Answer3.png")

=#
