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
using HDF5
@everywhere using LinearAlgebra             #@everywhere for distributed computing
include("juliaMC.jl");



# Function for Crossections, it is a linearly decreasing function with a gausian + cosine mixture
# A defines the mean of the gausian, shifts resonance peaks



## Two Crossections and 2 Nuclides. Crossections are properties of nuclides,
## and nuclides properties of the material class

numFe = 613;
numO = 641;
XSDIR = "../tendlHDF5"


energy = h5read("$XSDIR/Fe056/Fe056_0000.h5","energy");

scatFe = zeros(numFe, length(energy));
abspFe = zeros(numFe, length(energy));
totalFe = zeros(numFe, length(energy));

scatO = zeros(numO, length(energy));
abspO = zeros(numO, length(energy));
totalO = zeros(numO, length(energy));

for i =1:numFe
    if i <10
        nums = "000$i";
    elseif i<100
        nums = "00$i";
    else
        nums = "0$i";
    end
    files = "Fe056_$nums.h5";
    scatFe[i,:] = h5read("$XSDIR/Fe056/$files","elastic");
    abspFe[i,:] = h5read("$XSDIR/Fe056/$files","absorption");
    totalFe[i,:] = h5read("$XSDIR/Fe056/$files","total");
end

for i =1:numO
    if i <10
        nums = "000$i";
    elseif i<100
        nums = "00$i";
    else
        nums = "0$i";
    end
    files = "O016_$nums.h5";
    scatO[i,:] = h5read("$XSDIR/O016/$files","elastic");
    abspO[i,:] = h5read("$XSDIR/O016/$files","absorption");
    totalO[i,:] = h5read("$XSDIR/O016/$files","total");
end

println("<------Nuclear Data read------>")

totalMaxFe = zeros(length(energy));
totalMinFe = zeros(length(energy));
totalMaxO = zeros(length(energy));
totalMinO = zeros(length(energy));

for i = 1:length(energy)
    totalMaxFe[i]=findmax(totalFe[:,i])[1];
    totalMinFe[i]=findmin(totalFe[:,i])[1];
    totalMaxO[i]=findmax(totalO[:,i])[1];
    totalMinO[i]=findmin(totalO[:,i])[1];
end

total1 = Cross_section_Tendl(mt=1,reaction="total",energy_grid=energy, xs=totalFe);
totalMaxFe_xs = Cross_section(mt=1, reaction = "totalMax",energy_grid=energy, xs=totalMaxFe);
totalMinFe_xs = Cross_section(mt=1, reaction = "totalMin",energy_grid=energy, xs=totalMinFe);
scat1 = Cross_section_Tendl(mt=2,reaction="elastic_scatter",energy_grid=energy, xs=scatFe);
absp1 = Cross_section_Tendl(mt=27,reaction="absorption",energy_grid=energy, xs= abspFe);

total2 = Cross_section_Tendl(mt=1,reaction="total",energy_grid=energy, xs=totalO);
totalMaxO_xs = Cross_section(mt=1, reaction = "totalMax",energy_grid=energy, xs=totalMaxO);
totalMinO_xs = Cross_section(mt=1, reaction = "totalMin",energy_grid=energy, xs=totalMinO);
scat2 = Cross_section_Tendl(mt=2,reaction="elastic_scatter",energy_grid=energy, xs=scatO);
absp2 = Cross_section_Tendl(mt=27,reaction="absorption",energy_grid=energy, xs= abspO);

totalFe = 0;
totalMaxFe=0;
totalMinFe=0;
scatFe = 0;
abspFe = 0;

totalO = 0;
totalMaxO = 0;
totalMinO = 0;
scatO = 0;
abspO = 0;
GC.gc();

nuclide1 = Nuclide_Tendl(Name="Fe56",XS=[scat1,absp1],total_micro=total1,total_bounds=[totalMaxFe_xs,totalMinFe_xs], atomicWeight = 55.454479886265396);
nuclide2 = Nuclide_Tendl(Name="O16",XS=[scat2,absp2],total_micro=total2,total_bounds=[totalMaxO_xs,totalMinO_xs], atomicWeight = 15.857525022762783);

material1 = Material_Tendl(name="IronMixture", nucs=[nuclide1,nuclide2], atomic_density = [0.4,0.6], density = 0.1,id=1);

## Energy grid for tally
en_grid = [i for i=1e4:5e4:2e7];
tally1 = Flux_tally_pbox(energy_bins = en_grid)

## material and tally are properties of the juliaMC class

simulation1 = juliaMC(material=material1,n=10000, Tally_batch=tally1,n_batch=10 ,grid = energy)
#@time runMovie(simulation1)    ##option to make gif of sim, warning: takes for ever
println("Vanila MC")
a = @time runPar(simulation1,[1,1]);

simulation1 = juliaMC(material=material1,n=10000, Tally_batch=tally1,n_batch=10,grid = energy)
println("TMC")
b = @time runTotalMonteCarlo(simulation1,1000);

simulation1 = juliaMC(material=material1,n=100000, Tally_batch=tally1,n_batch=100,grid = energy)
println("FlySampling")
c = @time runFlySampling(simulation1);

#For creating plot
plotTally_pbox(b,c,a)

#For saving in csv format. HDF5 preferible but not currently compatible with risk cluster
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
