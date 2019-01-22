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
@everywhere using HDF5
@everywhere using LinearAlgebra             #@everywhere for distributed computing
include("juliaMC.jl");



# Function for Crossections, it is a linearly decreasing function with a gausian + cosine mixture
# A defines the mean of the gausian, shifts resonance peaks



## Two Crossections and 2 Nuclides. Crossections are properties of nuclides,
## and nuclides properties of the material class

numFe = 613;
numO = 641;

energy = h5read("../data_formated/Fe056/Fe056_0000.h5","energy");

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
    scatFe[i,:] = h5read("../data_formated/Fe056/$files","elastic");
    abspFe[i,:] = h5read("../data_formated/Fe056/$files","absorption");
    totalFe[i,:] = h5read("../data_formated/Fe056/$files","total");
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
    scatO[i,:] = h5read("../data_formated/O016/$files","elastic");
    abspO[i,:] = h5read("../data_formated/O016/$files","absorption");
    totalO[i,:] = h5read("../data_formated/O016/$files","total");
end

println("<------Nuclear Data read------>")

total1 = Cross_section_Tendl(mt=1,reaction="total",energy_grid=energy, xs=totalFe);
scat1 = Cross_section_Tendl(mt=2,reaction="elastic_scatter",energy_grid=energy, xs=scatFe);
absp1 = Cross_section_Tendl(mt=27,reaction="absorption",energy_grid=energy, xs= abspFe);

total2 = Cross_section_Tendl(mt=1,reaction="total",energy_grid=energy, xs=totalO);
scat2 = Cross_section_Tendl(mt=2,reaction="elastic_scatter",energy_grid=energy, xs=scatO);
absp2 = Cross_section_Tendl(mt=27,reaction="absorption",energy_grid=energy, xs= abspO);

totalFe = 0;
scatFe = 0;
abspFe = 0;

totalO = 0;
scatO = 0;
abspO = 0;
GC.gc();

nuclide1 = Nuclide(Name="Fe56",XS=[total1,scat1,absp1]);
nuclide2 = Nuclide(Name="O16",XS=[total2,scat2,absp2]);

material1 = Material(name="IronMixture", nucs=[nuclide1,nuclide2], atomic_density = [0.6,0.6], density = 4.0,id=1);

## Energy grid for tally
tally1 = Flux_tally(energy_bins = energy)

## material and tally are properties of the juliaMC class

simulation1 = juliaMC(material=material1,n=1000, Tally_batch=tally1,n_batch=10)
#@time runMovie(simulation1)    ##option to make gif of sim, warning: takes for ever
println("Vanila MC")
a = @time runPar(simulation1);

simulation1 = juliaMC(material=material1,n=100000, Tally_batch=tally1,n_batch=10)
println("TMC")
b = @time runTotalMonteCarlo(simulation1,100);

simulation1 = juliaMC(material=material1,n=100000, Tally_batch=tally1,n_batch=100)
println("FlySampling")
c = @time runFlySampling(simulation1);

#For creating plot
plotTally(b,c,a)

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
