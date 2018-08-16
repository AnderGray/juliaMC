include("juliaMC.jl");

en_grid = [i for i=1e6:2e5:3e8];
function xs(x,A) 
    if x<14.e7
        y=(2e8*cos.(3.0*x).*pdf.(Normal(A,1e7),x)-5e-8x +3)+20;
    else
        y= 1000
    end
end

scat1 = Cross_section(mt=1,reaction="elastic_scatter",energy_grid=en_grid, xs= xs.(en_grid,1e8));
absp1 = Cross_section(mt=2,reaction="absorption",energy_grid=en_grid, xs= xs.(en_grid,5e7));

scat2 = Cross_section(mt=1,reaction="elastic_scatter",energy_grid=en_grid, xs= xs.(en_grid,1.5e8));
absp2 = Cross_section(mt=2,reaction="absorption",energy_grid=en_grid, xs= xs.(en_grid,8e7));
#=
scat1 = Cross_section(mt=1,reaction="elastic_scatter",energy_grid=en_grid, xs= xs1.(en_grid));
absp1 = Cross_section(mt=2,reaction="absorption",energy_grid=en_grid, xs= xs1.(en_grid)*0.01);

scat2 = Cross_section(mt=1,reaction="elastic_scatter",energy_grid=en_grid, xs= xs1.(en_grid));
absp2 = Cross_section(mt=2,reaction="absorption",energy_grid=en_grid, xs= xs1.(en_grid)*0.01);
=#
nuclide1 = Nuclide(Name="Fe56",XS=[scat1,absp1]);
nuclide2 = Nuclide(Name="Fe54",XS=[scat2,absp2]);

material1 = Material(name="IronMixture", nucs=[nuclide1,nuclide2], atomic_density = [0.6,0.6], density = 4.0,id=1);

tally_grid = [i for i=1e6:2e6:3e8];

tally1 = Flux_tally(energy_bins = tally_grid, radius=100)

material1 = Material(name="IronMixture", nucs=[nuclide1,nuclide2], atomic_density = [0.6,0.6], density = 4.0,id=1);


simulation1 = juliaMC(material=material1,n=100000, Tally_batch=tally1,n_batch=10)

#@time runMovie(simulation1)


println("Vanila MC")
a = @time runPar(simulation1);

#plotTally(a)
tally1 = Flux_tally(energy_bins = tally_grid)
simulation1 = juliaMC(material=material1,n=100000, Tally_batch=tally1,n_batch=10)

println("TMC")
b = @time runTotalMonteCarlo(simulation1,10);
#b= @time runFlySampling(simulation1);

tally1 = Flux_tally(energy_bins = tally_grid)
simulation1 = juliaMC(material=material1,n=100000, Tally_batch=tally1,n_batch=10)


println("FlySampling")
c = @time runFlySampling(simulation1);

plotTally(b,c,a)

#plotTally(c)

#=
#a=[1,5,10,50,100,500,1000,5000]

a = [1000000,5000000,10000000,50000000,100000000,500000000]
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

answer = plot(a,tal, dpi=300, size=(1000,1000));

savefig(answer, "Answer3.png")
=#
