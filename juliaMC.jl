

@everywhere using Parameters;
@everywhere using Distributions;
#using PyPlot;
@everywhere using Plots;
gr()

@everywhere include("nuclearData.jl");
@everywhere include("material.jl");
@everywhere include("particle.jl");
@everywhere include("source.jl");
@everywhere include("physics.jl");
@everywhere include("maths.jl");
@everywhere include("tally.jl");
    

#=
Object that holds the simulation parameters, the constructor ensures that the dimension of tally is the same as the number of batches in the simulation.
The parameters are as follos:

n: Number of particles per batch

n_batch: The number of batches

source: The source of particles, set to default see source.jl
material: The material to be used in the simulation. Only one can be used for now and this is added to the particles when generated

N_bank: An array of particles that is used for the simulation. This is generated using the generate() function, see source.jl However this will probably need to be changed due to batching. Either 2 dim array of columbs = n_batch or that the particles are generated at the begining of run().

Tall_batch: The Tally to be used in the simulation

=#

@everywhere @with_kw struct juliaMC
    
    n :: Int64 = 1000                                                    
    n_batch :: Int64 = 10                                                
    source :: Source = Source()                                          
    material :: Material                                                
    #N_bank :: Array{Particle, 2} = generate(source, material, n, n_batch)
    Tally_batch :: Flux_tally
    
    function juliaMC(n :: Int64, n_batch :: Int64, source :: Source, material :: Material, Tally_batch :: Flux_tally)
        
        Tally = Flux_tally(n=n_batch, energy_bins = Tally_batch.energy_bins,radius = Tally_batch.radius)
        
        new(n, n_batch, source, material, Tally)
    end
    
end




@everywhere function run(sim :: juliaMC)
    
    en = (sim.Tally_batch.energy_bins[2:end]+sim.Tally_batch.energy_bins[1:end-1])/2
    
    Tal = Array{Float64,2}(length(sim.Tally_batch.energy_bins)-1,sim.n_batch)

    
    
    for j=1:sim.n_batch
        N_bank = generate(sim.source,sim.material,sim.n)
        localTal = zeros(length(sim.Tally_batch.energy_bins)-1)
        for i=1:sim.n                                
            while N_bank[i].alive == true
                N_bank[i] = transport(N_bank[i]);
                if norm(N_bank[i].xyz)>sim.Tally_batch.radius
                    N_bank[i].alive=false;
                else
                    k = findInter(N_bank[i].last_E, sim.Tally_batch.energy_bins);
                    localTal[k] += N_bank[i].last_d*N_bank[i].wgt
                end
                end
            end
        Tal[:,j] = localTal;
        println("Finished Batch $j")
        #println(localTal)
    end
    
    Tal=Tal./sim.Tally_batch.volume;
    #plot(en, Tal)
    Global_tally = Flux_tally(energy_bins=sim.Tally_batch.energy_bins,n=1);
    
    for i = 1:length(sim.Tally_batch.Tally[:,1])
        
        Global_tally.Tally[i] = mean(Tal[i,:]);
        Global_tally.std[i] = std(Tal[i,:]);
        
    end
    

    
    plotTally(Global_tally)
    #return Global_tally
end



function runPar(sim :: juliaMC)
    
    #en = (sim.Tally_batch.energy_bins[2:end]+sim.Tally_batch.energy_bins[1:end-1])/2
    
    Tal = SharedArray{Float64,2}(length(sim.Tally_batch.energy_bins)-1,sim.n_batch)
    counter = SharedArray{Float64,1}(1)
    counter = 0;
    
    @sync @parallel for j=1:sim.n_batch
        #N_bank = generate(sim.source,sim.material,sim.n)
        localTal = zeros(length(sim.Tally_batch.energy_bins)-1)
        N_bank = generate(sim.source,sim.material,10000)
        o = 1;
        for i=1:sim.n   
            if i % 10000 == 0
                N_bank = generate(sim.source,sim.material,10000)
                o = 1
            end
            while N_bank[o].alive == true
                N_bank[o] = transport(N_bank[o]);
                if norm(N_bank[o].xyz)>sim.Tally_batch.radius
                    N_bank[o].alive=false;
                else
                    k = findInter(N_bank[o].last_E, sim.Tally_batch.energy_bins);
                    localTal[k] += N_bank[o].last_d*N_bank[o].wgt
                end
                end
            o+=1;
            end
        counter +=1
        Tal[:,j] = localTal;
        println("Finished Batch $counter")
        #println(localTal)
    end
    
    Tal=Tal./(sim.Tally_batch.volume*sim.n);
    #plot(en, Tal)
    Global_tally = Flux_tally(energy_bins=sim.Tally_batch.energy_bins,n=1);
    
    for i = 1:length(sim.Tally_batch.Tally[:,1])
        
        Global_tally.Tally[i] = mean(Tal[i,:]);
        Global_tally.std[i] = std(Tal[i,:]);
        
    end
    
    #plotTally(Global_tally)
    return Global_tally
end




function runTotalMonteCarlo(sim :: juliaMC , n :: Int64)
    
    m = length(sim.Tally_batch.energy_bins)-1
    Global = Flux_tally(energy_bins=sim.Tally_batch.energy_bins,n=1)
    
    means = zeros(m,n)
    stds = zeros(m,n)
    
    for i =1:n
        simulation = deepcopy(sim)
        for j = 1:sim.material.n_nuclides
            for k =1:length(sim.material.nucs[j].XS)
                simulation.material.nucs[j].XS[k].xs=sim.material.nucs[k].XS[k].xs*rand(Truncated(Normal(1,0.5),0.3,100))
            end
        end
        Tally = runPar(simulation)
        
        println("<------Finished Iteration $i------>")
        
        means[:,i] = Tally.Tally
        stds[:,i] = Tally.std
    end
    for l=1:m
        Global.Tally[l,:] = mean(means[l,:])
        Global.std[l,:] = sqrt(mean(stds[l,:].^2) + std(means[l,:]).^2)
    end
    
    Tally=0;
    means=0;
    stds=0;
    simulation=0;
    sim=0;
    gc()
    
    return Global
    
end


function runFlySampling(sim :: juliaMC)
    en = (sim.Tally_batch.energy_bins[2:end]+sim.Tally_batch.energy_bins[1:end-1])/2
    
    Tal = SharedArray{Float64,2}(length(sim.Tally_batch.energy_bins)-1,sim.n_batch)
    
    
    @sync @parallel for j=1:sim.n_batch
        #simulation = deepcopy(sim)
        #N_bank = generate(sim.source,sim.material,sim.n)
        localTal = zeros(length(sim.Tally_batch.energy_bins)-1)
        
        N_bank = generate(sim.source,sim.material,10000)
        o = 1;
        for i=1:sim.n   
            if i % 10000 == 0
                N_bank = generate(sim.source,sim.material,10000)
                o = 1
            end
            perturb = rand(Truncated(Normal(1,2),0.3,100))
            while N_bank[o].alive == true
                N_bank[o] = transportUQ(N_bank[o],perturb);
                if norm(N_bank[o].xyz)>sim.Tally_batch.radius
                    N_bank[o].alive=false;
                else
                    k = findInter(N_bank[o].last_E, sim.Tally_batch.energy_bins);
                    localTal[k] += N_bank[o].last_d*N_bank[o].wgt
                end
                end
            o+=1
            end
        Tal[:,j] = localTal;
        println("Finished Batch $j")
        #println(localTal)
    end
    
    Tal=Tal./(sim.Tally_batch.volume*sim.n);
    Global_tally = Flux_tally(energy_bins=sim.Tally_batch.energy_bins,n=1)
    
    for i = 1:length(sim.Tally_batch.Tally[:,1])
        
        Global_tally.Tally[i] = mean(Tal[i,:])
        Global_tally.std[i] = std(Tal[i,:])
        
    end
    #plot!(en, Tal)
    #println(Tal)
    
    return Global_tally
end





function runMovie(sim :: juliaMC)
    
    en = (sim.Tally_batch.energy_bins[2:end]+sim.Tally_batch.energy_bins[1:end-1])/2
    
    Tally = zeros(length(en), sim.n_batch)
    o = 1;
    tic()
    time = 0;
    w=0
    FrameNum=0;
    plt2=plot()
    for j=1:sim.n_batch
        N_bank = generate(sim.source,sim.material,10000)
        for i=1:sim.n
            if i % 10000 == 0
                N_bank = generate(sim.source,sim.material,10000)
                o = 1
            end
            perturb = rand(Truncated(Normal(1,0.5),0.3,100))
            while N_bank[o].alive == true
                N_bank[o] = transportUQ(N_bank[o],perturb);
                if norm(N_bank[o].xyz)>sim.Tally_batch.radius
                    N_bank[o].alive=false;
                else
                    k = findInter(N_bank[o].last_E, sim.Tally_batch.energy_bins);
                    Tally[k,j] += N_bank[o].last_d*N_bank[o].wgt
                    t=toq();
                    time +=t;
                    if w%100000 == 0
                        if FrameNum<10
                            s ="0000$FrameNum"
                        elseif FrameNum<100
                            s="000$FrameNum"
                        elseif FrameNum<1000
                            s="00$FrameNum"
                        elseif FrameNum<10000
                            s="0$FrameNum"
                        else
                            s="$FrameNum"
                        end
                        #println(j)
                        if j<2
                            plt=plot(en, Tally[:,j], title= "Fly sampling. Batch = $j | T = $time s",linewidth = 2);
                        else
                            b = Tally[:,1]
                            for h = 2:j
                                #println(h)
                                b = hcat(b,Tally[:,h])
                            end
                            if j<3
                                plt=plot(en,b,title= "Fly sampling. Batch = $j | T = $time s",linewidth = 2)
                            else
                                c=zeros(length(en))
                                standardDiv=zeros(length(en))
                                for z=1:length(en)
                                    c[z] = mean(b[z,1:j-1])
                                    standardDiv[z] = std(b[z,1:j-1])
                                end
                                plot(en,c,ribbon=(standardDiv,standardDiv),linewidth = 2, fillalpha=0.5)
                                plt=plot!(en,Tally[:,j],linewidth = 2, title= "Fly sampling. Batch = $j | T = $time s");
                            end
                        end
                        #ylims!((0,5))
                        #yaxis!("arb FLux",:log10)
                        savefig(plt, "gifs/Frame""$s"".png");
                        FrameNum+=1
                    end
                    w+=1
                    tic()
                end
            end
            o+=1
        end
        println("Finished Batch $j")
        b = Tally[:,1]
        if j>1
            for h = j:-1:2
                b = hcat(b,Tally[:,h])
            end
        end
        plt2=plot(en,b)
        #sim.Tally_batch.Tally[:,j]=sim.Tally_batch.Tally[:,j]/sim.Tally_batch.volume;
    end
    #ylims!((0,5))
    savefig(plt2, "gifs/FrameFinal"".png");
    
    Global_tally = Flux_tally(energy_bins=sim.Tally_batch.energy_bins,n=1)
    
    for i = 1:length(Tally[:,1])
        
        Global_tally.Tally[i] = mean(Tally[i,:])
        Global_tally.std[i] = std(Tally[i,:])
        
    end
    
    println(Global_tally.Tally)
    
    #plotTally(Global_tally)
    return Global_tally
end


