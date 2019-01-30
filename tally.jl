####
#
#   Tally class and plotters
#
#   Julia Version: V1.0
####
@with_kw mutable struct Flux_tally

    n :: Int = 10                                           # This must be the number of batches in the simulation
    radius :: Float64 = 100
    centre :: Array{Float64,1} = [0,0,0]
    energy_bins :: Array{Float64,1} = [0.0,]                # The width of the bins to be tallied for flux
    Tally = zeros(length(energy_bins)-1)
    volume :: Float64 = 4/3*pi*radius^3
    std :: Array{Float64,1} = zeros(length(energy_bins)-1)

end

@with_kw mutable struct Flux_tally_pbox

    n :: Int = 10                                           # This must be the number of batches in the simulation
    radius :: Float64 = 100
    centre :: Array{Float64,1} = [0,0,0]
    energy_bins :: Array{Float64,1} = [0.0,]                # The width of the bins to be tallied for flux
    Tally = zeros(length(energy_bins)-1)
    Tally_bounds = zeros(2,length(energy_bins)-1)
    std :: Array{Float64,1} = zeros(length(energy_bins)-1)
    std_bounds = zeros(2,length(energy_bins)-1)
    volume :: Float64 = 4/3*pi*radius^3
end




# Function for easy plotting of tally, here we will find the centre of thebins to plot against

function plotTally(Tal :: Flux_tally)


    en = (Tal.energy_bins[2:end]+Tal.energy_bins[1:end-1])/2

    if sum(Tal.std) == 0.0
        plot!(en, Tal.Tally)
    else
        plot!(en,Tal.Tally, ribbon=(Tal.std,Tal.std), title="Flux", fillalpha = 0.5, linewidth = 2)
    end

end



function plotTally(Tal1 :: Flux_tally, Tal2 :: Flux_tally, Tal3 :: Flux_tally)


    # Tal1 == TMC
    # Tal2 == FLY
    # Tal3 == Vinila
    en = (Tal1.energy_bins[2:end]+Tal1.energy_bins[1:end-1])/2;


    for i = 1:length(Tal1.Tally)

        if Tal1.Tally[i] == 0.0
            Tal1.Tally[i] = eps();
        end
        if Tal2.Tally[i] == 0.0
            Tal2.Tally[i] = eps();
        end
        if Tal3.Tally[i] == 0.0
            Tal3.Tally[i] = eps();
        end

    end

    a=Tal2.Tally./Tal1.Tally;
    b=Tal2.Tally./Tal3.Tally;

    CoV1 = Tal1.std./Tal1.Tally
    CoV2 = Tal2.std./Tal2.Tally
    CoV3 = Tal3.std./Tal3.Tally

    for i =1:length(en)

        if isnan(a[i])
            a[i]=0;
            CoV1[i]=0;
        end

        if isnan(b[i])
            b[i]=0;
        end

        if isnan(CoV2[i])
            CoV2[i]=0;
        end
        if isnan(CoV3[i])
            CoV3[i]=0;
        end

    end

    stds1 = sqrt.(CoV1.^2 .+ CoV2.^2);
    stds2 = sqrt.(CoV3.^2 .+ CoV2.^2);
    #stds1 = 3.*(Tal1.std + Tal2.std)

    #stds = hcat(zeros(length(en)),zeros(length(en)),stds1,stds1)

    tals = hcat(Tal1.Tally,Tal2.Tally)
    vars = hcat(Tal1.std,Tal2.std)

    tals2 = hcat(Tal3.Tally,Tal2.Tally)
    vars2 = hcat(Tal3.std,Tal2.std)


    a1 = hcat(a,ones(length(en)))

    y = ones(length(en))

    plt1 = plot(en, y, ribbon=(stds1,stds1),fillalpha = 0.5, label="CoV + 1");
    plt1 = plot!(en,a, title="means FLY/TMC", linewidth = 2, label="means");
    #plt2 = plot(en,b, title="stds FLY/TMC");

    plt2 = plot(en,Tal1.Tally, ribbon=(Tal1.std,Tal1.std), title="Flux plots", fillalpha = 0.5, label = "TMC", linewidth = 2, xaxis = :log, yaxis = :log);
    plt2 = plot!(en,Tal2.Tally, linewidth = 2, ribbon=(Tal2.std,Tal2.std), fillalpha = 0.5, label="Fly", xaxis = :log, yaxis = :log);

    plt3 = plot(en, y, ribbon=(stds2,stds2),fillalpha = 0.5,label="CoV + 1");
    plt3 = plot!(en, b, title="means FLY/Vanila", linewidth = 2,label="means");


    plt4 = plot(en,Tal3.Tally, ribbon=(Tal3.std,Tal3.std), title="Flux plots", fillalpha = 0.5, label = "Vanila", linewidth = 2, xaxis = :log , yaxis = :log);
    plt4 = plot!(en,Tal2.Tally, linewidth = 2, ribbon=(Tal2.std,Tal2.std), fillalpha = 0.5, label="Fly", xaxis = :log , yaxis = :log);

    fig=plot(plt1,plt2,plt3,plt4, dpi=300, size=(1000,1000))

    display(fig)
    savefig(fig,"Result100_MIl_batch.png")

    #plot(plt1,plt2,plt3,plt4)
    #=
    savefig(plt1, "means.png")
    savefig(plt2, "stds.png")
    savefig(plt3, "TMC&FLY.png")
    savefig(plt4, "Vanila&FLY.png")
    =#
end

function plotTally_pbox(Tal1 :: Flux_tally_pbox, Tal2 :: Flux_tally_pbox, Tal3 :: Flux_tally)


    # Tal1 == TMC
    # Tal2 == FLY
    # Tal3 == Vinila
    en = (Tal1.energy_bins[2:end]+Tal1.energy_bins[1:end-1])/2;


    for i = 1:length(Tal1.Tally)

        if Tal1.Tally[i] == 0.0
            Tal1.Tally[i] = eps();
        end
        if Tal2.Tally[i] == 0.0
            Tal2.Tally[i] = eps();
        end
        if Tal3.Tally[i] == 0.0
            Tal3.Tally[i] = eps();
        end
        if Tal1.Tally_bounds[i] == 0.0
            Tal1.Tally_bounds[i] = eps();
        end
        if Tal2.Tally_bounds[i] == 0.0
            Tal2.Tally_bounds[i] = eps();
        end


    end

    a=Tal2.Tally./Tal1.Tally;
    b=Tal2.Tally./Tal3.Tally;

    CoV1 = Tal1.std./Tal1.Tally
    CoV2 = Tal2.std./Tal2.Tally
    CoV3 = Tal3.std./Tal3.Tally

    for i =1:length(en)

        if isnan(a[i])
            a[i]=0;
            CoV1[i]=0;
        end

        if isnan(b[i])
            b[i]=0;
        end

        if isnan(CoV2[i])
            CoV2[i]=0;
        end
        if isnan(CoV3[i])
            CoV3[i]=0;
        end

    end

    stds1 = sqrt.(CoV1.^2 .+ CoV2.^2);
    stds2 = sqrt.(CoV3.^2 .+ CoV2.^2);
    #stds1 = 3.*(Tal1.std + Tal2.std)

    #stds = hcat(zeros(length(en)),zeros(length(en)),stds1,stds1)

    tals = hcat(Tal1.Tally,Tal2.Tally)
    vars = hcat(Tal1.std,Tal2.std)

    tals2 = hcat(Tal3.Tally,Tal2.Tally)
    vars2 = hcat(Tal3.std,Tal2.std)


    a1 = hcat(a,ones(length(en)))

    y = ones(length(en))

    plt1 = plot(en, y, ribbon=(stds1,stds1),fillalpha = 0.5, label="CoV + 1");
    plt1 = plot!(en,a, title="means FLY/TMC", linewidth = 2, label="means");
    #plt2 = plot(en,b, title="stds FLY/TMC");

    plt2 = plot(en,Tal1.Tally, ribbon=(Tal1.Tally_bounds[1,:]+2*Tal1.std_bounds[1,:],Tal1.Tally_bounds[2,:]+2*Tal1.std_bounds[2,:]), title="Flux plots", fillalpha = 0.5, label = "TMC", linewidth = 1.2, xaxis = :log);
    plt2 = plot!(en,Tal2.Tally, linewidth = 2, ribbon=(Tal2.Tally_bounds[1,:]+2*Tal2.std_bounds[1,:],Tal2.Tally_bounds[2,:]+2*Tal2.std_bounds[2,:]), fillalpha = 0.5, label="Fly", xaxis = :log);

    plt3 = plot(en, y, ribbon=(stds2,stds2),fillalpha = 0.5,label="CoV + 1");
    plt3 = plot!(en, b, title="means FLY/Vanila", linewidth = 2,label="means");


    plt4 = plot(en,Tal3.Tally, ribbon=(Tal3.std,Tal3.std), title="Flux plots", fillalpha = 0.5, label = "Vanila", linewidth = 1.2, xaxis = :log);
    plt4 = plot!(en,Tal2.Tally, linewidth = 2, ribbon=(Tal2.Tally_bounds[1,:]+2*Tal2.std_bounds[1,:],Tal2.Tally_bounds[2,:]+2*Tal2.std_bounds[2,:]), fillalpha = 0.5, label="Fly", xaxis = :log );

    fig=plot(plt1,plt2,plt3,plt4, dpi=300, size=(1000,1000))

    display(fig)
    savefig(fig,"Result100_MIl_batch.png")

end



function plotTally_pbox(Tal1 :: Flux_tally_pbox, Tal2 :: Flux_tally_pbox)


    # Tal1 == TMC
    # Tal2 == FLY
    # Tal3 == Vinila
    en = (Tal1.energy_bins[2:end]+Tal1.energy_bins[1:end-1])/2;


    for i = 1:length(Tal1.Tally)

        if Tal1.Tally[i] == 0.0
            Tal1.Tally[i] = eps();
        end
        if Tal2.Tally[i] == 0.0
            Tal2.Tally[i] = eps();
        end
        if Tal1.Tally_bounds[i] == 0.0
            Tal1.Tally_bounds[i] = eps();
        end
        if Tal2.Tally_bounds[i] == 0.0
            Tal2.Tally_bounds[i] = eps();
        end


    end

    a=Tal2.Tally./Tal1.Tally;

    CoV1 = Tal1.std./Tal1.Tally
    CoV2 = Tal2.std./Tal2.Tally

    for i =1:length(en)

        if isnan(a[i])
            a[i]=0;
            CoV1[i]=0;
        end

        if isnan(b[i])
            b[i]=0;
        end

        if isnan(CoV2[i])
            CoV2[i]=0;
        end
    end

    stds1 = CoV1;

    #sqrt.(CoV1.^2 .+ CoV2.^2);
    #stds1 = 3.*(Tal1.std + Tal2.std)

    #stds = hcat(zeros(length(en)),zeros(length(en)),stds1,stds1)

    tals = hcat(Tal1.Tally,Tal2.Tally)
    vars = hcat(Tal1.std,Tal2.std)

    tals2 = hcat(Tal3.Tally,Tal2.Tally)
    vars2 = hcat(Tal3.std,Tal2.std)


    a1 = hcat(a,ones(length(en)))

    y = ones(length(en))

    plt1 = plot(en, y, ribbon=(stds1,stds1),fillalpha = 0.5, label="CoV + 1");
    plt1 = plot!(en,a, title="means FLY/TMC", linewidth = 2, label="means");
    #plt2 = plot(en,b, title="stds FLY/TMC");

    plt2 = plot(en,Tal1.Tally, ribbon=(Tal1.Tally_bounds[1,:]+2*Tal1.std_bounds[1,:],Tal1.Tally_bounds[2,:]+2*Tal1.std_bounds[2,:]), title="Flux plots", fillalpha = 0.5, label = "TMC", linewidth = 1.2, xaxis = :log);
    plt2 = plot!(en,Tal2.Tally, linewidth = 2, ribbon=(Tal2.Tally_bounds[1,:]+2*Tal2.std_bounds[1,:],Tal2.Tally_bounds[2,:]+2*Tal2.std_bounds[2,:]), fillalpha = 0.5, label="Fly", xaxis = :log);

    LowerFly = Normal(Tally2.Tally_bounds[1,end/2],Tally2.std_bounds[1,end/2]);

    plt3 = plot(en, y, ribbon=(stds2,stds2),fillalpha = 0.5,label="CoV + 1");
    plt3 = plot!(en, b, title="means FLY/Vanila", linewidth = 2,label="means");


    plt4 = plot(en,Tal3.Tally, ribbon=(Tal3.std,Tal3.std), title="Flux plots", fillalpha = 0.5, label = "Vanila", linewidth = 1.2, xaxis = :log);
    plt4 = plot!(en,Tal2.Tally, linewidth = 2, ribbon=(Tal2.Tally_bounds[1,:]+2*Tal2.std_bounds[1,:],Tal2.Tally_bounds[2,:]+2*Tal2.std_bounds[2,:]), fillalpha = 0.5, label="Fly", xaxis = :log );

    fig=plot(plt1,plt2,plt3,plt4, dpi=300, size=(1000,1000))

    display(fig)
    savefig(fig,"Result100_MIl_batch.png")

end
