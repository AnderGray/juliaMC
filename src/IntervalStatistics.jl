##
#   Written By Ander Gray. Method fort computing interval mean and variance from and interval data set.
#   As described in S. Ferson (2007), "Experimental Uncertainty Estimation and Statistics for Data Having Interval Uncertainty"
#
##

using Statistics;

function IntervalStatistics(samples)
    
    if isequal(samples[:,1],samples[:,2])
        means = mean(samples[:,1]);
        vars = var(samples[:,1]);
        intervalMean = [means, means];
        intervalVariance = [vars,vars];
        return [intervalMean, intervalVariance];
     end
 
     intervalMean = [mean(samples[:,1]),mean(samples[:,2])];
     intervalVariance = [0.0,0.0];

     samplesSorted = zeros(size(samples));
     
     samplesSorted[:,1] = sort(samples[:,1]);
     samplesSorted[:,2] = sort(samples[:,2]);

     Ys = reshape(samplesSorted',length(samples),1);

     Ys = sort(Ys[:]);
 
     overlap = zeros(length(Ys)-1,1);
     Ns      = zeros(length(overlap),1);
     Sks     = zeros(length(overlap),1);
     Ms      = zeros(length(overlap),1);
 
     N = length(samples)/2;
 
     VarNull = 0;
     for i=1:length(overlap)
         sorted = sort([Ys[i],Ys[i+1]]);
         if sorted[1]>intervalMean[2]
         elseif sorted[2]<intervalMean[1]
         else
             overlap[i]=1;
             leftBounds = samplesSorted[:,2] .<= Ys[i];
             rightBounds = samplesSorted[:,1] .>= Ys[i+1];
             Ns[i] = sum(leftBounds) + sum(rightBounds);
             if Ns[i] == 0
                 VarNull = 1;
                 break;
             end
             leftOk   = samplesSorted[:,2] .* leftBounds;
             rightOk  = samplesSorted[:,1] .* rightBounds;
 
             #left = leftOk(leftOk ~=0);
             #right = rightOk(rightOk ~=0);
             left = filter!(x->x≠0,leftOk)
             right = filter!(x->x≠0,rightOk)
             Sks[i] = sum(left) + sum(right);
             test = Sks[i]/Ns[i];
             if sorted[1]<= test && sorted[2] >= test
                 if intervalMean[1]<= test && intervalMean[2] >= test
                     Ms[i] = (sum(left.^2) + sum(right.^2))/(length(samplesSorted)/2);
                 end
             end 
         end
     end
     
     # Computing Lower Bound of Variance
     if VarNull == 0
         MsOks      = ones(length(Ms),1) .* Ms;
         SksFinal0  = Sks .* MsOks;
         NsFinal0   = Ns .* MsOks;
 
         #MsFinal1   = Ms(Ms ~= 0);
         #SksFinal1  = SksFinal0(SksFinal0 ~=0)./MsFinal1;
         #NsFinal1   = NsFinal0(NsFinal0 ~=0)./MsFinal1;
         
         MsFinal1 = filter!(x->x≠0,vec(Ms));
         SksFinal1  = filter!(x->x≠0,vec(SksFinal0))./MsFinal1;
         NsFinal1   = filter!(x->x≠0.0,vec(NsFinal0))./MsFinal1;

         vars = MsFinal1 - SksFinal1.^2 ./((N .* NsFinal1));
         intervalVariance[1] = minimum(vars) * (N/(N-1));
     end

     # Computing Upper Bound of Variance
     # Better way to cart prod julia?
     # p = num2cell(samples,2);
     # cart = cartprod(p{:});
     # a = collect(Iterators.product(skinny[1,:],skinny[2,:],skinny[3,:]))
     
     cart = cartesianProd(samples);

     vars = var.(cart');
     intervalVariance[2] = maximum(vars);

     return intervalMean , intervalVariance
 

end


# Cartesian Product algorithm
function cartesianProd(Sets1)

    results = [[]];

    SetNum = length(Sets1[:,1]);

    for i=1:SetNum 
        currentSubArray  = Sets1[i,:];
        temp = [];
        for j=1:length(results) 
            for k =1:length(currentSubArray)
                push!(temp, vcat(results[j],currentSubArray[k]));
            end 
        end 
        results = temp;
    end
    return results;
end
