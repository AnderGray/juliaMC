####
#
#   This script conatains an empirical cdf class and sampler + very simple search algorithm
#
#   Julia Version: V1.0
####


# Descrete CDF class for constructing and empirical CDF of reactions
# members are cdf elements and weights are probability measures
@with_kw mutable struct Descrete_CDF
    members :: Array{Any, 1}
    weights :: Array{Float64,1}
end

# A simple sampler for the cdf class. Returns one of the elements with probability of the wights
# This is a function that is unique to the Descrete_CDF class, in julia this is the only way to tie a
# function to an object.
function (obj :: Descrete_CDF)()

    sample = rand()

    i=findInter(sample,obj.weights)

    return obj.members[i+1], i+1

end

# A simple search algorithm for finding where val is in array A.
function findInter(val :: Float64, A :: Array{Float64,1})
    iterand=-1
    for i=1:length(A)
        if val<A[i]
            iterand=i-1;
            break
        end
    end

    if iterand==-1
        # launch error here
    end
    return iterand
end
