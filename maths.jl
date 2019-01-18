####
#
#   This script conatains an empirical cdf class and sampler + very simple search algorithm
#
#   Julia Version: V1.0
####


# Descrete CDF class for constructing and empirical CDF of reactions
# members are cdf elements and weights are probability measures
@with_kw mutable struct Descrete_CDF
    members :: Array{Any,1}
    weights :: Array{Float64,1}
    #Constructor
    #=
    function Descrete_CDF(mems, probs :: Array{Float64,1})

        new(mems, pushfirst!(probs,0))
    end
=#
end
# A simple sampler for the cdf class. Returns one of the elements with probability of the wights
# This is a function that is unique to the Descrete_CDF class, in julia this is the only way to tie a
# function to an object.
function (obj :: Descrete_CDF)()

    sample = rand()

    i=binarySearch(sample,obj.weights)

    return obj.members[i], i

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

function binarySearch(val :: Float64, A :: Array{Float64,1}, l :: Int, r :: Int)

    if (r>=l)
        mid = floor(Int,l+(r-l)/2)
        if (A[end-1] <= val && A[end]>val)
            return length(A)-1
        end
        mid = floor(Int,l+(r-l)/2)
        if (A[mid] <= val && A[mid+1]>val)
            return mid
        end
        if (A[mid] > val)
            return binarySearch(val, A , l, mid-1)
        end
        return binarySearch(val, A, mid+1,r)
    end
    return throw(UndefVarError(:binarySearchFailed))

end

function binarySearch(val :: Float64, A :: Array{Float64,1})

    return binarySearch(val,A,1,length(A))

end
