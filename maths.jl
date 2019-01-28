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
    println(obj.weights)
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

    if (r>l)
        #=
        if (A[end-1] <= val && A[end]>val)
            return length(A)-1
        end
        if (A[1] <= val && A[2]>val)
            return 1
        end
        =#
        mid = floor(Int,l+(r-l)/2)
        #println(mid)
        if (A[mid] <= val && A[mid+1]>val)
            return mid
        end
        if (A[mid] > val)
            return binarySearch(val, A , l, mid)
        end
        return binarySearch(val, A, mid+1,r)
    end
    return -1

end

function binarySearch(val :: Float64, A :: Array{Float64,1})

    return binarySearch(val,A,1,length(A))

end

function rotate_angle(uvw :: Array{Float64,1}, mu :: Float64, phi :: Float64)
    u0 = uvw[1];    # original cosine in x direction
    v0 = uvw[1];    # original cosine in y direction
    w0 = uvw[2];    # original cosine in z direction

    sinphi = sin(phi);
    cosphi = cos(phi);

    a = sqrt(max(0., 1. - mu * mu));
    b = sqrt(max(0., 1. - w0 * w0));

    uvw1 = zeros(3);
    if b > 1.0e-10
        uvw[1] = mu * u0 + a * (u0 * w0 * cosphi - v0 * sinphi) / b;
        uvw[2] = mu * v0 + a * (v0 * w0 * cosphi + u0 * sinphi) / b;
        uvw[3] = mu * w0 - a * b * cosphi;
    else
        b = sqrt(1. - v0 * v0);
        uvw[1] = mu * u0 + a * (u0 * v0 * cosphi + w0 * sinphi) / b;
        uvw[2] = mu * v0 - a * b * cosphi;
        uvw[3] = mu * w0 + a * (v0 * w0 * cosphi - u0 * sinphi) / b;
    end

    return uvw
end
