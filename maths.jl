@with_kw mutable struct Descrete_CDF
    
    members :: Array{Any, 1}
    weights :: Array{Float64,1}


end


function (obj :: Descrete_CDF)()
    
    sample = rand()
    
    i=findInter(sample,obj.weights)
    
    return obj.members[i+1], i+1
    
end



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