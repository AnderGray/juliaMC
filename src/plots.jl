
place = 1.1;



meanFly = Normal(c.Tally[floor(Int,end/place)],c.std[floor(Int,end/place)])
lowerFly = Normal(c.Tally_bounds[1,floor(Int,end/place)],c.std_bounds[1,floor(Int,end/place)])
upperFly = Normal(c.Tally_bounds[2,floor(Int,end/place)],c.std_bounds[2,floor(Int,end/place)])

meanFlyRange = collect(range(meanFly.μ - 3*meanFly.σ, length=500, stop=meanFly.μ + 3*meanFly.σ))
lowerFlyRange = collect(range(lowerFly.μ - 3*lowerFly.σ, length=500, stop=lowerFly.μ + 3*lowerFly.σ))
upperFlyRange = collect(range(upperFly.μ - 3*upperFly.σ, length=500, stop=upperFly.μ + 3*upperFly.σ))

meanFlyCdf = cdf.(meanFly,meanFlyRange)
lowerFlyCdf = cdf.(lowerFly,lowerFlyRange)
upperFlyCdf = cdf.(upperFly,upperFlyRange)

meanTMC = Normal(b.Tally[floor(Int,end/place)],b.std[floor(Int,end/place)])
lowerTMC = Normal(b.Tally_bounds[1,floor(Int,end/place)],b.std_bounds[1,floor(Int,end/place)])
upperTMC = Normal(b.Tally_bounds[2,floor(Int,end/place)],b.std_bounds[2,floor(Int,end/place)])

meanTMCRange = collect(range(meanTMC.μ - 3*meanTMC.σ, length=500, stop=meanTMC.μ + 3*meanTMC.σ))
lowerTMCRange = collect(range(lowerTMC.μ - 3*lowerTMC.σ, length=500, stop=lowerTMC.μ + 3*lowerTMC.σ))
upperTMCRange = collect(range(upperTMC.μ - 3*upperTMC.σ, length=500, stop=upperTMC.μ + 3*upperTMC.σ))

meanTMCCdf = cdf.(meanTMC,meanTMCRange)
lowerTMCCdf = cdf.(lowerTMC,lowerTMCRange)
upperTMCCdf = cdf.(upperTMC,upperTMCRange)

plt1=plot(meanFlyRange,meanFlyCdf,color="red")
plt1=plot!(lowerFlyRange,lowerFlyCdf,color="red")
plt1=plot!(upperFlyRange,upperFlyCdf,color="red")

plt1=plot!(meanTMCRange,meanTMCCdf,color="blue")
plt1=plot!(lowerTMCRange,lowerTMCCdf,color="blue")
plt1=plot!(upperTMCRange,upperTMCCdf,color="blue")

fig=plot(plt1, dpi=300, size=(1000,1000))
display(fig)
savefig(fig,"Pboxes.png")


#2 == 9.9e6
#4 == 5e6
#7 == 2.81e6
#120 == 110000.0
