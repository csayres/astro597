using Plots
using DelimitedFiles
using Statistics
using FFTW
using Optim


#################################################################################
####################################################################################
#TRANSIT Model stuff

struct Circle
    r::Float64
    x::Float64
    y::Float64
end

# https://stackoverflow.com/questions/4247889/area-of-intersection-between-two-circles
# this ones more stable
function circOverlap(A, B)
    r = A.r
    R = B.r
    # d = hypot(B.x - A.x, B.y - A.y)
    d = sqrt((B.x-A.x)^2 + (B.y - A.y)^2) # this is faster?

    if d >= r + R
        return 0
    end
    # check next case in which one circle is completely inside
    # the other
    if d <= abs(r-R)
        area1 = pi * r^2
        area2 = pi * R^2
        return minimum([area1,area2])
    end

    if R < r
        # swap
        r = B.r
        R = A.r
    end
    part1 = r*r*acos((d*d + r*r - R*R)/(2*d*r))
    part2 = R*R*acos((d*d + R*R - r*r)/(2*d*R))
    part3 = 0.5*sqrt((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R))

    area = part1 + part2 - part3;
    return area
end

function agolFlux(star, planet)
    # star and planet are a Circle struct
    # directly copied from analytic lightcurve paper 2008
    # put things in terms of their paper
    # returns relative flux F
    r_p = planet.r
    d = sqrt((star.x-planet.x)^2 + (star.y - planet.y)^2)
    z = d / star.r
    p = r_p / star.r
    if (1+p) < z
        lambdaE = 0
    elseif abs(1-p) < z & z<= 1+p
        k_1 = acos((1-p*p+z*z)/(2*z))
        K_o = acos((p*p + z*z - 1)/(2*p*z))
        lambdaE = 1/pi*(p*p*k_o+k_1-sqrt((4*z*z-(1+z*z-p*p)^2)/4))
    elseif z<=1-p
        lambdaE = p*p
    elseif z<=p-1
        lambaE = 1
    else
        error("unhandled case in agolFlux!!!")
    return 1 - lambdaE
end


starRad = 1
starArea = pi * starRad^2
planetRad = starRad / 10

star = Circle(starRad, 0, 0)
starY = [0, 0.5, 1]
starX = -2:.001:2
for sy in starY
    overlapA = zeros(length(starX))
    for (ii, x) in enumerate(starX)
        planet = Circle(planetRad, x, sy)
        overlapA[ii] = circOverlap(star, planet)
    end

    light = (starArea .- overlapA) ./ starArea

    plot(starX, light, ylabel="flux", xlabel="star X", alpha=0.5, markersize=1, dpi=150)
    savefig("modeltrans$sy")
end

function linLimb(rStar,  coeffs)
    # x is radial distance from center of star
    # assumed from 0-1
    # c1, c2 < 1
    mu = sqrt.(1 .- rStar.^2)
    I = (1 .- coeffs[1].*(1 .- mu))
    I = I ./ sum(I) # normalize to integrate to 1
    return I
end

function quadLimb(rStar,  coeffs)
    # x is radial distance from center of star
    # assumed from 0-1
    # c1, c2 < 1
    mu = sqrt.(1 .- rStar.^2)
    I = (1 .- coeffs[1].*(1 .- mu) .- coeffs[2].*(1 .- mu).^2)
    I = I ./ sum(I) # normalize to integrate to 1
    return I
end

function sqrtLimb(rStar,  coeffs)
    # x is radial distance from center of star
    # assumed from 0-1
    # c1, c2 < 1
    mu = sqrt.(1 .- rStar.^2)
    I = (1 .- coeffs[1].*(1 .- mu) .- coeffs[2].*(1 .- mu.^2))
    I = I ./ sum(I) # normalize to integrate to 1
    return I
end

function expLimb(rStar,  coeffs)
    # x is radial distance from center of star
    # assumed from 0-1
    # c1, c2 < 1
    mu = sqrt.(1 .- rStar.^2)
    I = (1 .- coeffs[1].*(1 .- mu) .- coeffs[2]./(1 .- exp.(mu)))
    I = I ./ sum(I) # normalize to integrate to 1
    return I
end

limbList = [linLimb, quadLimb, sqrtLimb, expLimb]


nDisks = 1000
limbRs = LinRange(0, starRad, nDisks) #.000001 to avoid zero division error
limbRs = limbRs[2:end] # first point is r=0, so don't use it
limbIs = quadLimb(limbRs, [0.5, 0.5])

# normalize Is to integrate to 1
# limbIs = limbIs ./ sum(limbIs)

# plot the limb dark model
plot(limbRs, limbIs, ylabel="I", xlabel="star radius", dpi=150)
savefig("quadlimb")


for sy in starY
    lightArr = zeros(length(starX))
    for (ii, x) in enumerate(starX)
        planet = Circle(planetRad, x, sy)
        light = 0
        for (starRad, limbI) in zip(limbRs, limbIs)
            star = Circle(starRad, 0, 0)
            overlap = circOverlap(star, planet)
            area = pi * starRad * starRad
            l = (area - overlap) / area * limbI
            light += l

        end
        lightArr[ii] = light
    end

    plot(starX, lightArr, ylabel="flux", xlabel="star X", alpha=0.5, markersize=1, dpi=150)
    savefig("limbmodeltrans$sy")
end


##########################################################################################
##########################


function foldAndSortByTime(modTime, ts, flux)#, fluxErr)
    # return a folded time series time, flux, fluxErr
    ts = ts .% modTime # wrap the
    sortInds = sortperm(ts)
    ts = ts[sortInds]
    flux = flux[sortInds]
    # fluxErr = fluxErr[sortInds]
    return ts, flux #, fluxErr
end

function sumMinDist(flux)
    v1 = flux[1:end-1]
    v2 = flux[2:end]
    dv = v2 .- v1
    # compute array of successive distances
    return sum(dv.^2)
end

function periodSearch(ts, flux, minTime, maxTime)#, fluxErr)
    # a brute force approach
    # returns: bestPeriod, folded time vector, velocity vector
    timeStep = .00001 # days
    foldTimes = minTime:timeStep:maxTime
    minMetric = 1e16
    bestP = -1 # period init
    bestT = ts # times
    bestF = flux # velocities
    for (i, tryP) in enumerate(foldTimes)
        tryT, tryV = foldAndSortByTime(tryP, ts, flux)
        metric = sumMinDist(tryV)
        if metric < minMetric
            minMetric = metric
            bestT = tryT
            bestF = tryV
            # bestE = tryE
            bestP = tryP
        end
    end
    return bestP, bestT, bestF
end

function binData(ts, flux, nBins)
    nPts = length(flux)
    # print("nPts $nPts\n")
    stepSize = Int(floor(nPts / nBins))
    edge1 = 1
    edge2 = stepSize
    binnedF = zeros(nBins)
    binnedT = zeros(nBins)
    for i in 1:nBins
        # print("edges $edge1 $edge2\n")
        binnedF[i] = mean(flux[edge1:edge2])
        binnedT[i] = mean(ts[edge1:edge2])
        edge1 += stepSize
        edge2 += stepSize
        # print("$i\n")
    end
    return binnedT, binnedF
end

function movingAverage(ts, flux, width)
    newFlux = zeros(length(flux)-width)
    newT = zeros(length(newFlux))
    for i in 1:length(newFlux)
        newT[i] = ts[i]
        edge1 = i
        edge2 = i+width
        newFlux[i] = mean(flux[edge1:edge2])
    end
    return newT, newFlux
end



mp = readdlm("mystery_planet02.txt")
ts = mp[:,1]
flux = mp[:,2]
fluxErr = mp[:,3]

# ts, flux = movingAverage(ts, flux, 50)
# ts, flux = movingAverage(ts, flux, 50)
# ts, flux = smoothData(ts, flux, 500)

plot(ts, flux, ylabel="flux", xlabel="time", seriestype=:scatter, alpha=0.5, markersize=1, dpi=150)
savefig("raw")

tsSmooth, fluxSmooth = movingAverage(ts, flux, 100)
plot(tsSmooth, fluxSmooth, ylabel="flux", xlabel="time", seriestype=:scatter, alpha=0.5, markersize=1, dpi=150)
savefig("rawSmoothed")

# tsBin, fluxBin = binData(ts, flux, 1000)
# plot(tsBin, fluxBin, ylabel="flux", xlabel="time", seriestype=:scatter, alpha=0.5, markersize=1, dpi=150)
# savefig("rawBinned")


period, tsFold, fluxFold = periodSearch(ts, flux, 9.5, 10.5)#, fluxErr)

# period = period * 2.0

print("best period $period days\n")

#9.799999999999999

# period = 9.811
tsFold, fluxFold = foldAndSortByTime(period, ts, flux)#, fluxErr)

# plot the folded timeseries
plot(tsFold, fluxFold, ylabel="flux (%?)", xlabel="time folded on period=$period", seriestype=:scatter, alpha=0.5, markersize=1, dpi=150)
savefig("periodFolded")

tsSmooth, fluxSmooth = movingAverage(tsFold, fluxFold, 25)

plot(tsSmooth, fluxSmooth, minorgrid=true, xminorticks=0:.25:10, xticks=0:1:10, ylabel="flux", xlabel="time folded on period=$period", seriestype=:scatter, alpha=0.5, markersize=1, dpi=150)
savefig("periodFoldedSmooth")

# tsBin, fluxBin = smoothData(ts, flux, 500)

# plot(tsBin, fluxBin, ylabel="flux (%?)", xlabel="days folded on period=$period", seriestype=:scatter, alpha=0.5, markersize=1, dpi=150)
# savefig("periodFoldedBin")



################ transit modeling ###################################################
# find period, depth, impact parameter, duration, density of the star
function transitModel(tVec, P, t_o, T, b, r_p, limbModel, coeffs)
    # tVec = vector of times at which to compute transit
    # P = period
    # t_o = midpoint in time of transit, from time=0 (modulo period!)
    # T = transit duration (1st to 4th contact)
    # b = impact parameter (fraction of star's radius)
    # r_p = planet's radius (as a fraction of star's radius)
    # limbModel = an integer from 0-4 indicating which limb profile to use
    #              0 = uniform
    #              1 = linear (coeff2 is ignored)
    #              2 = quadratic
    #              3 = sqrt
    #              4 = exponential
    # coeffs = list of coefficients to pass to limb darkening models (ignored for
    #           limbModel = 0, the uniform case)

    # returns vector of modeled flux values normalized to 1

    # first wrap tVec on period P
    tVec = tVec .% P

    starRad = 1 # radius of star is always taken to be 1
    nCakeLayers = 1000
    if limbModel > 0
        limbRs = LinRange(0, starRad, nCakeLayers)
        limbRs = limbRs[2:end] # ignore r=0, first element
        limbIs = limbList[limbModel](limbRs, coeffs) # these are normalized to sum to 1
    else
        # a uniform star, no limb darkening
        limbRs = [1]
        limbIs = [1]
    end

    # based on T, duration of transit, determine the scaling between
    # position of the planet and time
    xContact = sqrt((starRad+r_p)^2 - (starRad*b)^2)
    xLength = 2*xContact
    # turn this in to a velocity dx/dt
    v = xLength / T

    # rescale all points in tVec to x, and center at t_o (center of transit at x = 0)
    # now x vec is relative to star x center for a single transit
    xVec = (tVec .- t_o) .* v

    fluxVec = zeros(length(xVec))
    for (ii, starX) in enumerate(xVec)
        planet = Circle(r_p, starX, b)
        limbSum = 0
        for (rStar, limbI) in zip(limbRs, limbIs)
            # rStar is increasing cake radii
            star = Circle(rStar, 0, 0)
            overlap = circOverlap(star, planet)
            area = pi * starRad * starRad
            flux = (area - overlap) / area * limbI
            limbSum += flux
        end
        fluxVec[ii] = limbSum
    end

    return fluxVec
end

### next come up with some good guesses for our model fit
# set t_o to be the median value in the smoothed folded data
# for values < 0.998
transitInds = fluxSmooth .< 0.998
tZoom = tsSmooth[transitInds]
t_o = median(tZoom)

# initial guess for depth of transit, take the upper quartile
fZoom = fluxSmooth[transitInds]
depth = quantile(fZoom, 0.75)
# convert depth into initial guess for planet radius r_p (fraction of star radius=1)
rStar = 1
r_p = sqrt(rStar^2*(1 - depth))
print("transit depth $depth\n")
print("r_p $r_p\n")

# estimate transit duration from visual inspection of plot
T = 0.4 # total transit duration
tao = 0.1 # ingress/egress duration

# finally estimate the impact parameter
b_o = sqrt(1 - r_p*T/tao)

tModel = transitModel(tsSmooth, period, t_o, T, b_o, r_p, 0, nothing)


plot(tsSmooth, fluxSmooth, minorgrid=true, xminorticks=0:.25:10, xticks=0:1:10, ylabel="flux", xlabel="time folded on period=$period", seriestype=:scatter, alpha=0.5, markersize=1, dpi=150)
plot!(tsSmooth, tModel)
savefig("model-smooth")

# seems somewhat reasonable lets see how well it fits the real data

tModel = transitModel(tsFold, period, t_o, T, b_o, r_p, 0, nothing)

chi2 = sum((tModel .- flux).^2)

print("chi2 first guess: $chi2\n")

# plot the folded data on period
plot(tsFold, fluxFold, minorgrid=true, xminorticks=0:.25:10, xticks=0:1:10, ylabel="flux", xlabel="time folded on period=$period", seriestype=:scatter, alpha=0.5, markersize=1, dpi=150)
plot!(tsFold, tModel)
savefig("model")

# define a struct for keeping track of params
struct Params
    period
    t_o
    T
    b
    r_p
    limbModel
    limbModelCoeff
    chi2
end

periodRange = LinRange(period-0.2, period+0.2, 100)

# t_oRange = LinRange(t_o-0.2, t_o+0.2, 50)
# TRange = LinRange(0.2, 0.6, 50)
# b_oRange = LinRange(b_o/10, 1, 50)
# r_pRange = LinRange(r_p/10, 0.7, 50)

# after running once:
t_oRange = LinRange(5.7, 5.8, 15)
TRange = LinRange(0.38, 0.41, 10)
b_oRange = LinRange(0.8, 1, 25)
r_pRange = LinRange(0.03, 0.1, 25)


bestParams = Params(period,t_o,T,b_o,r_p,0,0,chi2)

print("start params $bestParams\n")

if true # do initial grid search
    for to in t_oRange
        print("---------to $to\n")
        for T in TRange
            print("    ---------T $T\n")
            for b in b_oRange
                for rp in r_pRange
                    tModel = transitModel(ts, period, to, T, b, rp, 0, nothing)
                    chi2 = sum((tModel .- flux).^2)
                    if chi2 < bestParams.chi2
                        global bestParams
                        bestParams = Params(period, to, T, b, rp, 0, nothing, chi2)
                        print("newFit $bestParams\n")
                    end
                end
            end
        end
    end
else
    # best fit (from grid search)
    bestParams = Params(9.80797, 5.75, 0.38, 0.9083333333333333, 0.05625000000000001, 0, nothing, 0.03194680026858922)
    tModel = transitModel(tsFold, bestParams.period, bestParams.t_o, bestParams.T, bestParams.b, bestParams.r_p, 0, nothing)
    chi2 = sum((tModel .- flux).^2)
    bestParams = Params(bestParams.period, bestParams.t_o, bestParams.T, bestParams.b, bestParams.r_p, 0, nothing, chi2)
    print("bestParams in if block $bestParams\n")

end

# plot the folded data on period
plot(tsFold, fluxFold, minorgrid=true, xminorticks=0:.25:10, xticks=0:1:10, ylabel="flux", xlabel="time folded on period=$period", seriestype=:scatter, alpha=0.5, markersize=1, dpi=250)
plot!(tsFold, tModel)
savefig("modelSearchedFolded")


tModel = transitModel(ts, bestParams.period, bestParams.t_o, bestParams.T, bestParams.b, bestParams.r_p, 0, nothing)
plot(ts, flux, minorgrid=true, ylabel="flux", xlabel="time", seriestype=:scatter, alpha=0.5, markersize=1, dpi=250)
plot!(ts, tModel)
savefig("modelSearched")

#### next grid for limb darkening
bestParams1 = deepcopy(bestParams)

# finially lets let the minimizer try to find an even better soln, with limb darkening this time!
# use a quadratic limb profile

function minimizeModel(x)
    period = abs(x[1])
    t_o = abs(x[2])
    T = abs(x[3])
    b = abs(x[4])
    if b > 1
        b = 1
    end
    r_p = abs(x[5])
    if r_p > 1
        r_p = 1
    end
    c1 = abs(x[6])
    if c1 > 1
        c1 = 1
    end
    c2 = abs(x[7])
    if c2 > 1
        c2 = 1
    end
    tModel = transitModel(ts, period, t_o, T, b, r_p, 2, [c1, c2])
    chi2 = sum((tModel .- flux).^2)
    return chi2
end

x_o = [bestParams.period, bestParams.t_o, bestParams.T, bestParams.b, bestParams.r_p, 0, 0]

print("begin minimizer\n")

out = optimize(minimizeModel, x_o)
print("minimizer output:\n $out\n")
bestFit = Optim.minimum(out)
print("best fit $bestFit\n")
bestFit = Optim.minimizer(out)
print("best fit $bestFit\n")

# fix t_o, that should be well constrained
# fix r_p that should be well constrained
# lessen T range
# allow b to vary
# print("par1 $bestParams\n")

# using Profile
# @profile tModel = transitModel(ts, bestParams1.period, bestParams1.t_o, bestParams1.T, bestParams1.b, bestParams1.r_p, 2, [0.1, 0.45])
# Profile.print()


# c1s = LinRange(0,1,25)
# c2s = LinRange(0,1,25)
# TRange = LinRange(.3, .5, 25)
# # keep same b_o range as previous
# r_pRange = LinRange(0.04, 0.07, 25)
# for c1 in c1s
#     for c2 in c2s
#         tModel = transitModel(ts, bestParams1.period, bestParams1.t_o, bestParams1.T, bestParams1.b, bestParams1.r_p, 2, [c1, c2])
#         chi2 = sum((tModel .- flux).^2)
#         print("$c1 $c2 $chi2\n")
#         if chi2 < bestParams.chi2
#             global bestParams
#             bestParams = Params(bestParams1.period, bestParams1.t_o, bestParams1.T, bestParams1.b, bestParams1.r_p, 2, [c1, c2], chi2)
#             print("newFit+limbd $bestParams\n")
#         end

#     end
# end









########### FFT stuff
function doFFT(ts, flux)
    subMean = flux .- mean(flux)
    ft = fft(subMean)
    pwr = real.(ft .* conj.(ft))
    nBins = length(pwr)
    dt = mean(ts[2:end] - ts[1:end-1])
    dt_sig = std(ts[2:end] - ts[1:end-1])
    print("dt $dt with std $dt_sig\n")

    freq = float.(1:nBins) ./ (nBins .* dt)

    # only keep the first half (to nyquest)
    # look in a reasonable range of frequencies
    pwr = pwr[1:Int(floor(nBins / 2))]
    freq = freq[1:Int(floor(nBins / 2))]
    freq1 = freq[1]
    freq2 = freq[end]
    print("freqs $freq1 $freq2\n")

    pwr = pwr[1:10]
    freq = freq[1:10]

    maxVal, maxInd = findmax(pwr)
    bestFreq = freq[maxInd]
    period = 1 / bestFreq

    print("best period $period at freq $bestFreq\n")

    plot(freq, pwr, ylabel="fft", xlabel="frequency", alpha=0.5, markersize=1, dpi=150)
    savefig("fft")
    # tryP = 37.1275
    # period = 37.1275

    # period = 9.811 #looks pretty good too
    # tryT, tryV, tryE = foldAndSortByTime(tryP, ts, flux, fluxErr)
    return period
end

# period = doFFT(tsSmooth, fluxSmooth)

# # sort by increasing time
# # set maximum time to ensure no fold yet
# # not all points are in order from file
# ts1, flux1 = foldAndSortByTime(period, ts, flux)
# # plot the folded timeseries
# plot(ts1, flux1, ylabel="flux (%?)", xlabel="days folded on period=$period", seriestype=:scatter, alpha=0.5, markersize=1, dpi=150)
# savefig("periodFoldedFFT")





