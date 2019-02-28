using Plots
using DelimitedFiles
using Statistics
using FFTW
using Optim

struct Circle
    r::Float64
    x::Float64
    y::Float64
end


function batman(star, planet)
    # inputs: star, planet Circle structures
    # returns fraction of flux blocked from planet
    # basically copied from the batman paper on astro-ph

    d = sqrt((star.x-planet.x)^2 + (star.y - planet.y)^2)
    x = star.r
    rp = planet.r
    global A
    if (rp - d < x) & (rp + d > x)
        u = (d*d + x*x - rp*rp)/(2*d*x)
        u = min(u, 1) # because computers aren't perfect
        v = (d*d + rp*rp - x*x)/(2*d*rp)
        v = min(v, 1) # because computers aren't perfect
        w = (-d + x + rp)*(d + x - rp)*(d - x + rp)*(d + x + rp)
        w = max(w, 0) # because computers aren't perfect
        A = x*x *acos(u) + rp*rp*acos(v) - 0.5*sqrt(w)
    elseif x<= rp - d
        A = pi*x*x
    elseif x>= rp + d
        A = pi*rp*rp
    end
    return A
end


starRad = 1;
planetRad = starRad / 10;

star = Circle(starRad, 0, 0);
starY = [0, 0.5, 1];
starX = -2:.001:2;

for sy in starY
    flux = zeros(length(starX))
    for (ii, x) in enumerate(starX)
        planet = Circle(planetRad, x, sy)
        flux[ii] = 1 - batman(star, planet)
    end
    if sy == 0
        plot(starX, flux, plot_title="r_p = 0.1*r_star, no limb darkening", ylabel="relative flux", xlabel="time", label="b = $sy", alpha=0.5, markersize=1, dpi=250)
    else
        plot!(starX, flux, label="b = $sy", alpha=0.5, markersize=1)
    end
end
savefig("modeltrans")

function linLimb(rStar,  coeffs)
    # x is radial distance from center of star
    # assumed from 0-1
    # c1, c2 < 1
    mu = sqrt.(1 .- rStar.^2)
    I = (1 .- coeffs[1].*(1 .- mu))
    return I
end

function quadLimb(rStar,  coeffs)
    # x is radial distance from center of star
    # assumed from 0-1
    # c1, c2 < 1
    mu = sqrt.(1 .- rStar.^2)
    I = (1 .- coeffs[1].*(1 .- mu) .- coeffs[2].*(1 .- mu).^2)
    return I
end

function sqrtLimb(rStar,  coeffs)
    # x is radial distance from center of star
    # assumed from 0-1
    # c1, c2 < 1
    mu = sqrt.(1 .- rStar.^2)
    I = (1 .- coeffs[1].*(1 .- mu) .- coeffs[2].*(1 .- mu.^2))
    return I
end

function expLimb(rStar,  coeffs)
    # x is radial distance from center of star
    # assumed from 0-1
    # c1, c2 < 1
    mu = sqrt.(1 .- rStar.^2)
    I = (1 .- coeffs[1].*(1 .- mu) .- coeffs[2]./(1 .- exp.(mu)))
    return I
end

limbList = [linLimb, quadLimb, sqrtLimb, expLimb];


nDisks = 1000;
limbRs = LinRange(0, starRad, nDisks);
limbRs = limbRs[2:end];
midPointRs = (limbRs[2:end] .+ limbRs[1:end-1]) ./ 2;
limbIs = quadLimb(midPointRs, [0.5, 0.5]);

# plot the limb dark model
plot(midPointRs, limbIs, label="quadratic profile c1=0.5, c2=0.5", ylabel="I(r)", xlabel="star radius", dpi=250)
savefig("quadlimb")


for sy in starY
    flux = zeros(length(starX))
    for (ii, x) in enumerate(starX)
        planet = Circle(planetRad, x, sy)
        shellFluxes = zeros(length(limbRs))
        for (i,starRad) in enumerate(limbRs)
            star = Circle(starRad, 0, 0)
            flux_i = batman(star, planet)
            shellFluxes[i] = flux_i
        end
        dShell = shellFluxes[2:end] - shellFluxes[1:end-1]
        flux[ii] = 1 - sum(limbIs .* dShell)
    end

    if sy == 0
        plot(starX, flux, plot_title="r_p = 0.1*r_star, quadratic limb darkening c1=0.5 c2=0.5", ylabel="relative flux", xlabel="time", label="b = $sy", alpha=0.5, markersize=1, dpi=250)
    else
        plot!(starX, flux, label="b = $sy", alpha=0.5, markersize=1)
    end

    # plot(starX, lightArr, ylabel="flux", xlabel="star X", alpha=0.5, markersize=1, dpi=250)
end
savefig("limbmodeltrans")


##########################################################################################
##########################


function foldAndSortByTime(modTime, ts, flux)
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

plot(ts, flux, label="raw time series", ylabel="relative flux", xlabel="time", seriestype=:scatter, alpha=0.5, markersize=1, dpi=250)
savefig("raw")

tsSmooth, fluxSmooth = movingAverage(ts, flux, 20)
plot(tsSmooth, fluxSmooth, label="smoothed time series", ylabel="relative flux", xlabel="time", seriestype=:scatter, alpha=0.5, markersize=1, dpi=250)
savefig("rawSmoothed")

period, tsFold, fluxFold = periodSearch(ts, flux, 9.5, 10.5)

print("best period $period days\n")

tsFold, fluxFold = foldAndSortByTime(period, ts, flux)#, fluxErr)

# plot the folded timeseries
plot(tsFold, fluxFold, ylabel="relative flux", xlabel="time folded on period=$period", seriestype=:scatter, alpha=0.5, markersize=1, dpi=250)
savefig("periodFolded")

tsSmooth, fluxSmooth = movingAverage(tsFold, fluxFold, 25)

plot(tsSmooth, fluxSmooth, minorgrid=true, xminorticks=0:.25:10, xticks=0:1:10, ylabel="smoothed relative flux", xlabel="time folded on period=$period", seriestype=:scatter, alpha=0.5, markersize=1, dpi=250)
savefig("periodFoldedSmooth")


################ transit modeling ###################################################
function transitModel(tVec, P, t_o, T, b, r_p, limbModel, coeffs)
    # tVec = vector of times at which to compute transit
    # P = period
    # t_o = midpoint in time of transit, from time=0 (modulo period!)
    # T = transit duration (1st to 4th contact)
    # b = impact parameter (fraction of star's radius)
    # r_p = planet's radius (as a fraction of star's radius)
    # limbModel = an integer from 0-4 indicating which limb profile to use
    #              0 = uniform (no limb darkening)
    #              1 = linear (coeff2 is ignored)
    #              2 = quadratic
    #              3 = sqrt
    #              4 = exponential
    # coeffs = list of coefficients to pass to limb darkening models (ignored for
    #           limbModel = 0, the uniform case)

    # returns vector of modeled flux values relative to 1

    # first wrap tVec on period P
    tVec = tVec .% P

    starRad = 1 # radius of star is always taken to be 1


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
    for (ii, planetX) in enumerate(xVec)
        planet = Circle(r_p, planetX, b)
        if limbModel == 0
            # no limb darkening
            star = Circle(1, 0, 0)
            fluxVec[ii] = 1 - batman(star, planet)
        else
            # yes limb darkening
            limbIs = limbList[limbModel](midPointRs, coeffs)
            shellFluxes = zeros(length(limbRs))
            for (i,starRad) in enumerate(limbRs) # avoid r=0
                star = Circle(starRad, 0, 0)
                shellFluxes[i] = batman(star, planet)
            end
            dShell = shellFluxes[2:end] - shellFluxes[1:end-1]
            fluxVec[ii] = 1 - sum(limbIs .* dShell)
        end
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

# tModel = transitModel(tsSmooth, period, t_o, T, b_o, r_p, 0, nothing)

# seems somewhat reasonable lets see how well it fits the real data

tModel = transitModel(tsFold, period, t_o, T, b_o, r_p, 0, nothing)

chi2 = sum((tModel .- flux).^2)

print("chi2 first guess: $chi2\n")

# plot the folded data on period
plot(tsFold, fluxFold, minorgrid=true, xminorticks=0:.25:10, xticks=0:1:10, ylabel="flux", xlabel="time folded on period=$period", seriestype=:scatter, alpha=0.5, markersize=1, dpi=250)
plot!(tsFold, tModel)
savefig("modelGuess")

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

# periodRange = LinRange(period-0.2, period+0.2, 100)

# t_oRange = LinRange(t_o-0.2, t_o+0.2, 50)
# TRange = LinRange(0.2, 0.6, 50)
# b_oRange = LinRange(b_o/10, 1, 50)
# r_pRange = LinRange(r_p/10, 0.7, 50)

# after running once:
t_oRange = LinRange(5.7, 5.8, 15)
TRange = LinRange(0.38, 0.41, 10)
b_oRange = LinRange(0.8, 1, 25)
r_pRange = LinRange(0.03, 0.1, 25)


bestParams = Params(period,t_o,T,b_o,r_p,0,0,chi2) # initial guess

print("start params $bestParams\n")

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

# plot the folded data on period
tModel = transitModel(tsFold, bestParams.period, bestParams.t_o, bestParams.T, bestParams.b, bestParams.r_p, 0, nothing)
plot(tsFold, fluxFold, minorgrid=true, xminorticks=0:.25:10, xticks=0:1:10, ylabel="flux", xlabel="time folded on period=$period", seriestype=:scatter, alpha=0.5, markersize=1, dpi=250)
plot!(tsFold, tModel)
savefig("modelSearchedFolded")


tModel = transitModel(ts, bestParams.period, bestParams.t_o, bestParams.T, bestParams.b, bestParams.r_p, 0, nothing)
plot(ts, flux, minorgrid=true, ylabel="flux", xlabel="time", seriestype=:scatter, alpha=0.5, markersize=1, dpi=250)
plot!(ts, tModel)
savefig("modelSearched")

#### next grid for limb darkening
# bestParams1 = deepcopy(bestParams)

# finially lets let the minimizer try to find an even better soln, with limb darkening this time!
# use a quadratic limb profile

limbModel = 2 # corresponds to quadratic

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
    tModel = transitModel(ts, period, t_o, T, b, r_p, limbModel, [c1, c2])
    chi2 = sum((tModel .- flux).^2)
    return chi2
end

x_o = [bestParams.period, bestParams.t_o, bestParams.T, bestParams.b, bestParams.r_p, 0, 0]

print("begin minimizer\n")

out = optimize(minimizeModel, x_o)
print("minimizer output:\n $out\n")
bestChi2 = Optim.minimum(out)
print("best fit $bestFit\n")
bestFit = Optim.minimizer(out)
print("best fit $bestFit\n")

best_period = bestFit[1]
best_t_o = bestFit[2]
best_T = bestFit[3]
best_b = bestFit[4]
best_r_p = bestFit[5]
c1 = abs(bestFit[6])
c2 = abs(bestFit[7])

bestParams = Params(best_period, best_t_o, best_T, best_b, best_r_p, bestChi2)

# refold on new period
tsFold, fluxFold = foldAndSortByTime(best_period, ts, flux)

tModel = transitModel(ts, bestParams.period, bestParams.t_o, bestParams.T, bestParams.b, bestParams.r_p, limbModel, [c1, c2])
plot(tsFold, fluxFold, minorgrid=true, xminorticks=0:.25:10, xticks=0:1:10, ylabel="flux", xlabel="time folded on period=$period", seriestype=:scatter, alpha=0.5, markersize=1, dpi=250)
plot!(tsFold, tModel)
savefig("finalModelFolded")


plot(ts, flux, minorgrid=true, ylabel="flux", xlabel="time", seriestype=:scatter, alpha=0.5, markersize=1, dpi=250)
plot!(ts, tModel)
savefig("finalModel")

# quadritic profile
# plot the limb dark model
limbIs = quadLimb(midPointRs, [c1, c2])
plot(midPointRs, limbIs, ylabel="I", xlabel="star radius", dpi=250)
savefig("finalQuadlimb")


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

# function doFFT(ts, flux)
#     subMean = flux .- mean(flux)
#     ft = fft(subMean)
#     pwr = real.(ft .* conj.(ft))
#     nBins = length(pwr)
#     dt = mean(ts[2:end] - ts[1:end-1])
#     dt_sig = std(ts[2:end] - ts[1:end-1])
#     print("dt $dt with std $dt_sig\n")

#     freq = float.(1:nBins) ./ (nBins .* dt)

#     # only keep the first half (to nyquest)
#     # look in a reasonable range of frequencies
#     pwr = pwr[1:Int(floor(nBins / 2))]
#     freq = freq[1:Int(floor(nBins / 2))]
#     freq1 = freq[1]
#     freq2 = freq[end]
#     print("freqs $freq1 $freq2\n")

#     pwr = pwr[1:10]
#     freq = freq[1:10]

#     maxVal, maxInd = findmax(pwr)
#     bestFreq = freq[maxInd]
#     period = 1 / bestFreq

#     print("best period $period at freq $bestFreq\n")

#     plot(freq, pwr, ylabel="fft", xlabel="frequency", alpha=0.5, markersize=1, dpi=250)
#     savefig("fft")
#     # tryP = 37.1275
#     # period = 37.1275

#     # period = 9.811 #looks pretty good too
#     # tryT, tryV, tryE = foldAndSortByTime(tryP, ts, flux, fluxErr)
#     return period
# end

# period = doFFT(tsSmooth, fluxSmooth)

# # sort by increasing time
# # set maximum time to ensure no fold yet
# # not all points are in order from file
# ts1, flux1 = foldAndSortByTime(period, ts, flux)
# # plot the folded timeseries
# plot(ts1, flux1, ylabel="flux (%?)", xlabel="days folded on period=$period", seriestype=:scatter, alpha=0.5, markersize=1, dpi=250)
# savefig("periodFoldedFFT")





