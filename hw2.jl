using Plots
using DelimitedFiles
using Statistics
using FFTW

function foldAndSortByTime(modTime, ts, flux)#, fluxErr)
    # return a folded time series time, flux, fluxErr
    ts = ts .% modTime # wrap the
    sortInds = sortperm(ts)
    ts = ts[sortInds]
    flux = flux[sortInds]
    # fluxErr = fluxErr[sortInds]
    return ts, flux #, fluxErr
end

function sumMinDist(flux) #, fluxErr)

    # nPts = length(flux)
    # # print("nPts $nPts\n")
    # nBins = 200
    # stepSize = Int(floor(nPts / nBins))
    # edge1 = 1
    # edge2 = stepSize
    # binned = zeros(nBins)
    # for i in 1:nBins
    #     # print("edges $edge1 $edge2\n")
    #     binned[i] = mean(flux[edge1:edge2])
    #     edge1 += stepSize
    #     edge2 += stepSize
    #     # print("$i\n")
    # end


    v1 = flux[1:end-1]
    v2 = flux[2:end]

    e1 = fluxErr[1:end-1]
    e2 = fluxErr[2:end]
    sigmas = e1.^2 .+ e2.^2
    dv = v2 .- v1
    # compute array of successive distances
    # scaled by the errors
    metric = sum(dv.^2) # ./ sigmas) # minimize this
    return metric
end

function periodSearch(ts, flux)#, fluxErr)
    # a brute force approach
    # returns: bestPeriod, folded time vector, velocity vector, error vector
    timeStep = .001 # days
    minTime = 7
    maxTime = 12 # by visual inspection period is obviously less than this
    foldTimes = minTime:timeStep:maxTime
    # minMetric = sumMinDist(flux, fluxErr)
    minMetric = 1e16
    bestP = -1 # period init
    bestT = ts # times
    bestF = flux # velocities
    # bestE = fluxErr # errors
    for (i, tryP) in enumerate(foldTimes)
        tryT, tryV = foldAndSortByTime(tryP, ts, flux)#, fluxErr)
        metric = sumMinDist(tryV) #, tryE)
        # metrics[i] = metric
        if metric < minMetric
            minMetric = metric
            bestT = tryT
            bestF = tryV
            # bestE = tryE
            bestP = tryP
        end
    end
    return bestP, bestT, bestF #, bestE
end

function smoothData(ts, flux, nBins)
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

plot(ts, flux, ylabel="flux (%?)", xlabel="smoothed", seriestype=:scatter, alpha=0.5, markersize=1, dpi=150)
savefig("smoothed")

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

pwr = pwr[1:50]
freq = freq[1:50]

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

# sort by increasing time
# set maximum time to ensure no fold yet
# not all points are in order from file
ts1, flux1 = foldAndSortByTime(period, ts, flux)
# plot the folded timeseries
plot(ts1, flux1, ylabel="flux (%?)", xlabel="days folded on period=$period", seriestype=:scatter, alpha=0.5, markersize=1, dpi=150)
savefig("periodFoldedFFT")

period, tsFold, fluxFold = periodSearch(ts, flux)#, fluxErr)

print("best period $period days\n")

# period = 9.811
tsFold, fluxFold = foldAndSortByTime(period, ts, flux)#, fluxErr)

# plot the folded timeseries
plot(tsFold, fluxFold, ylabel="flux (%?)", xlabel="days folded on period=$period", seriestype=:scatter, alpha=0.5, markersize=1, dpi=150)
savefig("periodFolded")

tsSmooth, fluxSmooth = movingAverage(ts, flux, 500)

plot(tsSmooth, fluxSmooth, ylabel="flux (%?)", xlabel="days folded on period=$period", seriestype=:scatter, alpha=0.5, markersize=1, dpi=150)
savefig("periodFoldedSmooth")
