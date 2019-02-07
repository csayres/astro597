using Plots
using DelimitedFiles

function foldAndSortByTime(modTime, ts, radVel, vErr)
    # return a folded time series time, radVel, vErr
    ts = ts .% modTime # wrap the
    sortInds = sortperm(ts)
    ts = ts[sortInds]
    radVel = radVel[sortInds]
    vErr = vErr[sortInds]
    return ts, radVel, vErr
end

function sumMinDist(radVel, vErr)
    v1 = radVel[1:end-1]
    v2 = radVel[2:end]
    e1 = vErr[1:end-1]
    e2 = vErr[2:end]
    sigmas = e1.^2 .+ e2.^2
    dv = v2 .- v1
    # compute array of successive distances
    # scaled by the errors
    metric = sum(dv.^2 ./ sigmas) # minimize this
    return metric
end

function periodSearch(ts, radVel, vErr)
    # a brute force approach
    # returns: bestPeriod, folded time vector, velocity vector, error vector
    timeStep = .001 # days
    minTime = 0
    maxTime = 1000 # by visual inspection period is obviously less than this
    foldTimes = minTime:timeStep:maxTime
    minMetric = sumMinDist(radVel, vErr)
    bestP = -1 # period init
    bestT = ts # times
    bestV = radVel # velocities
    bestE = vErr # errors
    for (i, tryP) in enumerate(foldTimes)
        tryT, tryV, tryE = foldAndSortByTime(tryP, ts, radVel, vErr)
        metric = sumMinDist(tryV, tryE)
        # metrics[i] = metric
        if metric < minMetric
            minMetric = metric
            bestT = tryT
            bestV = tryV
            bestE = tryE
            bestP = tryP
        end
    end
    return bestP, bestT, bestV, bestE
end

mp = readdlm("mystery_planet1.txt")
ts = mp[:,1]
flux = mp[:,2]
fluxError = mp[:,3]
tZero = ts[1]
ts = ts .- tZero # normalize time to begin at zero

maxErr = maximum(fluxErr)
minErr = minimum(fluxErr)
print("errs $maxErr $minErr\n")

# sort by increasing time
# set maximum time to ensure no fold yet
# not all points are in order from file
ts, radVel, vErr = foldAndSortByTime(maximum(ts)+100, ts, radVel, vErr)
period, tsFold, radVelFold, vErrFold = periodSearch(ts, radVel, vErr)

print("best period $period days\n")

# plot the folded timeseries
plot(tsFold, radVelFold, ylabel="rv (km/s)", xlabel="days folded on period=$period", markersize=1, markershape=:auto, line=:solid, linealpha=0.5, yerror=vErrFold, dpi=150)
savefig("periodFolded")
# plotTimeseries(tsFold, radVelFold, vErrFold, "periodFolded")


function keplerSolve(e,M)
    # returns a E for a given e, M
    # ensure M is wrapped between 0 and 2*pi
    while M < 0
        M += 2*pi
    end
    while M > 2*pi
        M -= 2*pi
    end
    E_o = M + 0.85 * e * sign(sin(M))
    # exit tolerance
    epsilon = 1e-12
    maxIter = 1000
    iter = 1

    while iter < maxIter
        global E_next # eventual return value
        global E_error
        iter += 1
        der = 1 - e * cos(E_o)
        if der == 0
            der = 1e-16 # avoid zero division
        end
        E_next = E_o - (E_o - e * sin(E_o) - M) / der
        E_error = abs(E_o-E_next)
        if E_error < epsilon
            break # stop iteration we've converged
        end
        # update value and continue looping
        E_o = E_next
    end

    return E_o
end

# make these containers outside the function for reusability
Es = zeros(length(ts))
A = zeros(3,3)
b = zeros(3)

# inverse sigma squred is used allover the place so
# just compute it here
invSigma2 = vErr .^ -2

function solveOrbit(e, period, t_p)
    # returns a chi2, h, c, and v_o for the given inputs
    # Ms = 2*pi/period .* (ts .- t_p)
    # Ms = Ms .% (2*pi)

    Ms = 2*pi/period .* (ts .- t_p)

    # # wrap Ms to 0-2*pi
    # for ind in findall(Ms.<0)
    #     print("updating index $ind\n")
    #     Ms[ind] = Ms[ind] + 2*pi
    # end


    # use keplers eqn solver to get the Es
    for (ii, M) in enumerate(Ms)
        # ensure M is wrapped correctly between 0-2pi
        ### cut and past in slover here to increase speed?
        Es[ii] = keplerSolve(e,M)
    end
    fs = 2 * atan.( ((1+e)/(1-e))^0.5 * tan.(Es/2) )
    # linear solver for h, c and v_0
    # precompute some stuff
    cosfs = cos.(fs) # compute once
    sinfs = sin.(fs) # compute once
    sincosfs = cosfs .* sinfs
    cos2fs = cosfs .^ 2
    sin2fs = sinfs .^ 2

    A[1,1] = sum(invSigma2 .* cos2fs )
    A[1,2] = sum(invSigma2 .* sincosfs)
    A[1,3] = sum(invSigma2 .* cosfs)

    A[2,1] = A[1,2]
    A[2,2] = sum(invSigma2 .* sin2fs)
    A[2,3] = sum(invSigma2 .* sinfs)

    A[3,1] = A[1,3]
    A[3,2] = A[2,3]
    A[3,3] = sum(invSigma2)

    b[1] = sum(invSigma2 .* radVel .* cosfs)
    b[2] = sum(invSigma2 .* radVel .* sinfs)
    b[3] = sum(invSigma2 .* radVel)

    h,c,v_o = A\b

    vModel = h .* cosfs .+ c .* sinfs .+ v_o

    chi2 = sum( (radVel .- vModel).^2 .* invSigma2 )
    return chi2, h, c, v_o
end

# brute force solve for e, and tp
eRange = 0.7:.01:0.99
tpRange = 0:0.1:period
bestChi2 = 1e16 # initialize to a big number
best_e = 0
best_tp = 0
best_h = 0
best_c = 0
best_vo = 0
for e in eRange
    for t_p in tpRange
        global bestChi2
        chi2, h, c, v_o = solveOrbit(e, period, t_p)
        if chi2 < bestChi2
            global best_e
            global best_tp
            global best_h
            global best_c
            global best_vo
            bestChi2 = chi2
            best_e = e
            best_tp = t_p
            best_h = h
            best_c = c
            best_vo = v_o
        end
    end
end

print("brute force e=$best_e and t_p=$best_tp\n")

# convert h, c, v_o to gamma, K, omega
gamma = best_vo - best_e * best_h # m/s
K = sqrt(best_h^2+best_c^2) # m/s
omega = acos(best_h / K) * 360 / (2*pi) # degrees
print("period = $period days\ne = $best_e\nt_p = $best_tp days after tZero = $tZero (first point in timeseries)\ngamma = $gamma km/s\nK = $K km/s\nomega = $omega degrees\n")


# plot answers
function plotOrbitSoln(e, P, t_p, h, c, v_o)
    tMax = maximum(tsFold)
    tMin = minimum(tsFold)
    tPts = tMin:0.01:tMax
    Ms = 2*pi/P .* (tPts .- t_p)

    # Ms = Ms .% (2*pi)
    Es = zeros(length(tPts))
    # use keplers eqn solver to get the Es
    for (ii, M) in enumerate(Ms)
        Es[ii] = keplerSolve(e,M)
    end
    fs = 2 * atan.( ((1+e)/(1-e))^0.5 * tan.(Es/2) )
    vModel = h .* cos.(fs) .+ c .* sin.(fs) .+ v_o
    plot(tsFold, radVelFold, ylabel="rv (km/s)", xlabel="days folded on period=$P", markersize=2, seriestype=:scatter, yerror=vErrFold, label="Folded Data", dpi=150)
    plot!(tPts, vModel, line=:solid, label="Model")
    savefig("model")
end

plotOrbitSoln(best_e, period, best_tp, best_h, best_c, best_vo)


