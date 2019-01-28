# using PyCall
# # required for my mac os to work
# @pyimport matplotlib
# matplotlib.use("TkAgg")

# using PyPlot

using Plots
using DelimitedFiles

# M = E - e * sin E
# M = 2*pi / period (t - t_o)

# g(E) --> 0 = E - e * sin E - M
function g(E::Float64, e::Float64, M::Float64)::Float64
    # E = Eccentric anomaly input in radians
    # e = eccentricity (0-1)
    # M = mean anomaly, input in radians

    # returns g as provided in HW
    return E - e * sin(E) - M

end

# time derivative of g(E) for use in newton solver
function dg_dt(E::Float64, e::Float64)::Float64
    # E = Eccentric anomaly input in radians
    # e = eccentricity (0 - 1)

    # returns dg/dt as provided in HW
    return 1 - e * cos(E)
end


# a newton solver for Keplers equation
function newtonSolver(e::Float64, M::Float64)::Float64

    # e = eccentricity (0-1)
    # M = mean anomaly, input in radians (0-2pi)

    # returns
    # E = Eccentric anomaly in radians, iterations took, and last error

    # starting point for solver, starting point suggested in lecture
    E_o = M + 0.85 * e * sign(sin(M))
    # exit tolerance
    epsilon = 1e-15
    maxIter = 1000000
    iter = 1

    while iter < maxIter
        global E_next # eventual return value
        global E_error
        iter += 1
        der = dg_dt(E_o, e)
        if der == 0
            der = 1e-16 # avoid zero division
        end
        E_next = E_o - g(E_o,e,M) / der
        E_error = abs(E_o-E_next)
        if E_error < epsilon
            break # stop iteration we've converged
        end
        # update value and continue looping
        E_o = E_next
    end

    return E_next
end


#### parse input file for mystery plant
# columns are time, rv (in m/s) and rv error

function sortByTime(ts, radVel, vErr)
    sortInds = sortperm(ts)
    ts = ts[sortInds]
    radVel = radVel[sortInds]
    vErr = vErr[sortInds]
    return ts, radVel, vErr
end

function foldByTime(ts, modTime)
    # return timeseries folded at modulo modTime
    return ts .% modTime
end

function plotTimeseries(ts, radVel, vErr, figname)
    # plot the time series
    plot(ts, radVel, markersize=1, markershape=:auto, line=:solid, linealpha=0.5, yerror=vErr, dpi=150)
    savefig(figname)
end

# sort by increasing time and normalize to t=0
# ts, radVel, vErr = sortByTime(ts, radVel, vErr)
# ts = ts .- ts[1]

# plotTimeseries(ts, radVel, vErr, "initTS")

function sumMinDist(radVel, vErr)
    v1 = radVel[1:end-1]
    v2 = radVel[2:end]
    e1 = vErr[1:end-1]
    e2 = vErr[2:end]
    sigmas = e1.^2 .+ e2.^2
    dv = v2 .- v1
    # compute array of successive distances
    # scaled by the errors
    # metric = sum(sqrt.(dt.^2 .+ dv.^2) ./ sigmas)

    metric = sum(dv.^2 ./ sigmas) # faster
    return metric

end

function periodSearch(ts, radVel, vErr)
    # a brute force approach
    # returns: bestPeriod, folded time vector, velocity vector, error vector
    timeStep = .001 # second
    minTime = 0
    maxTime = 1000 # by visual inspection
    foldTimes = minTime:timeStep:maxTime
    metrics = Array{Float64,1}(undef, length(foldTimes))
    minMetric = sumMinDist(radVel, vErr)
    bestP = -1 # period
    bestT = ts # times
    bestV = radVel # velocities
    bestE = vErr # errors
    for (i, tryP) in enumerate(foldTimes)
        tryT = foldByTime(ts, tryP)
        tryT, tryV, tryE = sortByTime(tryT, radVel, vErr)
        metric = sumMinDist(tryV, tryE)
        metrics[i] = metric
        if metric < minMetric
            minMetric = metric
            bestT = tryT
            bestV = tryV
            bestE = tryE
            bestP = tryP
        end
    end
    print("best period $bestP with metric $minMetric\n")
    plotTimeseries(bestT, bestV, bestE, "periodFolded")
    # plot metrics
    plot(foldTimes, metrics, dpi=150)
    savefig("metrics")
    return bestP, bestT, bestV, bestE
end

function vModel(t_i, P, e, t_p, gamma, omega, K)
    M = 2*pi/P*(t_i-t_p)
    E = newtonSolver(e, M)
    # print("E $E\n")
    f = 2 * atan( ((1+e)/(1-e))^0.5 * tan(E/2) )
    h = K * cos(omega)
    c = -K * sin(omega)
    v_o = gamma + K*e*cos(omega)
    v_rad = h*cos(f) + c*sin(f) + v_o
    return v_rad
end

mp = readdlm("mystery_planet1.txt")
ts = mp[:,1]
radVel = mp[:,2]
vErr = mp[:,3]

# sort by increasing time
ts, radVel, vErr = sortByTime(ts, radVel, vErr)
ts = ts .- ts[1] # normalize time to begin at zero
period, tsFold, radVelFold, vErrFold = periodSearch(ts, radVel, vErr)
nVals = length(ts)

Es = zeros(length(ts)) # reusable!
A = zeros(3,3) # reusable
b = zeros(3) # reusable
invSigma2 = vErr .^ -2 # just compute this once

function newton2(e,M)
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

function solveOrbit(e, period, t_p)
    # e = x[1]
    # period = x[2]
    # t_p = x[3]
    # period = 111.487

    Ms = 2*pi/period .* ts .- t_p
    Ms = Ms .% (2*pi)
    # use keplers eqn solver to get the Es
    for (ii, M) in enumerate(Ms)
        ### cut and past in slover here to increase speed?
        Es[ii] = newton2(e,M)
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

    # A[2,1] = sum(invSigma2 .* sincosfs)
    A[2,1] = A[1,2]
    A[2,2] = sum(invSigma2 .* sin2fs)
    A[2,3] = sum(invSigma2 .* sinfs)

    # A[3,1] = sum(invSigma2 .* cosfs)
    A[3,1] = A[1,3]
    # A[3,2] = sum(invSigma2 .* sinfs)
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

function minSolveOrbit(x)
    e = x[1]
    period = x[2]
    t_p = x[3]
    chi2, h, c, v_o = solveOrbit(e, period, t_p)
    return chi2
end

function minSolveFixP(x)
    e = x[1]
    t_p = x[2]
    chi2, h, c, v_o = solveOrbit(e, period, t_p)
    return chi2
end

# K = (maximum(radVel)-minimum(radVel)) / 2
# print("K $K\n")
x = zeros(3)
x[1] = 0.5 # e
x[2] = period
x[3] = period / 2 # t_p
print("x $x\n")

# bestChi2 = 1e16

using Optim
output = optimize(minSolveOrbit, x)
print("solver output $output\n")
eSolve, pSolve, t_pSolve = Optim.minimizer(output)
print("solved e $eSolve, solved P $pSolve\n")

# x2 = zeros(2)
# x[1] = 0.5
# x[2] = period / 3
# output = optimize(minSolveFixP, x2)
# print("solver output $output\n")
# eSolve, pSolve, t_pSolve = Optim.minimizer(output)
# print("solved e $eSolve, solved P $pSolve\n")

chi2, h, c, v_o = solveOrbit(eSolve, pSolve, t_pSolve)
print("chi2 $chi2\n")
print("h $h c $c v_o $v_o\n")

# h = K cos(omega)
# c = -K sin(omega)
# v_o = gamma + K * e cos(omega)

# convert h, c, v_o to gamma, K, omega
gamma = v_o - eSolve * h # m/s
K = sqrt(h^2+c^2) # m/s
omega = acos(h / K) * 360 / (2*pi) # degrees
print("period $pSolve e $eSolve, t_p $t_pSolve gamma $gamma K $K omega $omega\n")




# bestChi = 1e16
# bestOut = zeros(4)
# bestIn = zeros(2)
# for e in 0:.01:0.99
#     print("e: $e\n")
#     for t_p in 0:.001:period
#         global bestChi
#         out = solveOrbit(e, period, t_p)
#         if out[1] < bestChi
#             global bestIn
#             global bestOut
#             global bestChi
#             bestChi = chi2
#             bestOut = out
#             print("e $e\n")
#             bestIn[1] = e
#             bestIn[2] = t_p
#         end
#     end
# end


# print("bestIn $bestIn bestOut $bestOut\n")





# function minChi(x)
#     # x contains array of values to find
#     # x = e, t_p, gamma, K
#     e = x[1]
#     t_p = x[2]
#     gamma = x[3]
#     omega = x[4]
#     K = x[5]
#     chi2 = 0
#     for i in 1:nVals
#         v_i = radVel[i]
#         t_i = ts[i]
#         e_i = vErr[i]
#         chi2 += (v_i - vModel(t_i, P, e, t_p, gamma, omega, K))^2 / e_i^2
#     end
#     return chi2
# end

# x0 = [0, ts[5], 0, 0, initK]
# using Optim
# optimize(minChi, x0)
# print("x0 $x0\n")


# def
# M_planet =







# HD 80606 b



