import numpy
import multiprocessing
from scipy import optimize
import time
import matplotlib.pyplot as plt
import rebound


AU = 149597870700  # m
secondsPerDay = 86400
mSun = 1.989e30  # kg
rSun = 6.95700e8 # m
mJupiter = 1.898e27  # kg
mEarth = 5.972e24  # kg
mStar = 0.32 * mSun  # (paper)
nTimesteps = 10000

G = 6.67408e-11  # m^3 kg^-1 s^-2


def getMassPlanet(x, K, P, e, i):
    # use this in a minimizer to solve for mass
    # must be in meters seconds radians kgs units and angles in radians!
    mp = x[0]
    a = 2 * numpy.pi * G / (mStar + mp)**2
    b = (mp * numpy.sin(i))**3
    c = P * (1 - e**2)**(3 / 2.)
    return K - (a * b / c)**(1 / 3.)


class Planet(object):
    def __init__(self, name, P, a, K, e, omega, MA, mPublished):
        # inputs in seconds, kg, radians, meters
        # P period in seconds
        # a semi-major axis in m
        # K velocity semi amplitude m/s
        # e eccentricity
        # omega: radians
        # MA: radians
        # mPublished = published mass value
        self.name = name
        self.P = P
        self.a = a
        self.K = K
        self.e = e
        self.omega = omega
        self.MA = MA
        self.mPublished = mPublished

        # lists for each model run
        self.i = []
        self.m = []
        self.model_a = []
        self.model_e = []
        self.model_x = []
        self.model_y = []

    def addInclination(self, i):
        # update inclination and mass
        self.i.append(i)
        sol = optimize.root(getMassPlanet, [self.mPublished], args=(self.K, self.P, self.e, i))
        self.m.append(float(sol.x[0]))
        self.model_a.append(numpy.zeros(nTimesteps))
        self.model_e.append(numpy.zeros(nTimesteps))
        self.model_x.append(numpy.zeros(nTimesteps))
        self.model_y.append(numpy.zeros(nTimesteps))


planetD = Planet(
    name="d",
    P=1.937780 * secondsPerDay,
    a=0.02080665 * AU,
    K=6.56,
    e=0.207,
    omega=numpy.radians(234),
    MA=numpy.radians(355),
    mPublished=6.83 * mEarth,
)

planetC = Planet(
    name="c",
    P=30.0881 * secondsPerDay,
    a=0.129590 * AU,
    K=88.34,
    e=0.25591,
    omega=numpy.radians(48.76),
    MA=numpy.radians(294.59),
    mPublished=0.7142 * mJupiter,
)

planetB = Planet(
    name="b",
    P=61.1166 * secondsPerDay,
    a=0.208317 * AU,
    K=214.00,
    e=0.0324,
    omega=numpy.radians(50.3),
    MA=numpy.radians(325.7),
    mPublished=2.2756 * mJupiter,
)

planetE = Planet(
    name="e",
    P=124.26 * secondsPerDay,
    a=0.3343 * AU,
    K=3.42,
    e=0.055,
    omega=numpy.radians(239),
    MA=numpy.radians(335),
    mPublished=14.6 * mEarth,
)

planets = [planetD, planetC, planetB, planetE]

inclinationList = numpy.arange(1, 21)
endTime = planetD.P * 1000  # 1000 orbits of shortest period planet
times = numpy.linspace(0, endTime, nTimesteps)
for inc in inclinationList:
    # print("simulating inclination %.2f" % inc)
    for planet in planets:
        planet.addInclination(numpy.radians(inc))
    sim = rebound.Simulation()
    sim.units = ('s', 'm', 'kg')
    sim.integrator = "whfast"
    sim.dt = planetD.P / 100.0  # shortest period

    sim.add(m=mStar)
    for planet in planets:
        sim.add(
            m=planet.m[-1],
            P=planet.P,
            M=planet.MA,
            omega=planet.omega,
            e=planet.e
        )
    sim.move_to_com()

    for i, t in enumerate(times):
        sim.integrate(t, exact_finish_time=0)

        for j, planet in enumerate(planets):
            planet.model_x[-1][i] = sim.particles[j+1].x
            planet.model_y[-1][i] = sim.particles[j+1].y
            planet.model_a[-1][i] = sim.particles[j+1].a
            planet.model_e[-1][i] = sim.particles[j+1].e

# plot results
fig, axGrid = plt.subplots(int(len(inclinationList)//4), 4, figsize=(14, 12))
axList = axGrid.flatten()
for i, (inc, ax) in enumerate(zip(inclinationList, axList)):
    for planet in planets:
        ax.plot(times, planet.model_a[i] / AU, label="%s"%planet.name)
        if inc < 17:
            ax.get_xaxis().set_ticks([])
        else:
            ax.set_xlabel("time (seconds)\n1000 orbits of planet d")
        if i % 4 == 0:
            ax.set_ylabel("semi-major axis (AU)")
        ax.set_ylim([0, planetE.a*1.5 / AU])
        ax.set_title("i = %i degrees"%inc)

axList[0].legend()
# figure out reasonable labels
plt.savefig("inc_vs_a.png", dpi=250)
plt.close(fig)

### upper limit masses
ind = numpy.argwhere(inclinationList==13)[0][0] # index where inclination == 13
print("\nupper limit masses at inclination 13 deg:")
for planet in planets:
    print("planet %s mass: %.2f x (Earth Mass)"%(planet.name, planet.m[ind] / mEarth))


## impact parameter from  seager:
# b = a cos(i) / R_star * (1-e**2)/(1+esin(omega))
rStar = 0.3 * rSun

print("\nTransit inclinations:")
for planet in planets:
    i = numpy.degrees(numpy.arccos(rStar / planet.a * (1 + planet.e * numpy.sin(planet.omega)) / (1 - planet.e**2)))
    print("planet %s transits at inclination: %.2f" % (planet.name, i))


i = numpy.radians(86.66)
print("\nlower limit masses at inclination 86.66 deg:")
for planet in planets:
    sol = optimize.root(getMassPlanet, [planet.mPublished], args=(planet.K, planet.P, planet.e, i))
    m = float(sol.x[0])
    print("planet %s mass: %.2f x (Earth Mass)"%(planet.name, m / mEarth))


################ compute radial velocity ##########################################
# we look from +z at system
for inc in [0, 59]:
    # publishedInclination = numpy.radians(59)
    sim = rebound.Simulation()
    sim.units = ('s', 'm', 'kg')
    sim.integrator = "whfast"
    sim.dt = planetD.P / 100.0  # shortest period

    sim.add(m=mStar)
    for planet in [planetC, planetB]:
        sim.add(
            m=planet.mPublished,
            P=planet.P,
            M=planet.MA,
            omega=planet.omega,
            e=planet.e,
            inc=numpy.radians(inc),
        )
    sim.move_to_com()

    fig = rebound.OrbitPlot(sim, trails=True, periastron=True)
    plt.title("inclination = %i degrees"%inc)
    plt.savefig("orbit_%i.png"%inc)
    plt.close(fig)

    radialVx = numpy.zeros(nTimesteps)
    radialVy = numpy.zeros(nTimesteps)
    radialVz = numpy.zeros(nTimesteps)

    for i, t in enumerate(times):
        sim.integrate(t, exact_finish_time=0)
        radialVx[i] = sim.particles[0].vx
        radialVy[i] = sim.particles[0].vy
        radialVz[i] = sim.particles[0].vz
    plt.figure()
    # scale the radial velocity by the inclination angle
    plt.plot(times, radialVx, label="vx")
    plt.plot(times, radialVy, label="vy")
    plt.plot(times, radialVz, label="vz - radial")
    plt.title("inclination = %i degrees"%inc)
    plt.xlabel("time (seconds)")
    plt.ylabel("velocity (m/s)")
    plt.ylim([-400,400])
    plt.legend()
    plt.savefig("rv_%i.png"%inc)
    plt.close()


########## fit 4 keplarians to vx ################

def kepEq(x, M, e):
    E = x[0]
    return M - (E - e * numpy.sin(E))


def vRad(K, omega, e, Mlist):
    # ignore barycentric velocity
    h = K * numpy.cos(omega)
    c = -K * numpy.sin(omega)
    v_o = K * e * numpy.cos(omega)
    Es = numpy.zeros(len(Mlist))
    for ii, M in enumerate(Mlist):
        EGuess = [M + 0.85 * e * numpy.sin(numpy.sin(M))]
        sol = optimize.root(kepEq, EGuess, args=(M, e))
        E = float(sol.x[0])
        Es[ii] = E

    f = 2 * numpy.arctan(((1 + e) / (1 - e)**0.5 * numpy.tan(Es / 2)))
    return h * numpy.cos(f) + c * numpy.sin(f) + v_o

# def solveKeplarians(x, modelV):
#     # run this through a minimizer to try and tweak the orbital paramters such
#     # that the best match to the rebound simulated rv is obtained by a
#     # series of 4 keplarians.

#     # reshape x in to rows and columns
#     # columns are P, MA, omega, e, K (at time 0)
#     # rows are planet d c b e
#     x = numpy.asarray(x).reshape(2,5)
#     vRads = []
#     for P, MA, omega, e, K in x:
#         P = abs(P)
#         e = abs(e)
#         K = abs(K)
#         MList = MA + 2 * numpy.pi / P * (times - times[0])
#         v = -1*vRad(K, omega, e, MList)
#         vRads.append(v)
#     vRads += vRadConst # const from planet parameters we aren't varying
#     vSum = numpy.sum(vRads, axis=0)
#     chi2 = numpy.sum((vRads-vSum)**2)
#     return chi2


def getF(e, P, M):
    ### solve keplers eqn with newtons method
    Ms = M + 2 * numpy.pi / P * (times - times[0])
    Es = Ms + 0.85 * e * numpy.sign(numpy.sin(Ms))
    epsilon = 1e-8
    _iter = 0
    while True:
        _iter += 1
        if _iter > 1000:
            print("max iter reached!!!")
            raise RuntimeError("gah!")
        dEs = 1 - e * numpy.cos(Es)
        dEs[dEs==0] = 1e-16 # protect against zero division
        E_next = Es - (Es - e * numpy.sin(Es) - Ms) / dEs
        E_error = numpy.max(numpy.abs(Es - E_next))
        if E_error < epsilon:
            break
        Es = E_next
    fs = 2 * numpy.arctan( ((1+e)/(1-e))**0.5 * numpy.tan(Es/2.0) )
    return fs

def fitKeplers(x, nPlanets, measV):
    x = numpy.asarray(x).reshape(nPlanets,3)
    fs = [] # holds a list of fs for each planet
    for e, P, M in x:
        fs.append(getF(e, P, M))
    fs = numpy.asarray(fs)
    # linear solver for h, c, and v_o
    cosfs = numpy.cos(fs)
    sinfs = numpy.sin(fs)

    stack = numpy.vstack((cosfs, sinfs, numpy.ones(len(fs[0]))))

    # linear solver for h, c, and v_o
    # build array A
    A = []
    for element in stack:
        A.append(numpy.sum(element*stack, axis=1))
    A = numpy.asarray(A)

    # multiply every row of stack by measV
    measV = numpy.array([measV]*len(stack))
    b = numpy.sum(stack*measV, axis=1)

    coeffs = numpy.linalg.solve(A, b)
    rv = numpy.zeros(len(fs[0]))
    for coeff, row in zip(coeffs, stack):
        rv += coeff*row
    # next compute the chi2 and return it
    chi2 = numpy.sum((rv - measV)**2)
    return chi2, rv

### solve radial velocities for lest massive planets
# we wont vary these
# vRadConst = []
# for planet in [planetD, planetE]:
#     vRadConst.append(getRadVel(planet.e, planet.P, planet.MA))


# def solveKeplarians2(x, vRadMeas):
#     vRads = []
#     x = numpy.asarray(x).reshape(2,3)
#     for e, P, M in x:
#         vRads.append(getRadVel(e, P, M))
#     vRads += vRadConst # const from planet parameters we aren't varying




vRads = []
plt.figure(figsize=(10,10))
for planet in [planetC, planetB]:
    MList = planet.MA + 2 * numpy.pi / planet.P * (times - times[0])
    v = vRad(planet.K, planet.omega, planet.e, MList)
    v = -1*v # this was necessary
    plt.plot(times, v, label="planet %s"%planet.name)
    vRads.append(v)
plt.legend()
plt.xlabel("time (seconds)")
plt.ylabel("radial velocity (m/s)")
plt.savefig("vr_components.png", dpi=150)
plt.close()

vRads = numpy.asarray(vRads)
vRadSum = numpy.sum(vRads, axis=0)
plt.figure(figsize=(10,10))
plt.plot(times, radialVz, label="rebound (WH-Fast) vz")
plt.plot(times, vRadSum, label="sum of keplarians")
plt.legend()
plt.xlabel("time (seconds)")
plt.ylabel("radial velocity (m/s)")
plt.savefig("reboundVsKepler.png")
plt.close()

### look at residuals
plt.figure(figsize=(10, 10))
plt.plot(times / planetD.P, radialVz - vRadSum, label="WH-Fast - sum(Keplers)")
p = plt.axhspan(-10, 10, facecolor='red', alpha=0.2, label="+/- 10 m/s")
p = plt.axhspan(-5, 5, facecolor='green', alpha=0.2, label="+/- 5 m/s")
plt.legend()
plt.ylabel("radial velocity residulas (m/s)")
plt.xlabel("innermost planet orbits")
plt.savefig("residuals.png", dpi=150)
plt.close()


#### try to tweak orbital params for a better fit
initialParams = []
for planet in [planetC, planetB]:
    initialParams += [planet.e, planet.P, planet.MA]
initialParams = numpy.array(initialParams)


nPlanets = 2
t1 = time.time()
chi2, rv = fitKeplers(initialParams, nPlanets, radialVz)

print("init chi2 %.5f"%chi2)

print("kep fit took %.2f"%(time.time()-t1))
plt.figure(figsize=(10,10))
plt.plot(times, radialVz, label="rebound (WH-Fast) vz")
plt.plot(times, rv, label="simultaneously fit keplerains")
plt.legend()
plt.xlabel("time (seconds)")
plt.ylabel("radial velocity (m/s)")
plt.savefig("simul.png")
plt.close()

### look at residuals
plt.figure(figsize=(10, 10))
plt.plot(times / planetD.P, radialVz - rv, label="WH-Fast - simul")
p = plt.axhspan(-10, 10, facecolor='red', alpha=0.2, label="+/- 10 m/s")
p = plt.axhspan(-5, 5, facecolor='green', alpha=0.2, label="+/- 5 m/s")
plt.legend()
plt.ylabel("radial velocity residulas (m/s)")
plt.xlabel("innermost planet orbits")
plt.savefig("simulresiduals.png", dpi=150)
plt.close()

nParams = len(initialParams)

varyMs = numpy.radians(numpy.linspace(0,359.9,45))
varyEs = numpy.linspace(1/2.0, 2.0, 10)
gridParams = []
for p1M in varyMs:
    for p1e in varyEs:
        p1 = [planetC.e*p1e, planetC.P, p1M]
        for p2M in varyMs:
            for p2e in varyEs:
                gridParams.append(p1+[planetB.e*p2e, planetB.P, p2M])
gridParams = numpy.array(gridParams)



def parallelGridSearch(gridInds):
    bestChi2 = 1e16
    bestRV = None
    bestGridInd = None
    for ii, gridInd in enumerate(gridInds):
        # scalings = numpy.array(scalingGrid[gridInd])
        params = gridParams[gridInd]
        if ii % 1000 == 0:
            print("%.0f percent done"%(100*ii/float(len(gridInds))))
        chi2, rv = fitKeplers(params, nPlanets, radialVz)
        if chi2 < bestChi2:
            bestChi2 = chi2
            bestRV = rv
            bestGridInd = ii
    return bestChi2, bestGridInd, bestRV


allGridInds = numpy.arange(len(gridParams))
splitGridInds = numpy.array_split(allGridInds, 12)
p = multiprocessing.Pool(12)
t1 = time.time()
res = p.map(parallelGridSearch, splitGridInds)
minInd = numpy.argmin([x[0] for x in res])
bestRV = res[minInd][-1]
bestVals = gridParams[res[minInd][1]]
print("grid search found chi2 (%.2f): "%res[minInd][0])
print("planet c | e=%.2f P=%.4f (days) MA=%.4f (degrees)"%(bestVals[0], bestVals[1]/secondsPerDay, numpy.degrees(bestVals[2])))
print("planet b | e=%.2f P=%.4f (days) MA=%.4f (degrees)"%(bestVals[0+3], bestVals[1+3]/secondsPerDay, numpy.degrees(bestVals[2+3])))

print("grid search took : %.2f"%(time.time()-t1))


### look at residuals
plt.figure(figsize=(10, 10))
plt.plot(times / planetD.P, radialVz - bestRV, label="WH-Fast - simul")
p = plt.axhspan(-10, 10, facecolor='red', alpha=0.2, label="+/- 10 m/s")
p = plt.axhspan(-5, 5, facecolor='green', alpha=0.2, label="+/- 5 m/s")
plt.legend()
plt.ylabel("radial velocity residulas (m/s)")
plt.xlabel("innermost planet orbits")
plt.savefig("bestsimulresiduals.png", dpi=150)
plt.close()

def kepMinimizer(params, nPlanets, measV):
    chi2, rv = fitKeplers(params, nPlanets, measV)
    return chi2


### use a minimizer!
initialGuess = gridParams[res[minInd][1]]
# initialGuess = [2.55910000e-01, 2.59961184e+06, 2.12989332e-01, 6.48000000e-02, 5.28047424e+06, 1.38443066e+00]
paramLowerBounds = [0, planetC.P/200., 0.1, 0, planetB.P/200., 0.1]
paramUpperBounds = [1, planetC.P*200, numpy.radians(359.9), 1, planetB.P*200, numpy.radians(359.9)]

bounds = optimize.Bounds(paramLowerBounds, paramUpperBounds)
res = optimize.minimize(kepMinimizer, initialParams, args=(nPlanets, radialVz), method='trust-constr', options={'verbose': 1}, bounds=bounds)
bestVals = res.x

chi, bestRV = fitKeplers(bestVals, nPlanets, radialVz)
print("minimizer chi2: %.5f"%chi)
print("minimizer found: ")
print("planet c | e=%.2f P=%.4f (days) MA=%.4f (degrees)"%(bestVals[0], bestVals[1]/secondsPerDay, numpy.degrees(bestVals[2])))
print("planet b | e=%.2f P=%.4f (days) MA=%.4f (degrees)"%(bestVals[0+3], bestVals[1+3]/secondsPerDay, numpy.degrees(bestVals[2+3])))

plt.figure(figsize=(10, 10))
plt.plot(times / planetD.P, radialVz - bestRV, label="WH-Fast - minimize")
p = plt.axhspan(-10, 10, facecolor='red', alpha=0.2, label="+/- 10 m/s")
p = plt.axhspan(-5, 5, facecolor='green', alpha=0.2, label="+/- 5 m/s")
plt.legend()
plt.ylabel("radial velocity residulas (m/s)")
plt.xlabel("innermost planet orbits")
plt.savefig("bestminimizedresiduals.png", dpi=150)
plt.close()

################## TTV ########################
secondsPerYear = 3.154e7
endTime = 4 * secondsPerYear  # 1000 orbits of shortest period planet
plt.figure(figsize=(10,10))
for sm in [mStar/2, mStar, 2*mStar]:
    sim = rebound.Simulation()
    sim.units = ('s', 'm', 'kg')
    sim.integrator = "whfast"
    sim.dt = planetD.P / 100.0  # shortest period
    sim.add(m=sm)
    for planet in planets:
        sim.add(
            m=planet.mPublished,
            P=planet.P,
            M=planet.MA,
            omega=planet.omega,
            e=planet.e,
        )
    sim.move_to_com()

    transitTimes = []
    tStep = planetC.P / 100.0 # choose something less than C's period
    p = sim.particles
    # try catch crossings when the y coordinate moves from negative to positive
    while sim.t < endTime:
        # compare y coords of star and planet c
        # assumes the observer is at +x
        y_old = p[2].y - p[0].y # planet c - star
        t_old = sim.t
        sim.integrate(sim.t + tStep)
        t_new = sim.t
        if y_old*(p[2].y-p[0].y)<0 and p[2].x-p[0].x>0:
            # sign changed (y_old*y<0), planet in front of star (x>0)
            while t_new - t_old > 1e-7:
                if y_old*(p[2].y-p[0].y)<0:
                    t_new = sim.t
                else:
                    t_old = sim.t
                sim.integrate((t_new+t_old)/2.)
            transitTimes.append(sim.t)
            sim.integrate(sim.t+tStep)

    transitTimes = numpy.asarray(transitTimes)

    nTransits = len(transitTimes)
    transitNumber = range(nTransits)
    p = numpy.polyfit(transitNumber, transitTimes, deg=1)
    fitLine = p[0] * transitNumber + p[1]
    plt.plot(transitTimes / secondsPerYear, (transitTimes - fitLine)/3600.0, '-o', alpha=0.7, label="star mass = %.2f x (star mass)"%(sm / mStar))
plt.xlabel("time (years)")
plt.ylabel("TTV (hours)")
plt.legend()
plt.savefig("TTV.png")
plt.close()


# # ffmpeg -r 10 -f image2 -i step_%06d.png -pix_fmt yuv420p orbit.mp4



