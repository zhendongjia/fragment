#!/usr/bin/env python
import random
import math
import os
import shutil


fragment_location = os.getenv('HOME')+'/fragment/fragment'
runfile_location = os.getenv('HOME')+'/fragment/runfile'

run_start = 2
run_end = 2

zmoon = 7.34767309e22/1.9891e30
mj = 0.0009546
aj = 5.2
ej = 0.08

for loop in range (run_start - 1, run_end):
    if not os.path.exists('%d'%(loop+1)):
        os.mkdir('%d'%(loop+1))
    os.chdir('%s/%d'%(runfile_location,(loop+1)))
    print '%s/%d'%(runfile_location,(loop+1))

    kstart = 1
    tcomp = 100000.0

    n = 2
    nbtot = 100
    nbpert = 10
    nrpert = 2
    nrand = 100
    nran = 1

    eta = 0.001
    deltat = 100.0
    tcrit = 1e6
    dfmin = 0.001
    rin = 1.0
    rout = 1.2
    scale = 3
    zmax = 0.0

    kz = [None] + [1,1,0,0,2,1,0,0,2,1,0,0,0,0,1,0,0,0,0,0]

    add_j = 1
    gas_p = 1
    gas_d = 0

    if add_j == 1:
        kz[3] = 1
        kz[17] = 1
    if gas_p == 1:
        kz[18] = 1
    if gas_d == 1:
        kz[19] = 1

        
    ci = 0.7
    cii = 0.5
    vc = 5.5e1
    zkm = 1.0e-7
    si = 3.0e7
    cej = 3.0e6

    zm1 = 1.0
    zm2 = 170.0

    rscale = 1.0

    np = n 
    mp = [zmoon]*np
#    radius = [rin + (rout - rin)*random.uniform(0,1) for i in range(np)]
    radius = [1,3]
    ecc = [0]*np
    a = [radius[i]/(1-ecc[i]) for i in range(np)]
    theta = [random.uniform(0,2*math.pi) for i in range(np)]
    phi = [random.uniform(0,2*math.pi) for i in range(np)]
    x = [radius[i]*math.cos(theta[i]) for i in range(np)]
    y = [radius[i]*math.sin(theta[i]) for i in range(np)]
    z = [zmax*math.cos(i) for i in phi]
    vp = [math.sqrt((1+mp[i])/math.sqrt(x[i]**2+y[i]**2)) for i in range(np)] 
    vx = [-vp[i]*math.sin(theta[i]) for i in range(np)]
    vy = [vp[i]*math.cos(theta[i]) for i in range(np)]
    vz = [-zmax*math.sin(i) for i in phi]

    t_dep = 1e5
    r_edge = 12
    r_in = 4
    dens0 = 200
    dens0 = dens0*1.125e-7
    
    f = open('input','w')
    print >> f, kstart, tcomp
    print >> f, n, nbtot, nbpert, nrpert, nrand, nran
    print >> f, eta, deltat, tcrit, dfmin, rin, rout, scale, zmax
    for i in kz[1:]: print >> f, i,
    print >> f
    print >> f, ci, cii, vc, zkm, si, cej
    if kz[15] != 0: print >> f, rscale
    for i in range(np): print >> f, mp[i], x[i], y[i], z[i], vx[i], vy[i], vz[i]
    if kz[13] != 0: print >> f, zm1, zm2
    if kz[3] == 1: print >> f, mj, aj, ej
    if kz[18] == 1: print >> f, t_dep, r_edge, r_in, dens0
    f.close()

    os.system("%s<input>output" % fragment_location)

    os.chdir('%s'%runfile_location)
