#!/usr/bin/env python

import Gnuplot
Gnuplot.GnuplotOpts.default_term ='post enh color eps'
import math
import os
import resource

begin = 10
end = 10
runfile_location=os.getenv('HOME')+'/fragment/runfile'

class Orbit:
    def __init__(self, line):
        seg = line.split()
        start = seg.index(':') + 1
        self.name = int(seg[start])
        (self.m, self.semi, self.ecc, self.x, self.y, self.z, 
         self.vx, self.vy, self.vz, self.t1, self.t2, self.w) = [
            float(seg[start+i+1]) for i in range(12)]

    def __str__(self):
        return 'name=%d, m=%f, semi=%f, ecc=%f, w=%f' % (
            self.name, self.m, self.semi, self.ecc, self.w) 


def getstyle(name):
    h = '%06X' % hash(str(name))
    return 'l ls 1 lc rgb "#%s"' % h[-6:]


class PlotUtil:
    def __init__(self, name, xmax, ymax):
        self.prev = { }
        self.g = Gnuplot.Gnuplot()
        self.g('set term post enh color eps')
        self.g('set output "%s.eps"' % name)
        self.g('unset key')
        self.g('set log')
        self.g.set_range('xrange', (100, xmax))
        self.g.set_range('yrange', (0.9, ymax))



    def finish(self):
        for name, values in self.prev.items():
            p = Gnuplot.PlotItems.Data(values, with_=getstyle(name),
                                       title=str(name))
            self.g._add_to_queue([p])
        self.g.refresh()
        self.g.close()

    def add(self, name, time, value, prev=[]):
        for i in prev:
            if i.name==name: continue
            p = Gnuplot.PlotItems.Data(
                (self.prev[i.name][-1], (time,value)),
                with_=getstyle(i.name))
            self.g._add_to_queue([p])
        if name not in self.prev: self.prev[name]=[]
        self.prev[name].append((time,value))

# Class to store the collision event.
class Collision:
    def __init__(self, output, i):
        seg = output[i].split()
        self.time = float(seg[3])
        self.src = [Orbit(output[i+1]), Orbit(output[i+2])]
        i = i + 3
        while output[i].startswith('   DIAG'): i=i+1
        if output[i].startswith('     MERGE:'):
            while not output[i].startswith('     MERGER :'): i+=1
            self.dst = [Orbit(output[i])]
        elif output[i].startswith('     FRAGMT:'):
            self.dst = [ ]
            while not output[i].startswith('     FRAGM. :'): i+=1
            while output[i].startswith('     FRAGM. :'):
                self.dst.append(Orbit(output[i]))
                i += 1
            while not output[i].startswith('     MERGER :'): i+=1
            self.dst.append(Orbit(output[i]))
        else: raise(InputError, output[i], 'Invalid line %d' % i)


#resource.setrlimit(resource.RLIMIT_NOFILE, (2048, 8192))
for loop in range(begin-1, end):
    os.chdir('%s/%d'%(runfile_location,(loop+1)))
    collisions = [ ]
    output = open('output').readlines()
    for i in range(len(output)):
        if output[i].startswith('   COLLISION'):
            collisions.append(Collision(output, i))

    fc = open('collisions', 'w')
    for c in collisions:
        print >>fc, 'Collision @ %f' % c.time
        for i in c.src: print >>fc, 'src: ', i
        for i in c.dst: print >>fc, 'dst: ', i
        print >>fc

    tcrit = float(output[12].split()[2])
    plot_semi = PlotUtil('semi', tcrit, 5.2)
    plot_ecc = PlotUtil('ecc', tcrit, 1)

    cid = 0
    time = 0.0

    for line in open('fort.3').readlines():
        if line.startswith('N T SCALE DUMMY'):
            time = float(line.split()[5])
            while cid < len(collisions) and collisions[cid].time < time:
                ctime = collisions[cid].time
                for i in collisions[cid].src:
                    plot_semi.add(i.name, ctime, i.semi)
                    plot_ecc.add(i.name, ctime, i.ecc)
                for i in collisions[cid].dst:
                    plot_semi.add(i.name, ctime, i.semi, collisions[cid].src)
                    plot_ecc.add(i.name, ctime, i.ecc, collisions[cid].src) 
                cid += 1
        elif line.startswith('M R V'):
            orbit = Orbit(line)
            plot_semi.add(orbit.name, time, orbit.semi)
            plot_ecc.add(orbit.name, time, orbit.ecc)
    plot_semi.finish()
    plot_ecc.finish()

    os.chdir(runfile_location)

