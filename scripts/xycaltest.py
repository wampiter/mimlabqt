from threading import Thread
from approachimagingparamspark import *

ypoints = 128

class linescanthread(Thread):
    '''
    Takes linescan as thread.
    Use with linescanthread.start()
    
    Args:
    afm: pass in afm instrument
    lines: list of lines to scan
    '''
    def __init__(self, afm, lines):
        Thread.__init__(self)
        self.afm = afm
        self.lines = lines
    def run(self):
        for line in self.lines:
            self.currentline = line
            self.afm.LineScan(line, True)

afm.set_fRate(.5)

xycalstart = []
ls = linescanthread(afm, [0])
ls.start()
while ls.is_alive():
    xy = [getattr(daq, 'get_ai%i' % XYCHANS[0])(), 
          getattr(daq, 'get_ai%i' % XYCHANS[1])()]
    xycalstart.append(xy)
    
xycalend = []
ls = linescanthread(afm, [ypoints])
ls.start()
while ls.is_alive():
    xy = [getattr(daq, 'get_ai%i' % XYCHANS[0])(), 
          getattr(daq, 'get_ai%i' % XYCHANS[1])()]
    xycalend.append(xy)

#Calibrate XY reading:
xmin = min([i[0] for i in xycalend])
xmax = max([i[0] for i in xycalstart])
ymin = min([i[1] for i in xycalstart])
ymax = max([i[1] for i in xycalend])

p1 = plot([i[0] for i in xycalstart])
print xmin
print xmax
print ymin
print ymax
