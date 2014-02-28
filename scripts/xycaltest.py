from threading import Thread
from approachimagingparamspark import *

xycalstart = []
ls = linescanthread(afm, [0])
ls.start()
while ls.is_alive():
    xy = [getattr(daq, 'get_ai%i' % XYCHANS[0]), 
          getattr(daq, 'get_ai%i' % XYCHANS[1])]
    xycalstart.append(xy)
    
xycalend = []
ls = linescanthread(afm, [ypoints])
ls.start()
while ls.is_alive():
    xy = [getattr(daq, 'get_ai%i' % XYCHANS[0]), 
          getattr(daq, 'get_ai%i' % XYCHANS[1])]
    xycalend.append(xy)

#Calibrate XY reading:
xmin = min(xycalend)
xmax = max(xycalstart)
ymin = min(xycalstart)
ymax = max(xycalend)

class linescanthread(Thread):
    '''
    Takes linescan as thread.
    Use with linescanthread.start()
    
    Args:
    afm: pass in afm instrument
    lines: list of lines to scan
    '''
    def __init__(self, afm, lines):
        Thread.__init__self()
        self.afm = afm
        self.lines = lines
    def run(self):
        for line in self.lines:
            self.currentline = line
            self.afm.LineScan(line)