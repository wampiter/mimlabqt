import numpy as np
import msvcrt
import qt
import logging

import approachimagingparamspark
reload (approachimagingparamspark)
from approachimagingparamspark import * #used like C constants, always in CAPS
import taskclasses as tc
reload(tc)
import measurement_general as mg
from threading import Thread

def measure(feedback=False, xsize = 0, ysize = 0, angle = 0,
			xpoints = 128, ypoints = 128):
    '''
    Takes a topo measurement. Aborts with 'q' (or stop)
    in non-feedback, 'u' and 'd' will move stage up and down (backwards)
    
	Arg:
	xer: object for XER control window

    Kwargs:
    feedback -- whether or not to use feedback
    xsize -- scan size in um
    angle -- measurement ange (in radians)
    '''
    #Check that ysize not specified without xsize:
    if xsize == 0 and not ysize == 0:
        logging.error('Cannot have slow size but no fast size')
        raise ValueError(xsize, ysize)
        
    #set up qt meaurement flow
    qt.flow.connect('measurement-end', tc.kill_all_tasks)
    qt.mstart()
    #get instruments (or initialize if not already)
    if qt.instruments.get_instruments_by_type('NI_DAQ'):
        daq = qt.instruments.get_instruments_by_type('NI_DAQ')[0]
    else:
        daq = qt.instruments.create('daq', 'NI_DAQ', id = 'Dev1')
    if qt.instruments.get_instruments_by_type('xer'):
        afm = qt.instruments.get_instruments_by_type('xer')[0]
    else:
        afm = qt.instruments.create('afm', 'xer')

    #Store current height of scanner
    zstart = np.zeros(1)    
    zstart[0] = getattr(daq,'get_ao%i' % DCZCHAN[0])()
    if np.isnan(zstart[0]):
        daq.set('ao%i' % DCZCHAN[0], 0.0)
        zstart[0] = 0.0 

    #Prepare approach data structure to store absolutely all raw data:
    approach_data = qt.Data(name='approach_curves')
    approach_data.add_coordinate('sample')
    approach_data.add_value('optical')
    approach_data.add_value('z [mV]')
    approach_data.add_value('x [um]')
    approach_data.add_value('y [um]')
    approach_plot = qt.Plot2D(approach_data, name = 'approach_trace')
    approach_data.create_file()
    approach_data.copy_file('parkapproachtoposcan.py')
    approach_data.copy_file('approachimagingparamspark.py')
    #Prepare spatial data structures for "normal" scan data:
    if xsize > 0:
        spatial_data_right = qt.Data(name = 'spatial_data_right')
        spatial_data_left = qt.Data(name = 'spatial_data_left')
        spatial_data = [spatial_data_right, spatial_data_left]
        for data_obj in spatial_data:
            data_obj.add_value('x [mV]')
            data_obj.add_value('y [mV]')
            data_obj.add_value('z [mV]')
            data_obj.topoplot2d = qt.Plot2D(data_obj, coordim = 0, 
                                name = 'topography linecuts %s' % data_obj)
            if ysize > 0:
                data_obj.topoplot3d = qt.Plot3D(spatial_data, 
                                                name = 'topography')
            data_obj.create_file()
    else:
        spatial_data = False
    
    #Set up AFM:
    afm.set_fSizeX(xsize)
    afm.set_fSizeY(ysize)
    if xsize > 0:
        afm.set_fRotation(angle)
        afm.set_fRate(float(SAMPLERATE)/float(SAMPLES*xpoints))
        
        #take calibration scan:
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
        xminv = min(xycalend)
        xmaxv = max(xycalstart)
        yminv = min(xycalstart)
        ymaxv = max(xycalend)
    
    def voltagetoposition(vx,vy):
        x = ((xsize*np.cos(angle)+ysize*np.sin(angle))*(v-xminv)/(xmaxv-xminv)
            - ysize*np.sin(angle))
        y = (xsize*np.sin(angle)+ysize*np.cos(angle))*(v-yminv)/(ymaxv-yminv)
        return [x,y]

    zscan = afm.get_zScanner()
    if not zscan:
        zscan = 0
    
    #Set up tasks:
    zactask = tc.AcOutTask(DEV, ACZCHAN, SAMPLES, SAMPLE_RATE, sync = True)    
    zacdata = GenSineWave(SAMPLES, AMPLITUDE, PHASE)
    zactask.set_signal(zacdata)
    if feedback:
        ztask = tc.DcOutTask(DEV, [DCCHANS[2]])
    else:
        ztask = False
    maintask = mimCallbackTask(DEV, TOPOCHAN + XYCHANS, SAMPLES, SAMPLE_RATE,
                    afm, feedback, ztask, approach_data, spatial_data, zstart,
                    xsize, angle, xpoints, ypoints)
    
    maintask.userin = False
    zactask.StartTask()
    maintask.StartTask()
    
    while maintask.userin != 'q':
        maintask.userin = False
        if msvcrt.kbhit():
            maintask.userin = msvcrt.getch()
        try:
            qt.msleep(.2)
        except:
            logging.warning('Live plotting glitch (No big deal)')
        if not feedback:
	    if maintask.userin == 'w':
		maintask.z[0] += Z_STEP
		getattr(daq, 'set_ao%i' % DCZCHAN[0])(maintask.z[0])
		print maintask.z[0]
	    elif maintask.userin == 's':
                maintask.z[0] -= Z_STEP
        	getattr(daq, 'set_ao%i' % DCZCHAN[0])(maintask.z[0])
        	print maintask.z[0]
            elif maintask.userin == 'a':
                zscan -= Z_SCANNER_STEP
                afm.set_zScanner(zscan)
                print zscan
            elif maintask.userin == 'd':
                zscan += Z_SCANNER_STEP
                afm.set_zScanner(zscan)
                print zscan


    print('%i approaches completed' % maintask.callcounter)
    #stop tasks before error
    qt.mend()
    qt.msleep(.1)
    getattr(daq,'set_ao%i' % ACZCHAN[0])(0)
    #close all 3 data objects    
    approach_data.close_file()
    if spatial_data:
        for data_obj in spatial_data:
            data_obj.close_file()

class mimCallbackTask(tc.AnalogInCallbackTask):
    '''
    Creates (but does not start) task which runs continuously, executing 
    callback code every time an approach curve is taken

    Args:
    channels -- list of input channels to use, first C then R
    samples --
    samperate --
    zstart -- z-voltage at which to begin
    feedback -- whether or not to use topography (z) feedback
    '''
    def __init__(self, dev, channels, samples, samplerate, afm, feedback, 
                 ztask, approach_data, spatial_data, zstart, xsize, angle,
                 xpoints, ypoints):
        tc.AnalogInCallbackTask.__init__(self, dev, channels, samples, samplerate)
        self.z = zstart
        self.feedback = feedback
        self.userin = False # stores whether task has been interrupted by user
        self.callcounter = 0 # iterates each approach curve (callback)
        self.ztask = ztask
        self.approach_data = approach_data
        self.spatial_data = spatial_data
        self.xsize = xsize
        if xsize > 0:
            self.ls = linescanthread(afm, np.arange(ypoints))
            self.ls.start()
        self.currentline = 0
        self.state = 'start'
        self.xpoints = xpoints
        if angle > 90:
            self.rotated = True
        else:
            self.rotated = False

    def EveryNCallback(self):
        tc.AnalogInCallbackTask.EveryNCallback(self)    
        if self.xsize > 0:   
            if not self.ls.is_alive():
                self.userin = 'q'
                print 'Completed scan'
                self.StopTask()
                return
            else:
                if self.currentline != self.ls.currentline:
                    self.currentline = self.ls.currentline
                    print('on line %i' % self.currentline) 
                    for data_obj in self.spatial_data:
                        data_obj.new_block()
                    self.state = 'start'
                    prev_xy = False
                    
                        
        #separate out channels in most recent trace
        odata = self.split_data[0]
        if self.feedback:
           #For feedback, compare sample height around point of mass force to equilibrium value:
            contactheight = np.mean(odata[CONTACT_START:CONTACT_STOP])
            eqheight = np.mean(odata[EQ_START:EQ_STOP])
            displace = contactheight - eqheight
            delta = LOOP_GAIN * (displace - TARGET_DISPLACE)# Check sign
            if abs(delta) <.01:
                self.z[0] += delta
            else:
                self.z[0] += .01 * np.sign(delta)
            self.ztask.set_voltage(self.z)
        #Check that Z is within limits
        if self.z[0] > Z_MAX:
            self.z[0] = Z_MAX
            logging.warning('Reached maximum allowable value: %f' % Z_MAX)
        elif self.z[0] < Z_MIN:
            self.z[0] = Z_MIN
            logging.warning('Reached minimum allowable value: %f' % Z_MIN)
        #record raw data:
        if self.xsize > 0:
            [xpos, ypos] = voltagetoposition(np.mean(self.split_data[1]),
                                            np.mean(self.split_data[2]))
        else:
            xpos = ypos = 0
        self.approach_data.add_data_point(
                np.arange(SAMPLES), odata, self.z * 1e3 * np.ones(SAMPLES),
                            xpos*np.ones(SAMPLES), ypos*np.ones(SAMPLES))
        self.approach_data.new_block()

        #state checks
        if self.xsize > 0:
            if not self.rotated and prev_xy:
                frac_change = (xpos - prev_xy[0])/(self.xsize)
            elif prev_xy:
                frac_change = (ypos - prev_xy[1])/(self.xsize)
            if self.state == 'start' and self.prev_xy:
                if frac_change > 0.5/xpoints:
                    self.state = 'out'
            elif self.state == 'out':
                if frac_change < 0:
                    self.state = 'back'
            elif self.state == 'back':
                if frac_change > -0.5/xpoints:
                    self.state = 'done'
                    
            if self.state == 'out':
                spatial_data_current = self.spatial_data[0]
            elif self.state == 'back':
                spatial_data_current = self.spatial_data[1]
            else:
                spatial_data_current = False
            if spatial_data_current:
                spatial_data_current.add_data_point(
                    xpos, ypos, self.z[0] * 1e3)
            self.prev_xy = [xpos,ypos]
        self.callcounter += 1
        
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
            self.afm.LineScan(line)
