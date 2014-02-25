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
import parkstructs
import ctypes as c

def measure(xer, feedback=False, xsize = 0, ysize = 0, angle = 0,
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

    #Store current position of scanner
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
            data_obj.add_coordinate('x [mV]')
            data_obj.add_coordinate('y [mV]')
            data_obj.add_value('z [mV]')
            data_obj.topoplot2d = qt.Plot2D(
                    data_obj, coordim = 0, name = 'topography linecuts %s' % data_obj)
            if ysize > 0:
                data_obj.topoplot3d = qt.Plot3D(spatial_data, 
                                                name = 'topography')
            data_obj.create_file()
    else:
        spatial_data = False

    #Set up tasks:
    zactask = tc.AcOutTask(DEV, ACZCHAN, SAMPLES, SAMPLE_RATE, sync = True)    
    zacdata = GenSineWave(SAMPLES, AMPLITUDE, PHASE)
    zactask.set_signal(zacdata)
    if feedback:
        ztask = tc.DcOutTask(DEV, [DCCHANS[2]])
    else:
        ztask = False
    maintask = mimCallbackTask(afm, TOPOCHAN, SAMPLES, SAMPLE_RATE, feedback,
                                ztask, approach_data, spatial_data, zstart)
    
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
	    	if maintask.userin == 'u':
				maintask.z[0] += Z_STEP
	    	elif maintask.userin == 'd':
        		maintask.z[0] -= Z_STEP
	    	getattr(daq, 'set_ao%i' % DCZCHAN[0])(maintask.z[0])


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
    def __init__(self, afm, channels, samples, samplerate, feedback, ztask,
                 approach_data, patial_data, zstart):
        tc.AnalogInCallbackTask.__init__(self, dev, channels, samples, samplerate)
        self.z = np.zeros(1)
        self.z[0] = start_position[2]
        self.feedback = feedback
        self.userin = False # stores whether task has been interrupted by user
        self.callcounter = 0 # iterates each approach curve (callback)
        self.ztask = ztask
        self.approach_data = approach_data
        self.fastvec = fastvec
        self.slowvec = slowvec
        self.xytask = xytask
        self.spatial_data = spatial_data
        self.repeat = repeat

    def EveryNCallback(self):
        tc.AnalogInCallbackTask.EveryNCallback(self)      
        if self.callcounter % len(self.fastvec) == 0:
            slowcounter = self.callcounter/len(self.fastvec)
            if slowcounter >= len(self.slowvec) and not self.repeat:
                self.userin = 'q'
                print 'Completed scan'
                self.StopTask()
                return
            else:
                #For non-yscan casese the following line will remain yvec[0]
                self.slowpos = self.slowvec[slowcounter % len(self.slowvec)]
                print('slow direction set to') 
                print(self.slowvec[slowcounter % len(self.slowvec)])
                if self.spatial_data:
                    for data_obj in self.spatial_data:
                        data_obj.new_block()
                        
        #separate out channels in most recent trace
        odata = self.data
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
        else: # if not in feedback mode
            pass
        #Check that Z is within limits
        if self.z[0] > Z_MAX:
            self.z[0] = Z_MAX
            logging.warning('Reached maximum allowable value: %f' % Z_MAX)
        elif self.z[0] < Z_MIN:
            self.z[0] = Z_MIN
            logging.warning('Reached minimum allowable value: %f' % Z_MIN)
        #Calculate xy poitionfor ext point
        fastindex = self.callcounter % len(self.fastvec)
        fastpos = self.fastvec[fastindex]
        xy = fastpos + self.slowpos
        #Finally, do the deed:
        self.xytask.set_voltage(xy)
        #record raw data:
        self.approach_data.add_data_point(
                np.arange(SAMPLES), odata, self.z * 1e3 * np.ones(SAMPLES))
        self.approach_data.new_block()

        #If we're in the right-going segmet
        if len(self.fastvec) > 2:
            if fastindex < len(self.fastvec)/2:
                spatial_data_current = self.spatial_data[0] #right
            else:
                spatial_data_current = self.spatial_data[1] #left
            spatial_data_current.add_data_point(
                xy[0] * 1e3, xy[1] * 1e3, self.z[0] * 1e3)
        
        self.callcounter += 1
        print self.z[0]
