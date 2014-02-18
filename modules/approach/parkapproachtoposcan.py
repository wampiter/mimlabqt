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

def measure(feedback=False, xvec=np.zeros(1), yvec=np.zeros(1), 
            angle=0.0):
    '''
    Takes a topo measurement. Aborts with 'q' (or stop)
    in non-feedback, 'u' and 'd' will move stage up and down (backwards)
    
    Kwargs:
    feedback -- whether or not to use feedback
    xvec -- fast scan vector (numpy array)for 1D and 2D scans, false for 0D
    yvec -- slow scan vector (numpy array)for 2D scans, false for 1D and 0D
    angle -- measurement ange (in radians)
    '''
    #Check that xvec is specified if yvec is
    if len(yvec) > 1 and not len(xvec) > 1:
        logging.error('Cannot have slow but not fast scan vector')
        raise ValueError(len(xvec),len(yvec))
    #check dimensions of scan
    scanx=False;scany=False
    if len(yvec) > 1:
        scany = True
    if len(xvec) > 1:
        scanx = True
        
    #set up qt meaurement flow
    qt.flow.connect('measurement-end', tc.kill_all_tasks)
    qt.mstart()
    #Reset instrument (or initialize if not already)
    if qt.instruments.get_instruments_by_type('NI_DAQ'):
        daq = qt.instruments.get_instruments_by_type('NI_DAQ')[0]
        daq.reset()
    else:
        daq = qt.instruments.create('daq', 'NI_DAQ', id = 'Dev1')

    #Store current position of scanner
    start_position = np.zeros(3)    
    for i, chan in enumerate(DCCHANS):
        start_position[i] = getattr(daq,'get_ao%i' % chan)()
        if np.isnan(start_position[i]):
            daq.set('ao%i' % chan, 0.0)
            start_position[i] = 0.0 
    #if we're not already at the right xy position, lift z and move there:
    if (abs(start_position[0] - xvec[0]) > 0.001 or 
            abs(start_position[1] - yvec[0]) > 0.001):
        mg.ramp(daq, 'ao %i' % DCCHANS[2], start_position[2] + Z_LIFT)     
        mg.ramp(daq, 'ao%i' % DCCHANS[0], xvec[0])
        mg.ramp(daq, 'ao%i' % DCCHANS[1], yvec[0])
        mg.ramp(daq, 'ao %i' % DCCHANS[2], start_position[2])
        
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
    if scanx:
        spatial_data_right = qt.Data(name = 'spatial_data_right')
        spatial_data_left = qt.Data(name = 'spatial_data_left')
        spatial_data = [spatial_data_right, spatial_data_left]
        for data_obj in spatial_data:
            data_obj.add_coordinate('x [mV]')
            data_obj.add_coordinate('y [mV]')
            data_obj.add_value('z [mV]')
            data_obj.topoplot2d = qt.Plot2D(
                    spatial_data, name = 'topography linecuts')
            if scany:
                data_obj.topoplot3d = qt.Plot3D(spatial_data, 
                                                name = 'topography')
            data_obj.create_file()
    else:
        spatial_data = False
    
    # Scan out and back
    xvec = np.append(xvec,xvec[::-1]) 
    fastvec = np.zeros((len(xvec),2))
    slowvec = np.zeros((len(yvec),2))
    for i, x in enumerate(xvec):
        fastvec[i] = np.array(
                [np.cos(angle) * xvec[i], np.sin(angle) * xvec[i]])
    for i, y in enumerate(yvec):
        slowvec[i] = np.array(
                [-np.sin(angle) * yvec[i], np.cos(angle) * yvec[i]])

    #Set up tasks:
    zactask = tc.AcOutTask(DEV, ACZCHAN, SAMPLES, SAMPLE_RATE, sync = True)    
    zacdata = GenSineWave(SAMPLES, AMPLITUDE, PHASE)
    zactask.set_signal(zacdata)
    if feedback:
        ztask = tc.DcOutTask(DEV, [DCCHANS[2]])
    else:
        ztask = False
    xytask = tc.DcOutTask(DEV, DCCHANS[0:2])
    maintask = mimCallbackTask(DEV, TOPOCHAN, SAMPLES, SAMPLE_RATE, feedback,
                                ztask, approach_data, fastvec, slowvec, xytask, start_position, spatial_data)
    
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
	    getattr(daq, 'set_ao%i' % DCCHANS[2])(maintask.z[0])


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
    def __init__(self, dev, channels, samples, samplerate, feedback, ztask,
                 approach_data, fastvec, slowvec, xytask, start_position, spatial_data):
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

    def EveryNCallback(self):
        tc.AnalogInCallbackTask.EveryNCallback(self)      
        if self.callcounter % len(self.fastvec) == 0:
            slowcounter = self.callcounter/len(self.fastvec)
            if slowcounter >= len(self.slowvec) and not repeat:
                self.userin = 'q'
                print 'Completed scan'
                self.StopTask()
                return
            else:
                #For non-yscan casese the following line will remain yvec[0]
                self.slow_pos = self.slowvec[slowcounter % len(self.slowvec)]
                print('slow direction set to %f volts.' 
                      % self.slowvec[slowcounter % len(self.slowvec)])
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
        xytask.set_voltage(xy)
        #record raw data:
        self.approach_data.add_data_point(
                np.arange(SAMPLES), odata, self.z * 1e3 * np.ones(SAMPLES))
        self.approach_data.new_block()

        #If we're in the right-going segmet
        if len(self.fastvec) > 1:
            if fastindex < len(self.fastvec)/2:
                spatial_data_current = self.spatial_data[0] #right
            else:
                spatial_data_current = self.spatial_data[1] #left
            spatial_data_current.add_data_point(
                self.fastvec[self.callcounter] * 1e3, self.y * 1e3, self.z * 1e3)#NEEDS FIX
        
        self.callcounter += 1
        print self.z[0]
