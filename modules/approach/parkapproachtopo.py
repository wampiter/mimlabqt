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

def measure(feedback=False):
    #set up qt meaurement flow
    qt.flow.connect('measurement-end', tc.kill_all_tasks)
    qt.mstart()
    #Reset instrument (or initialize if not already)
    if qt.instruments.get_instruments_by_type('NI_DAQ'):
        daq = qt.instruments.get_instruments_by_type('NI_DAQ')[0]
        daq.reset()
    else:
        daq = qt.instruments.create('daq', 'NI_DAQ', id = 'Dev1')
              
    #Prepare approach data structure to store absolutely all raw data:
    approach_data = qt.Data(name='approach_curves')
    approach_data.add_coordinate('sample')
    approach_data.add_value('optical')
    approach_data.add_value('z [mV]')

    approach_plot = qt.Plot2D(approach_data, name = 'approach_trace')

    approach_data.create_file()
    approach_data.copy_file('parkapproachtopo.py')
    approach_data.copy_file('approachimagingparamspark.py')

    #get initial position:
    zstart = np.zeros(1)
    zstart[0] = getattr(daq,'get_ao%i' % DCCHANS[2])()
    if not zstart[0]:
        zstart[0] = 0.0
    #Set up tasks:
    zactask = tc.AcOutTask(DEV, ACZCHAN, SAMPLES, SAMPLE_RATE, sync = True)    
    zacdata = GenSineWave(SAMPLES, AMPLITUDE, PHASE)
    zactask.set_signal(zacdata)
    if feedback:
        ztask = tc.DcOutTask(DEV, [DCCHANS[2]])
    else:
        ztask = False
    maintask = mimCallbackTask(DEV, TOPOCHAN, SAMPLES, SAMPLE_RATE, zstart, feedback,
                                ztask, approach_data)
    
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
    def __init__(self, dev, channels, samples, samplerate, zstart, feedback, ztask,
                 approach_data):
        tc.AnalogInCallbackTask.__init__(self, dev, channels, samples, samplerate)
        self.z = zstart
        self.feedback = feedback
        self.userin = False # stores whether task has been interrupted by user
        self.callcounter = 0 # iterates each approach curve (callback)
        self.ztask = ztask
        self.approach_data = approach_data

    def EveryNCallback(self):
        tc.AnalogInCallbackTask.EveryNCallback(self)      
                    
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
            
        #record raw data:
        self.approach_data.add_data_point(
                np.arange(SAMPLES), odata, self.z * 1e3 * np.ones(SAMPLES))
        self.approach_data.new_block()
            
        self.callcounter += 1
        print self.z[0]
