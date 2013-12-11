import numpy as np
import msvcrt
import qt

from approachimagingparamspark import * #used like C constants, always in CAPS
import taskclasses as tc
import measurement_general as mg


def measure(feedback=False, xvec=np.zeros(1), yvec=np.zeros(1), 
            angle=0.0, repeat=False):
    '''
    Takes an appraoch curve topo only measurement with laser feedback. Aborts with 'q' (or stop)
    in non-feedback, 'u' and 'd' will move stage up and down
    
    Kwargs:
    feedback -- whether or not to use feedback
    xvec -- fast scan vector (numpy array)for 1D and 2D scans, false for 0D
    yvec -- slow scan vector (numpy array)for 2D scans, false for 1D and 0D
    angle -- measurement ange (in radians)
    repeat -- whether to continuously repeat measurement until 'q' is pressed
    '''
    #Check that xvec is specified if yvec is
    if len(yvec) > 1 and not len(xvec) > 1:
        logging.error('Cannot have slow but not fast scan vector')
        raise ValueError(len(xvec),len(yvec))
    reload(approachimagingparams)
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
        start_position[i] = daq.get_parameters['ao%i' % chan]['value']
        if np.isnan(start_position[i]):
            daq.set('ao%i' % chan, 0.0)
            start_position[i] = 0.0 
    #if we're not already at the right xy position, lift z and move there:
    if (abs(start_position[0] - xvec(0)) > 0.001 or 
            abs(start_position[1] - yvec(0)) > 0.001):
        mg.ramp(daq, 'ao %i' % DCCHANS[2], start_position[2] + Z_LIFT)     
        mg.ramp(daq, 'ao%i' % DCCHANS[0], xvec[0])
        mg.ramp(daq, 'ao%i' % DCCHANS[1], yvec[0])
        mg.ramp(daq, 'ao %i' % DCCHANS[2], start_position[2])
        
    #Prepare approach data structure to store absolutely all raw data:
    approach_data = qt.Data(name='approach_curves')
    approach_data.add_coordinate('sample')
    approach_data.add_coordinate('x [mV]')
    approach_data.add_coordinate('y [mV]')
    approach_data.add_value('z [mV]')
    approach_data.add_value('MIM-C [V]')
    approach_data.add_value('MIM-R [V]')
    approach_data.mimcplot = qt.Plot2D(approach_data, name='MIM-C [V]', 
                                       coorddim = 0, valdim = 3)
    appraoch_data.mimrplot = qt.Plot2D(approach_data, name='MIM-R [V]',
                                       coorddim = 0, valdim = 4)
    approach_data.create_file()
    approach_data.copy_file('appraochimaging.py')
    approach_data.copy_file('approachimagingparams.py')
    #Prepare spatial data structures for "normal" scan data:
    if scanx:
        spatial_data_right = qt.Data(name = 'spatial_data_right')
        spatial_data_left = qt.Data(name = 'spatial_data_left')
        spatial_data = [spatial_data_right, spatial_data_left]
        for data_obj in spatial_data:
            data_obj.add_coordinate('x [mV]')
            data_obj.add_coordinate('y [mV]')
            data_obj.add_value('z [mV]')
            data_obj.add_value('MIM-C [V]')
            data_obj.add_value('MIM-R [V]')
            data_obj.topoplot2d = qt.Plot2D(
                    spatial_data, name = 'topography linecuts')
            data_obj.mimcplot2d = qt.Plot2D(
                    spatial_data, name = 'MIM-C linecuts', valdim = 3)
            data_obj.mimrplot2d = qt.Plot2D(
                    spatial_data, name = 'MIM-R linecuts', valdim = 4)
            if scany:
                data_obj.topoplot3d = qt.Plot3D(spatial_data, 
                                                name = 'topography')
                data_obj.mimcplot3d = qt.Plot3D(spatial_data, 
                                                name = 'MIM-C', valdim = 3)
                data_obj.mimrplot3d = qt.Plot3D(spatial_data, name = 'MIM-R', 
                                                valdim = 4)
            data_obj.create_file()
    
    # Scan out and back
    xvec = np.append(xvec,xvec[::-1]) 
    fastvec = np.zeros((len(xvec),2))
    slowvec = np.zeros((len(yvec),2))
    for i, x in enumerate(xvec):
        fastvec[i] = np.array(
                [np.cos(angle) * xvec[i], np.sine(angle) * xvec[i]])
    for i, y in enumerate(yvec):
        slowvec[i] = np.array(
                [-np.sine(angle) * yvec[i], np.cos(angle) * yvec[i]])

    #Set up tasks:
    zactask = tc.AcOutTask(DEV, ACZCHAN, SAMPLES, SAMPLERATE, sync = True)    
    zacdata = GenSineWave(DEV, samples,amplitude,phase)
    zactask.set_signal(zacdata)
    xytask = tc.DcOutTask(DEV, DCCHAN[0:2])
    ztask = tc.DcOutTask(DEV, [DCCHAN[2]])
    maintask = mimCalbackTask(DEV, MIMCHANS, SAMPLES, SAMPLERATE, zstart, feedback)
    
    maintask.userin = False
    while maintask.userin != 'q':
        if msvcrt.kbhit():
            maintask.userin = msvcrt.getch()
        try:
            qt.msleep(.5)
        except:
            logging.warning('Live plotting glitch (No big deal)')

    print('%i approaches completed' % maintask.callcounter)
    #stop tasks before error    
    qt.mend()
    qt.msleep(.1)
    #close all 3 data objects    
    approach_data.close_file()
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
    def __init__(self, channels, samples, samplerate, zstart, feedback):
        tc.AnalogInCallbackTask.__init__(self, channels, samples, samplerate)
        self.z = zstart
        self.feedback = feedback
        self.userin = False # stores whether task has been interrupted by user
        self.callcounter = 0 # iterates each approach curve (callback)

    def EveryNCallback(self):
        tc.AnalogInCallbackTask.EveryNCallback(self)      
        if self.callcounter % len(xvec) == 0:
            slowcounter = self.callcounter/len(fastvec)
            if slowcounter >= len(yvec) and not repeat:
                self.userin = 'q'
                print 'Completed scan'
                self.StopTask()
                return
            else:
                #For non-yscan casese the following line will remain yvec[0]
                self.slow_pos = slowvec[slowcounter % len(slowvec)]
                print('slow direction set to %f volts.' 
                      % yvec[slowcounter % len(yvec)])
                for data_obj in spatial_data:
                    data_obj.new_block()
                    
        #separate out channels in most recent trace
        cdata = self.data[0:SAMPLES]
        rdata = self.data[SAMPLES:2*SAMPLES]
        if self.feedback:
            #take derivative and smooth:
            self.datadiff = np.diff(mg.movingaverage(cdata,INNER_WINDOW)) 
            self.datadiff = mg.movingaverage(self.datadiff,OUTER_WINDOW)
            self.datadiff[0:OUTER_WINDOW/2] = 0 #kill ragged edges
            self.datadiff[len(self.datadiff) - OUTER_WINDOW/2
                          : len(self.datadiff)] = 0#..
            #Sample at which contact occurs            
            minarg = np.argmin(self.datadif) 
            #increment z proportional to offset of contact to make negfeedback:
            if (minarg > LOWER_SAMPLE_THRESHOLD
                and minarg < UPPER_SAMPLE_THRESHOLD #ensure minarg in window
                and np.amin(self.datad)<MAGNITUDE_THRESHOLD):
                    self.z += (minarg - CONTACT_SAMPLE)*FEEDBACK_GAIN
        else: # if not in feedback mode
            if self.userin == 'u':
                self.z += Z_STEP #step up
                self.userin = False
            elif self.userin == 'd':
                self.z -= Z_STEP #step down
                self.userin = False 
        #Check that Z is within limits
        if self.z > Z_MAX:
            self.z = Z_MAX
            logging.warning('Reached maximum allowable value: %f' % Z_MAX)
        elif self.z < Z_MIN:
            self.z = Z_MIN
            logging.warning('Reached minimum allowable value: %f' % Z_MIN)
        #Calculate xy poitionfor ext point
        fastindex = self.callcounter % len(fastvec)
        fastpos = fastvec[fastindex]
        xy = fastpos + self.slowpos
        #Finally, do the deed:
        xytask.set_voltage(xy)
        ztask.set_voltage(self.z)
        
        #record raw data:
        try:
            approach_data.add_data_point(
                    arange(SAMPLES), xy[0] * 1e3 * ones(SAMPLES), 
                    xy[1] * 1e3 * ones(SAMPLES), self.z * 1e3 * ones(SAMPLES), 
                    cdata, rdata)
            approach_data.new_block()
        except:
            logging.warning('Failed to record approach curve')
        #record processed data:
        try:
            mimCabs = np.mean(cdata[FAR_FIRST_SAMP:FAR_LAST_SAMP])\
                - np.mean(cdata[CLOSE_FIRST_SAMP:CLOSE_LAST_SAMP])
            mimRabs = np.mean(rdata[FAR_FIRST_SAMP:FAR_LAST_SAMP])\
                - np.mean(rdata[CLOSE_FIRST_SAMP:CLOSE_LAST_SAMP])
            
            #If we're in the right-going segmet
            if fastindex < len(fastvec)/2:
                spatial_data_current = spatial_data_right
            else:
                spatial_data_current = spatial_data_left
                
            spatial_data_current.add_data_point(
                    xvec[self.callcounter] * 1e3, self.y * 1e3, self.z * 1e3,
                    mimCabs, mimRabs)
        except:
            logging.warning('Failed to record processed point data')
            
        self.callcounter += 1
