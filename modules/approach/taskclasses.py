from PyDAQmx import Task
from PyDAQmx.DAQmxConstants import *
from PyDAQmx.DAQmxTypes import *
import numpy as np
import qt

#I assume that the device is 'Dev2'.

#This variable is a list of all the current task objects.
running_tasks = []

def kill_all_tasks():
    '''
    Stop and Clear all tasks
    '''
    for task in running_tasks:
        try:
            task.StopTask()
        except:
            logging.warning('Task already stopped')
        task._ClearTaskOld()
    del running_tasks[:]

qt.flow.connect('measurement-end',kill_all_tasks)
    
class BaseTask(Task):
    def __init__(self):
        Task.__init__(self)
        running_tasks.append(self)
        
    def ClearTask(self):
        '''Clear task from DAQ and running_tasks'''
        Task.ClearTask(self)
        del running_tasks[running_tasks.index(self)]
        
    def _ClearTaskOld(self):
        Task.ClearTask(self)

    def kill_task(self):
        '''Stop and Clear task'''
        self.StopTask()
        self.ClearTask()
        

class DcOutTask(BaseTask):
    '''
    Creates a task for outputting DC Voltages ASAP.
    Each write (set_voltage) operation will take ~1ms PER CHANNEL
    
    Args:
        Channels: List of input channels (ints). Order matters!
    '''
    def __init__(self, dev, channels, start_voltage = False):
        if len(channels) < 1:
            logging.error('Channels must be list with at least one channel')
            raise ValueError(len(channels))
        BaseTask.__init__(self)
        for chan in channels:
            self.CreateAOVoltageChan("Dev%i/ao%i" % (dev,chan), "", -10.0, 10.0,
                                     DAQmx_Val_Volts, "")
        self._voltage = [False]
        
    def set_voltage(self,voltage):
        '''
        Set output voltage
        
        Args:
            voltage: numpy array of voltges
        '''
        self.WriteAnalogF64(1, 1, 10.0, DAQmx_Val_GroupByChannel, voltage,
                            None,None)
        self._voltage = voltage
        
    def get_voltage(self):
        '''Get current output voltge'''
        if self._voltage[0]:
            return self._voltage
        else: 
            logging.error('Voltage not yet definied for this DcOutTask.')
            
class AcOutTask(BaseTask):
    '''
    Create a task for outputting an AC signal periodicaly.
    
    Args:
        channels: list of output channels
        samples: samples in the signal
        samplerate: sampling rate
    
    Kwargs:
        sync: synchronize with analog input
    '''
    def __init__(self, dev, channels, samples, samplerate, sync = False):
        if len(channels) < 1:
            logging.error('Channels must be list with at least one channel')
            raise ValueError(len(channels))
        BaseTask.__init__(self)
        for chan in channels:
            self.CreateAOVoltageChan("Dev%i/ao%i" % (dev,chan), "", -10.0, 10.0,
                                     DAQmx_Val_Volts, None)
        self.CfgSampClkTiming("",samplerate,DAQmx_Val_Rising,
                                DAQmx_Val_ContSamps, samples)
        self.samples = samples
        if sync:
            self.CfgDigEdgeStartTrig('/Dev%i/ai/StartTrigger' % dev, 
                                       DAQmx_Val_Rising)
        self._signal = [False]
        
    def set_signal(self, signal):
        '''
        Sets Signal to be outputted.
        
        Args:
            signal: numpy array of float64 values of length samples * channels
            Group by channel
        '''
        self.WriteAnalogF64(self.samples, False, 10.0, 
                              DAQmx_Val_GroupByChannel, signal, None, None)
        self._signal = signal
        
    def get_signal(self):
        if self._signal.any():
            return self._signal
        else:
            logging.error('no signal definied')
            
def AnologInCallbackTask(BaseTask):
    '''
    Creates a continuous analog read task which after set of samples overwrites
    its data to self.data. Other things can be done on callback by adding to
    EveryNCallback(self)
    
    Args:
        chan: list of input channels
        samples: number of samples
        samplerate: sampling rate
    '''
    def __init__(self, dev, channels, samples, samplerate):
        BaseTask.__init__(self)
        self.samples = samples
        self.samplerate = samplerate
        self.data_samples = samples * len(channels)
        self.data = np.zeros(data_samples)
        for chan in channels:
            self.CreateAIVoltageChan("Dev%i/ai%i" % (dev,chan), "", DAQmx_Val_RSE, 
                                     -10.0, 10.0, DAQmx_Val_Volts, None)
        self.CfgSampClkTiming("", self.samplerate, DAQmx_Val_Rising, 
                              DAQmx_Val_ContSamps, self.samples)
        self.AutoRegisterEveryNSamplesEvent(DAQmx_Val_Acquired_Into_Buffer, 
                                            self.samples, 0)
        self.AutoRegisterDoneEvent(0)
        
    def EveryNCallback(self):
        read = int32()
        self.ReadAnalogF64(self.samples, 10.0, DAQmx_Val_GroupByChannel, 
                           self.data, self.data_samples, byref(read), None)
                           
    def DoneCallback(self, status):
        print "Status",status.value
        return 0
