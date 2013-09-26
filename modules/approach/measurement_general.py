import qt
import logging
import numpy as np

def ramp(instrument, parameter, target_value, valstep=False, timestep=0.01):
    '''
    Ramps parameter on insrument from its current value to target_value
    one valstep every timestep.
    '''
    start_value = instrument.get('parameter')
    if not start_value:
        start_value = instrumet.get_parameters()[parameter]['value']
        if start_value:
            print 'so using last set value of %f' % start_value
        else:
            raise TypeError("and this parameter as not yet been set")
        
    if not valstep:
        valstep = (target_value - start_value)/100.0
    vec = np.arange(start_value, target_value,
                 abs(valstep) * sign(target_value - start_value))
    for value in vec:
        instrument.set(parameter, value)
        qt.msleep(timestep)
        
def movingaverage(interval, window_size):
    '''
    Moving average of interval. Returns array of ame size as interval,
    though you may want to ignore data withing window_size/2 of either end.
    '''
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')