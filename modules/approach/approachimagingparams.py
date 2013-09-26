import numpy as np

#OUTPUT CHANNEL NUMBERS
DCCHANS = [1,2,0] #list of DC channel numbers [x,y,z]
ACZCHAN = 3
MIMCHANS = [4,5] #Anaog input channels [C,R]

SAMPLES = 200 #Samples per approach curve
SAMPLE_RATE = 2.0e4
AMPLITUDE = 2.5
PHASE = np.pi

INNER_WINDOW = 12 # %i Smoothing window used on raw input data
OUTER_WINDOW = 12 # %i Smoothing window used on derivative of smoothed data

LOWER_SAMPLE_THRESHOLD = 95 # %i bottom of window in which feedback occurs
UPPER_SAMPLE_THRESHOLD = 145 # %i top of window in which feedback occurs
MAGNITUDE_THRESHOLD = -1.0e-3 # %f max amplitude of the min of the derivative
                                #at which feedback occurs
CONTACT_SAMPLE = 135 # %i Sample at which feedback holds point of contact
FEEDBACK_GAIN = 0.8e-4 #%f feedback gain in V/sample

Z_STEP = 2.0e-3 # %f User controlled step in non-feedback (approach) mode

Z_MAX = 0.35 #Maximum Z voltage. VERY IMPORTANT TO NOT CRASH TIP
Z_MIN = 0.0  #Minimum Z voltage. Prevent runaway in other (less bad) direction
Z_LIFT = 0.01

FAR_FIRST_SAMP = 60 #Defines window for far MIM point
FAR_LAST_SAMP = 80
CLOSE_FIRST_SAMP = 160 #Defines window for close MIM point
CLOSE_LAST_SAMP = 180

def GenSineWave(elements, amplitude , phase):
    wave = np.zeros(elements)
    for i in np.arange(elements):
        wave[i] = amplitude * np.sin(phase+(2*np.pi*i/elements))
    return wave