import numpy as np

#OUTPUT CHANNEL NUMBERS
DEV = 2
DCCHANS = [1,2,3] #list of DC channel numbers [x,y,z]
ACZCHAN = [0] #list of AC Z channel
MIMCHANS = [4,5] #Anaog input channels [C,R]
TOPOCHAN = [0] #Analog input for laser feedback

SAMPLES = 200 #Samples per approach curve
SAMPLE_RATE = 2.0e4
AMPLITUDE = .7
PHASE = np.pi

Z_STEP = 2.0e-2 # %f User controlled step in non-feedback (approach) mode

Z_MAX = 1.0 #Maximum Z voltage. VERY IMPORTANT TO NOT CRASH TIP
Z_MIN = -0.5  #Minimum Z voltage. Prevent runaway in other (less bad) direction
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
