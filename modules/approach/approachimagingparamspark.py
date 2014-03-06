import numpy as np

#OUTPUT CHANNEL NUMBERS
DEV = 1
DCZCHAN = [3] #list of DC channel numbers [x,y,z]
ACZCHAN = [0] #list of AC Z channel

#INPUT CHANNELS
MIMCHANS = [4,5] #Anaog input channels [C,R]
TOPOCHAN = [0] #Analog input for laser feedback
XYCHANS = [1,2]

SAMPLES = 200 #Samples per approach curve
SAMPLE_RATE = 2.0e3
AMPLITUDE = .8
PHASE = np.pi

Z_STEP = 2.0e-2 # %f User controlled step in non-feedback (approach) mode
Z_SCANNER_STEP = 2.0e-2

Z_MAX = 2.0 #Maximum Z voltage. VERY IMPORTANT TO NOT CRASH TIP
Z_MIN = -3.0  #Minimum Z voltage. Prevent runaway in other (less bad) direction
Z_LIFT = 0.01

CONTACT_CENTER = 55
CONTACT_SPAN = 10
EQ_CENTER = 10
EQ_SPAN = 10

LOOP_GAIN = 0.01
TARGET_DISPLACE = 0.05

CONTACT_START = CONTACT_CENTER - CONTACT_SPAN/2
CONTACT_STOP = CONTACT_CENTER + CONTACT_SPAN/2
EQ_START = EQ_CENTER - EQ_SPAN/2
EQ_STOP = EQ_CENTER + EQ_SPAN/2


def GenSineWave(elements, amplitude , phase):
    wave = np.zeros(elements)
    for i in np.arange(elements):
        wave[i] = amplitude * np.sin(phase+(2*np.pi*i/elements))
    return wave
