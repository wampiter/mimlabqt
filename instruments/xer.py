from instrument import Instrument
import types
import logging
import numpy
import ctypes as c
import qt

class ScanConfigParam(c.Structure):
    _fields_ = [("nWidth",c.c_int),
                ("nHeight",c.c_int),
                ("fOverScan",c.c_float),
                ("bSineScan",c.c_int),
                ("bDetectorDriven",c.c_int)]

class ImageParam(c.Structure):
    _fields_ = [("fSizeX",c.c_float),
                ("fSizeY",c.c_float),
                ("fOffsetX",c.c_float),
                ("fOffsetY",c.c_float),
                ("fRotation",c.c_float),
                ("fRate",c.c_float),
                ("bSwap",c.c_int)]

class xer(Instrument):
    '''
    '''

    def __init__(self, name):
        '''
        '''
        # Initialize wrapper functions
        logging.info('Initializing instrument xer')
        Instrument.__init__(self, name, tags=['physical'])

	#Create all important object:
	self.raw = c.cdll.LoadLibrary("d:\\Dropbox\\Users\\Scott\\XER\\bin\\XERemote.dll")
	self.raw.AddWndListener()
	self.raw.Connect()
	self.scp = ScanConfigParam()
	self.ip = ImageParam()

        # Add parameters to wrapper
        self.add_parameter('nWidth',
            flags=Instrument.FLAG_GETSET,
            units='', type=types.IntType)

	self.add_parameter('nHeight',
	    flags=Instrument.FLAG_GETSET,
            units='', type=types.IntType)

	self.add_parameter('fOverScan',
	    flags=Instrument.FLAG_GETSET,
	    units='', type=types.FloatType)

	self.add_parameter('fRate',
	    flags=Instrument.FLAG_GETSET,
	    units='Hz', type=types.FloatType)

	self.add_parameter('fSizeX',
	    flags=Instrument.FLAG_GETSET,
	    units='um', type=types.FloatType)

	self.add_parameter('fSizeY',
	    flags=Instrument.FLAG_GETSET,
	    units='um', type=types.FloatType)

	self.add_parameter('fOffsetX',
	    flags=Instrument.FLAG_GETSET,
	    units='um', type=types.FloatType)

	self.add_parameter('fOffsetY',
	    flags=Instrument.FLAG_GETSET,
	    units='um', type=types.FloatType)

	self.add_parameter('fRotation',
	    flags=Instrument.FLAG_GETSET,
	    units='degrees', type=types.FloatType)

	self.scp.nWidth = c.c_int(256)
	self.scp.nHeight = c.c_int(256)
        self.scp.fOverScan = c.c_float(2)
        self.ip.fRate = c.c_float(2)
        self.ip.fSizeX = c.c_float(4.0)
        self.ip.fSizeY = c.c_float(4.0)
        self.raw.ScanConfig(c.byref(self.scp))
        self.raw.ImageConfig(c.byref(self.ip))

# --------------------------------------
#           functions
# --------------------------------------

    def XYScannerMove(self, x, y):
        self.raw.XYScannerMove(c.c_double(x),c.c_double(y))

    def ZScannerMove(self, z):
        self.raw.ZScannerMove(c.c_double(z))

    def Approach(self):
        self.raw.Approach()

    def LineScan(self, line, wait):
        self.raw.LineScan(c.c_int(line), c.c_bool(wait))

    def Image(self, wait):
        self.raw.Image(c.byref(self.ip), c.c_bool(wait))

    def Abort(self):
        self.raw.Abort()


# --------------------------------------
#           parameters
# --------------------------------------
    def do_get_nWidth(self):
	return int(self.scp.nWidth)
	
    def do_set_nWidth(self, width):
        self.scp.nWidth = c.c_int(width)
        self.raw.ScanConfig(c.byref(self.scp))

    def do_get_nHeight(self):
        return int(self.scp.nHeight)

    def do_set_nHeight(self, height):
        self.scp.nHeight = c.c_int(height)
        self.raw.ScanConfig(c.byref(self.scp))

    def do_get_fOverScan(self):
        return float(self.scp.fOverScan)

    def do_set_fOverScan(self, f):
        self.scp.fOverScan = c.c_float(f)
        self.raw.ScanConfig(c.byref(self.scp))

    def do_get_fRate(self):
        return float(self.ip.fRate)

    def do_set_fRate(self, f):
        self.ip.fRate = c.c_float(f)
        self.raw.ImageConfig(c.byref(self.ip))

    def do_get_fSizeX(self):
        return float(self.ip.fSizeX)

    def do_set_fSizeX(self, l):
        self.ip.fSizeX = c.c_float(l)
        self.raw.ImageConfig(c.byref(self.ip))

    def do_get_fSizeY(self):
        return float(self.ip.fSizeY)

    def do_set_fSizeY(self, l):
        self.ip.fSizeY = c.c_float(l)
        self.raw.ImageConfig(c.byref(self.ip))

    def do_get_fOffsetX(self):
        return float(self.ip.fOffsetX)

    def do_set_fOffsetX(self, l):
        self.ip.fOffsetX = c.c_float(l)
        self.raw.ImageConfig(c.byref(self.ip))

    def do_get_fOffsetY(self):
        return float(self.ip.fOffsetY)

    def do_set_fOffsetY(self, l):
        self.ip.fOffsetY = c.c_float(l)
        self.raw.ImageConfig(c.byref(self.ip))

    def do_get_fRotation(self):
        return float(self.ip.fRotation)

    def do_set_fRotation(self, d):
        self.ip.fRotation = c.c_float(d)
        self.raw.ImageConfig(c.byref(self.ip))
