class simulations():
 
 	def __init__(self):

		self.Tcmb = 2.73 #K

		self.degrees2radians = np.pi / 180.
		self.arcmins2radians = self.degrees2radians / 60.

	def fn_get_lxly_az_angle(self,lx,ly): #get szimuthal angle here

		return 2*np.arctan2(lx, -ly)
		#return 2.*np.arctan2(ly, lx)

	def is_seq(self,o):
		return hasattr(o, '__len__')

	def is_seq_of_seq(self,o):
		if not self.is_seq(o):
			return False
		for s in o:
			if not self.is_seq(s):
				return False
			return True

	def get_lxly(self,mapparams): #get lx, ly

		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians

		lx, ly = np.meshgrid( np.fft.fftfreq( nx, dx ), np.fft.fftfreq( ny, dy ) )
		lx *= 2* np.pi
		ly *= 2* np.pi

		return lx, ly

	def Cls2CLS(self,Cls,mapparams): #convert Cls1d to CLS 2D

		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians

		# find how many different Cls were passed
		noofrows = len(Cls)
		if self.is_seq_of_seq(Cls):
			noofcols = np.shape(Cls)[1]
			Cls_tot = noofcols - 1 #first is els, then Cls
		else:
			els = np.arange(noofrows)
			Cls = np.array( [els,Cls] ).T
			Cls_tot = 1

		lx, ly = self.get_lxly(mapparams)
		L = np.sqrt(lx**2. + ly**2.)

		# processing Cls now
		CLS = np.zeros( (Cls_tot,L.shape[0],L.shape[1]) )
		for clcnt in range(Cls_tot):
			CLS[clcnt,:,:] = np.interp(L.flatten(), Cls[:,0], abs(Cls[:,clcnt+1]), right = 0.).reshape(L.shape) 
		return CLS

	def Dls2map(self, Dls, mapparams, passing_Dls = 1, CMB_outputscale = 1):

		#takes Dls as input and converts into CAMBMAPFFT
		###########################################################################
		#first check if only T or P is supplied as well; also obtain els
		noofrows = len(Dls)
		#get Dls in the correct format: noofels x noof_power_spectra (TT/EE/BB etc.)
		if self.is_seq_of_seq(Dls):
			els = np.asarray(Dls[:,0])
			noofcols = np.shape(Dls)[1]
			Dls_tot = noofcols - 1 #first is els, then Dls
		else:
			els = np.arange(2,noofrows)
			Dls = [els,Dls]

		###########################################################################
		#Dls to Cls
		if passing_Dls:
			if CMB_outputscale == 1:
				Cls = ( self.Tcmb**2. * Dls * 2 * np.pi ) / ( els[:,None] * (els[:,None] + 1) )
			else:
				Cls = ( Dls * 2 * np.pi ) / ( els[:,None] * (els[:,None] + 1) )
			Cls[:,0] = els
		else:
			Cls = Dls

		###########################################################################
		#Cls2map
		nx, ny, dx, dy = mapparams
		dx *= self.arcmins2radians
		dy *= self.arcmins2radians

		#scalefac = np.sqrt((nx * ny) / (dx * dy))
		scalefac = np.sqrt(1./ (dx * dy))

		CLS = self.Cls2CLS(Cls,mapparams)
		self.CLS = CLS #just store for future use
		self.scalefac = scalefac

		CLS = np.sqrt(CLS)
		for ccc in range(len(CLS)):
			CLS[ccc][np.isnan(CLS[ccc])] = 0.

		lx, ly = self.get_lxly(mapparams)
		angle = self.fn_get_lxly_az_angle(lx,ly)
		########################################################################
		########################################################################
		# processing Cls now
		if Dls_tot>1:
			TFFT, EFFT, BFFT, TEFFT = CLS
		else:
			TFFT = CLS

		if Dls_tot>1:
			FORQ = ( np.cos(angle) * EFFT - np.sin(angle) * BFFT )
			FORU = ( np.sin(angle) * EFFT + np.cos(angle) * BFFT )

			CAMBMAP_FFT = np.asarray([TFFT, EFFT, BFFT, TEFFT]) * scalefac
		else:
			CAMBMAP_FFT = np.asarray(TFFT) * scalefac

		self.CAMBMAP_FFT = CAMBMAP_FFT
		###########################################################################
		########################################################################
		
	def fn_make_cmb_sims(self, mapparams, noofsims = 100, in_uk = 1):

		CAMBMAP_FFT = self.CAMBMAP_FFT
		nx, ny, dx, dy = mapparams
		lx, ly = self.get_lxly(mapparams)
		angle = self.fn_get_lxly_az_angle(lx,ly)

		SIMS = []
		for n in range(noofsims):

			#make Gaussian realisation first
			GAUSS_REALS = np.asarray( [np.fft.fft2( np.random.randn(nx,ny) ) for teb in range(len(CAMBMAP_FFT) - 1)] )


			#add power based on CAMB Cls now
			tqulen = 5 #T,E,B,Q,U
			CMB = np.zeros( (tqulen, nx, ny) )
			for teb in range( len(GAUSS_REALS) ): #T,E,B
				if teb == 1: #for E include correlation between T
					t1 = GAUSS_REALS[0] * self.CLS[3] / self.CLS[0]**0.5
					t2 = GAUSS_REALS[1] * ( self.CLS[1] - (self.CLS[3]**2. /self.CLS[0]) )**0.5
					Emap_fft = (t1 + t2) * self.scalefac
					Emap_fft[np.isnan(Emap_fft)] = 0.
					CMB[teb] = np.fft.ifft2( Emap_fft ).real
				else:
					CMB[teb] = np.fft.ifft2( np.copy( CAMBMAP_FFT[teb] ) * GAUSS_REALS[teb] ).real

				#subtract mean from maps
				CMB[teb] = CMB[teb] - np.mean(CMB[teb])

			if tqulen>1:
				E_FFT, B_FFT = np.fft.fft2(CMB[1]),np.fft.fft2(CMB[2])
				CMB[3] = np.fft.ifft2( np.cos(angle) * E_FFT - np.sin(angle) * B_FFT ).real #Q
				CMB[4] = np.fft.ifft2( np.sin(angle) * E_FFT + np.cos(angle) * B_FFT ).real #U

			SIMS.append(CMB)

		SIMS = np.asarray(SIMS)

		if in_uk:
			SIMS *= 1e6

		return SIMS


########################################################################
########################################################################
########################################################################

import numpy as np, os, glob, os, sys, argparse
from pylab import *
from cosmologydigitisation.Polarisation.cosmdigitclasses import *

cosm = Cosmologist()
sims = simulations()



############################################################
############################################################

parser = argparse.ArgumentParser(description='')
parser.add_argument('-howmanysims', dest='howmanysims', action='store', help='howmanysims', type=int, default=10) #def 5
parser.add_argument('-boxsize_am', dest='boxsize_am', action='store', help='boxsize in arcmins', type=float, default=512.) #def 100
parser.add_argument('-dx', dest='dx', action='store', help='dx - pixel size in arcmins', type=int, default=1.0) #def 0.5
parser.add_argument('-Dlfile_len', dest='Dlfile_len', action='store', help='CAMB lensed Dls file', type=str, default='output_planck_r_0.0_2015_cosmo_lensedCls.dat')
parser.add_argument('-passing_Dls', dest='passing_Dls', action='store', help='passing Dls or Cls?', type=int, default=1)
parser.add_argument('-CMB_outputscale', dest='CMB_outputscale', action='store', help='CMB outputscale from CAMB - multiply Dls by T_cmb', type=int, default=1)
parser.add_argument('-in_uk', dest='in_uk', action='store', help='CMB map units in uK', type=int, default=1)

args = parser.parse_args()
args_keys = args.__dict__
for kargs in args_keys:
        param_value = args_keys[kargs]
        if isinstance(param_value, str):
                cmd = '%s = "%s"' %(kargs, param_value)
        else:
                cmd = '%s = %s' %(kargs, param_value)
        exec(cmd)

############################################################
############################################################
Dls_len = np.loadtxt(Dlfile_len,usecols=[0,1,2,3,4])
l = np.asarray(Dls_len[:,0])

nx = int( boxsize_am/dx )
mapparams = [nx, nx, dx, dx] #square!

#initialise CAMB params for sims
sims.Dls2map(Dls_len, mapparams, passing_Dls = passing_Dls, CMB_outputscale = CMB_outputscale)
CMB_SIMS = sims.fn_make_cmb_sims(mapparams, noofsims = howmanysims, in_uk = in_uk) #in uK


showmaps = False
showps = True

if showmaps:
	#show map plots now
	titles = ['T', 'E', 'B', 'Q', 'U']
	for n in range(howmanysims):
		CURRSIM = CMB_SIMS[n]
		for iii in range(len(CURRSIM)):
			ax = subplot(1,len(CURRSIM), iii + 1)
			imshow(CURRSIM[iii]);colorbar()
			title(titles[iii])
		show()
	#sys.exit()

print("plot ps")
if showps:
	# calculate and show ps
	titles = ['T', 'E', 'B', 'Q', 'U']
	cols = ["k", "r", "b", "g", "y"]
	ps = [0.0 for m in titles]
	for n in range(howmanysims):
		print(n)
		CURRSIM = CMB_SIMS[n]
		k = 0
		for iii in range(1):#len(CURRSIM)
			k, p, perr, h = cosm.sriniPowerSpectrum(mapparams, CURRSIM[iii])
			p = p * dx * dx * arcmins2radians * arcmins2radians * k * ( k + 1.0 ) / (2.0 * pi)
			ps[iii] += p
	ps = [ps[j]/float(howmanysims) for j in range(len(titles))]
	for iii in range(1):
		plot(k, ps[iii], label = titles[iii], color = "r")
	for iii in range(1):
		plot(l, (sims.Tcmb**2.0) * Dls_len[:, iii + 1] * 10**12, color = "k")

	legend()
	xscale("linear")
	yscale("linear")
	xlim((amin(k), 2500))
	legend(loc = "lower left")
	show()