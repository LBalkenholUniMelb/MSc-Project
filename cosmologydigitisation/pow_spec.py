def fn_radial_profile(Z, XY, noofbins = 100, bin_size = None, minbin = 0., maxbin = 10.):

	"""
	calculates the radial profile of power spectrum or any other quantity
	in case of power spectrum, the outputs are k, p_k, p_k_err
	"""

	Z = np.asarray(Z)
	if XY == None:
		Y, X = np.indices(Z.shape)
	else:
		X, Y = XY

	RADIUS = (X**2. + Y**2.) ** 0.5

	# Overwriting Maxbin Here!
	maxbin = amax(RADIUS)

	if bin_size == None:
		binarr=np.linspace(minbin,maxbin,noofbins)
	else:
		binarr=np.arange(minbin,maxbin,bin_size)

	#print(binarr)

	radprf=np.zeros((len(binarr),4))

	hit_count=[]

	for b,bin in enumerate(binarr):
		ind=np.where((RADIUS>=bin) & (RADIUS<bin+bin_size))
		radprf[b,0]=(bin+bin_size/2.)
		hits = len(np.where(abs(Z[ind])>0.)[0])

		if hits>0:
			radprf[b,1]=np.sum(Z[ind])/hits
			radprf[b,2]=np.std(Z[ind])
		hit_count.append(hits)

	hit_count=np.asarray(hit_count)
	std_mean=np.sum(radprf[:,2]*hit_count)/np.sum(hit_count)
	errval=std_mean/(hit_count)**0.5
	radprf[:,2]=errval
	radprf[:,3]=hit_count

	#print("pow")
	#print(sum(Z))
	#print(sum(radprf[:,1]*hit_count))

	return radprf

def fn_get_lxly(mapparams):

	"""
	mapparams: [nx, ny, dx, dy]: nx , ny = dimension of image; dx, dy = pixel resolution in arcmins
	just computes the fourier wave vector
	"""

	nx, ny, dx, dy = mapparams
	dx *= arcmins2radians
	dy *= arcmins2radians

	lx, ly = np.meshgrid( np.fft.fftfreq( nx, dx ), np.fft.fftfreq( ny, dy ) )

	#fourier wavevector k = 2*pi/lambda
	lx *= 2* np.pi
	ly *= 2* np.pi

	return lx, ly

def fn_plot_pow_spec(mapparams, MAP1, MAP2 = None, binsize = None):

	"""
	computes 2d power spectrum of image and then the radial bining of it. Returns the 1d power spectrum.
	mapparams: [nx, ny, dx, dy]: nx , ny = dimension of image; dx, dy = pixel resolution in arcmins
	MAP1: first image
	MAP2: second image: If none, then compute auto power spectrum of MAP1. else cross power specctrum of MAP1 and MAP2
	binsize: for 1d power spectrum
	"""
	nx, ny, dx, dy = mapparams
	dx_rad = dx * arcmins2radians

	if MAP2 == None: #compute auto power spectra
		cosm = cosmask([nx, ny])
		w = (float(1.0/(nx*ny))*sum(cosm*cosm))
		#print("PS REAL")
		#print(sum(MAP1*MAP1)*w)
		#print(sum(MAP1*cosm*MAP1*cosm))
		#print(sum(zeropad(MAP1*cosm)*zeropad(MAP1*cosm)))
		# MAP_F = np.fft.fft2(MAP1*cosm)
		MAP_F = np.fft.fft2(zeropad(MAP1*cosm))
		MAP_PSD = ( MAP_F*conjugate(MAP_F) )/(4.0 * nx * ny * w)

		#abs( np.fft.fft2(MAP1) * dx_rad)** 2 / (nx * ny)
		# Do zero padding and apply mask here

	else: #compute cross power spectra between 2 maps
		#subplot(121);imshow(MAP1);colorbar();subplot(122);imshow(MAP2);colorbar();show();sys.exit()
		MAP_PSD = np.fft.fft2(MAP1) * dx_rad * np.conj( np.fft.fft2(MAP2) ) * dx_rad / (nx * ny)

	lx, ly = fn_get_lxly([nx*2, ny*2, dx, dy])
	if binsize == None:
		binsize = (lx.ravel()[1] - lx.ravel()[0])  # * 10
	# if np.max(lx)>1e5: binsize *= 2 #just increasing binsize

	#subplot(121);imshow(np.fft.fftshift(MAP_PSD.real));title("2D transform");colorbar();show();
	pow_spec_1d = fn_radial_profile(MAP_PSD, (lx,ly), bin_size = binsize, minbin = np.min(abs(lx)), maxbin = np.max(lx))

	return pow_spec_1d #contains k, p_k, p_k_err, hitcount


################################################################################
################################################################################
################################################################################

################################################################################
################################################################################
import numpy as np, pickle, gzip
from matplotlib.pyplot import *
from numpy import savetxt, sum, conjugate, shape, amax, asarray, cos, pi, append, tile, transpose, zeros, real

arcmins2radians = np.radians(1./60.)

def cosmask(dimensions):
    # get cos array for y
    p_y = asarray(range(int(dimensions[0]/2)))
    cosy = cos(pi*p_y/dimensions[0])
    cosy = append(cosy[::-1], cosy)
    cosy = transpose(tile(cosy, (dimensions[1], 1)))
    # get cos array for x
    p_x = asarray(range(int(dimensions[1]/2)))
    cosx = cos(pi*p_x/dimensions[1])
    cosx = append(cosx[::-1], cosx)
    cosx = tile(cosx, (dimensions[0], 1))
    # return mask
    mask = cosy*cosx
    return mask


def zeropad(M):
    dim = shape(M)
    M_pad = zeros((2*asarray(dim)))
    M_pad[0:dim[0], 0:dim[1]] = M
    return M_pad