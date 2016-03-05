# A 165 LINE TOPOLOGY OPTIMIZATION CODE BY NIELS AAGE AND VILLADS EGEDE JOHANSEN, JANUARY 2013
# EXTENDED BY A COUPLE OF LINES FOR LOADING LOAD- AND SUPPORT VIA IMAGES BY MATHIAS BERNHARD, DECEMBER 2015
from __future__ import division
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from matplotlib import colors
import matplotlib.pyplot as plt
from skimage import io
# MAIN DRIVER
def main(volfrac,penal,rmin,ft,chg,fldr):
	print("Minimum compliance problem with OC")

	# load images first to automatically set nelx and nely
	samples = ['mbb','short_cantilever','beam','hook','chair']
	folder = fldr+'/' #samples[0]+'/'
	try:
		sup = io.imread(folder+'support.png')
		loa = io.imread(folder+'load.png')
	except Exception as e:
		print 'image missing'
		exit()

	nelx = sup.shape[1]
	nely = sup.shape[0]

	# ------------------------------------------
	# define passive (void or solid) regions
	# ------------------------------------------
	passive = np.zeros(nely*nelx)
	pas = None
	try:
		pas = io.imread(folder+'passive.png')
	except Exception as e:
		pass
	if not pas is None:
		for i in xrange(pas.shape[0]):
			for j in xrange(pas.shape[1]):
				if pas[i,j,0]>5:
					passive[j*nely+i] = 1
				if pas[i,j,1]>5:
					passive[j*nely+i] = 2


	print("nodes: " + str(nelx) + " x " + str(nely))
	print("volfrac: " + str(volfrac) + ", rmin: " + str(rmin) + ", penal: " + str(penal))
	print("Filter method: " + ["Sensitivity based","Density based"][ft])
	# user defined loop parameters
	tol = chg
	maxloops = 150
	# Max and min stiffness
	Emin=1e-9 # Young's modulus for void-like material, not zero to avoid singularity
	Emax=1.0  # Young's modulus for solid material
	# nu: poisson's ratio
	# dofs:
	ndof = 2*(nelx+1)*(nely+1) # number of degrees of freedom
	# Allocate design variables (as array), initialize and allocate sens.
	x=volfrac * np.ones(nely*nelx,dtype=float)
	xold=x.copy()
	xPhys=x.copy()
	g=0 # must be initialized to use the NGuyen/Paulino OC approach
	#dc=np.zeros((nely,nelx), dtype=float)
	# FE: Build the index vectors for the for coo matrix format.
	KE=lk()
	edofMat=np.zeros((nelx*nely,8),dtype=int)
	for elx in range(nelx):
		for ely in range(nely):
			el = ely+elx*nely
			n1=(nely+1)*elx+ely
			n2=(nely+1)*(elx+1)+ely
			edofMat[el,:]=np.array([2*n1+2, 2*n1+3, 2*n2+2, 2*n2+3,2*n2, 2*n2+1, 2*n1, 2*n1+1])
	# Construct the index pointers for the coo format
	iK = np.kron(edofMat,np.ones((8,1))).flatten()
	jK = np.kron(edofMat,np.ones((1,8))).flatten()
	# Filter: Build (and assemble) the index+data vectors for the coo matrix format
	nfilter=nelx*nely*((2*(np.ceil(rmin)-1)+1)**2)
	iH = np.zeros(nfilter)
	jH = np.zeros(nfilter)
	sH = np.zeros(nfilter)
	cc=0
	for i in range(nelx):
		for j in range(nely):
			row=i*nely+j
			kk1=int(np.maximum(i-(np.ceil(rmin)-1),0))
			kk2=int(np.minimum(i+np.ceil(rmin),nelx))
			ll1=int(np.maximum(j-(np.ceil(rmin)-1),0))
			ll2=int(np.minimum(j+np.ceil(rmin),nely))
			for k in range(kk1,kk2):
				for l in range(ll1,ll2):
					col=k*nely+l
					fac=rmin-np.sqrt(((i-k)*(i-k)+(j-l)*(j-l)))
					iH[cc]=row
					jH[cc]=col
					sH[cc]=np.maximum(0.0,fac)
					cc=cc+1
	# Finalize assembly and convert to csc format
	H=coo_matrix((sH,(iH,jH)),shape=(nelx*nely,nelx*nely)).tocsc()
	Hs=H.sum(1)

	# ------------------------------------------
	# boundary conditions and support
	# even numbers: freedom in x
	# odd numbers: freedom in y
	# ------------------------------------------
	dofs=np.arange(2*(nelx+1)*(nely+1))
	fixed = []
	for i in xrange(sup.shape[0]):
		for j in xrange(sup.shape[1]):
			if sup[i,j,0]>2:
				fixed.append(2*(j*(nely+1)+i)+1)
			if sup[i,j,1]>2:
				fixed.append(2*(j*(nely+1)+i))

	free=np.setdiff1d(dofs,fixed)

	# Solution and RHS (right-hand-side) vectors
	f=np.zeros((ndof,1)) # forces
	u=np.zeros((ndof,1)) # deformation matrix

	# ------------------------------------------
	# Set load
	# ------------------------------------------
	#f[1,0]=1 # MBB beam
	for i in xrange(loa.shape[0]):
		for j in xrange(loa.shape[1]):
			vertl = loa[i,j,0]-128
			horil = loa[i,j,1]-128
			if abs(vertl)>2:
				f[2*(j*(nely+1)+i)+1] = vertl/128.0
			if abs(horil)>2:
				f[2*(j*(nely+1)+i)] = horil/128.0

	# Initialize plot and plot the initial design
	plt.ion() # Ensure that redrawing is possible
	fig,ax = plt.subplots()
	im = ax.imshow(-xPhys.reshape((nelx,nely)).T, cmap='gray',\
	interpolation='nearest',norm=colors.Normalize(vmin=-1,vmax=0))
	fig.show()
   	# Set loop counter and gradient vectors
	loop=0
	change=1
	dv = np.ones(nely*nelx)
	dc = np.ones(nely*nelx)
	ce = np.ones(nely*nelx)

	while change>tol and loop<maxloops:
		loop=loop+1
		# Setup and solve FE problem
		sK=((KE.flatten()[np.newaxis]).T*(Emin+(xPhys)**penal*(Emax-Emin))).flatten(order='F')
		K = coo_matrix((sK,(iK,jK)),shape=(ndof,ndof)).tocsc()
		# Remove constrained dofs from matrix
		K = K[free,:][:,free]
		# Solve system
		u[free,0]=spsolve(K,f[free,0])
		# Objective and sensitivity
		ce[:] = (np.dot(u[edofMat].reshape(nelx*nely,8),KE) * u[edofMat].reshape(nelx*nely,8) ).sum(1)
		obj=( (Emin+xPhys**penal*(Emax-Emin))*ce ).sum()
		dc[:]=(-penal*xPhys**(penal-1)*(Emax-Emin))*ce
		dv[:] = np.ones(nely*nelx)
		# Sensitivity filtering:
		if ft==0:
			dc[:] = np.asarray((H*(x*dc))[np.newaxis].T/Hs)[:,0] / np.maximum(0.001,x)
		elif ft==1:
			dc[:] = np.asarray(H*(dc[np.newaxis].T/Hs))[:,0]
			dv[:] = np.asarray(H*(dv[np.newaxis].T/Hs))[:,0]
		# Optimality criteria
		xold[:]=x
		(x[:],g)=oc(nelx,nely,x,volfrac,dc,dv,g)
		# Filter design variables
		if ft==0:   xPhys[:]=x
		elif ft==1:	xPhys[:]=np.asarray(H*x[np.newaxis].T/Hs)[:,0]

		# passive elements
		xPhys[passive==1] = 0
		xPhys[passive==2] = 1

		# Compute the change by the inf. norm
		change=np.linalg.norm(x.reshape(nelx*nely,1)-xold.reshape(nelx*nely,1),np.inf)
		# Plot to screen
		im.set_array(-xPhys.reshape((nelx,nely)).T)
		# try to visualize u-matrix (deformation?)
		#im.set_array(u[::2].reshape((nelx+1,nely+1)).T)
		fig.canvas.draw()
		# Write iteration history to screen (req. Python 2.6 or newer)
		print("it.: {0} , obj.: {1:.3f} Vol.: {2:.3f}, ch.: {3:.3f}".format(\
					loop,obj,(g+volfrac*nelx*nely)/(nelx*nely),change))
	# Make sure the plot stays and that the shell remains
	plt.show()
	raw_input("Press any key...")
#element stiffness matrix
def lk():
	E=1.0
	nu=0.3
	k=np.array([1/2-nu/6,1/8+nu/8,-1/4-nu/12,-1/8+3*nu/8,-1/4+nu/12,-1/8-nu/8,nu/6,1/8-3*nu/8])
	KE = E/(1-nu**2)*np.array([
	[k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
	[k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
	[k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
	[k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
	[k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
	[k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
	[k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
	[k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]] ]);
	return (KE)
# Optimality criterion
def oc(nelx,nely,x,volfrac,dc,dv,g):
	l1=0
	l2=1e9
	move=0.2
	# reshape to perform vector operations
	xnew=np.zeros(nelx*nely)
	while (l2-l1)/(l1+l2)>1e-3:
		lmid=0.5*(l2+l1) # Lagrange multiplier
		xnew[:]= np.maximum(0.0,np.maximum(x-move,np.minimum(1.0,np.minimum(x+move,x*np.sqrt(-dc/dv/lmid)))))
		gt=g+np.sum((dv*(xnew-x)))
		if gt>0 :
			l1=lmid
		else:
			l2=lmid
	return (xnew,gt)
# The real main driver
if __name__ == "__main__":
	# Default input parameters
	volfrac=0.4
	rmin=5.4 # lower number: more branching (initial: 5.4, try 2.0 or 1.5) proposal: 0.04 * nelx
	penal=3.0 # ensure black and white solution
	ft=1 # ft==0 -> sensitivity filtering, ft==1 -> density filtering
	chg=0.1
	folder = 'mbb'
	import sys
	if len(sys.argv)>1: volfrac=float(sys.argv[1])
	if len(sys.argv)>2: rmin   =float(sys.argv[2])
	if len(sys.argv)>3: chg   =float(sys.argv[3])
	if len(sys.argv)>4: folder  =str(sys.argv[4])
	if len(sys.argv)>5: penal  =float(sys.argv[5])
	if len(sys.argv)>6: ft     =int(sys.argv[6])
	main(volfrac,penal,rmin,ft,chg,folder)
