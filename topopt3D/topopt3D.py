# A 165 LINE TOPOLOGY OPTIMIZATION CODE BY NIELS AAGE AND VILLADS EGEDE JOHANSEN, JANUARY 2013
from __future__ import division
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from matplotlib import colors
import matplotlib.pyplot as plt
# MAIN DRIVER
def main(nelx,nely,nelz,volfrac,penal,rmin,ft):
	print("Minimum compliance problem with OC")
	print("nodes: " + str(nelx) + " x " + str(nely) + " x " + str(nelz))
	print("volfrac: " + str(volfrac) + ", rmin: " + str(rmin) + ", penal: " + str(penal))
	print("Filter method: " + ["Sensitivity based","Density based"][ft])
	# user defined loop parameters
	tol = 0.1
	maxloops = 150
	# Max and min stiffness
	Emin=1e-9 # Young's modulus for void-like material, not zero to avoid singularity
	Emax=1.0  # Young's modulus for solid material
	# total number of elements
	nele = nelx * nely * nelz
	# dofs:
	ndof = 3*(nelx+1)*(nely+1)*(nelz+1) # number of degrees of freedom
	# Allocate design variables (as array), initialize and allocate sens.
	x=volfrac * np.ones(nele,dtype=float)
	xold=x.copy()
	xPhys=x.copy()
	g=0 # must be initialized to use the NGuyen/Paulino OC approach
	#dc=np.zeros((nely,nelx,nelz), dtype=float) # not necessary, assigned in loop
	# FE: Build the index vectors for the for coo matrix format.
	KE=lk()
	edofMat=np.zeros((nele,24),dtype=int)
	xylayer = (nelx+1)*(nely+1)
	for elz in range(nelz):
		for elx in range(nelx):
			for ely in range(nely):
				el = elz*(nelx*nely) + elx*nely + ely
				#n1=(nely+1)*elx+ely
				#n2=(nely+1)*(elx+1)+ely
				#edofMat[el,:]=np.array([2*n1+2, 2*n1+3, 2*n2+2, 2*n2+3,2*n2, 2*n2+1, 2*n1, 2*n1+1])
				#n1 = elz*(nelx+1)*(nely+1) + elx*(nely+1) + ely
				#n2 = (elz+1)*(nelx+1)*(nely+1) + elx*(nely+1) + ely
				#n3 = elz*(nelx+1)*(nely+1) + (elx+1)*(nely+1) + ely
				nodes = [elz*xylayer+elx*(nely+1)+ely+1,
					   elz*xylayer+(elx+1)*(nely+1)+ely+1,
					   elz*xylayer+(elx+1)*nely+ely+1,
					   elz*xylayer+(elx+1)*nely+ely,
					   (elz+1)*xylayer+elx*(nely+1)+ely+1,
	   				   (elz+1)*xylayer+(elx+1)*(nely+1)+ely+1,
					   (elz+1)*xylayer+(elx+1)*nely+ely+1,
					   (elz+1)*xylayer+(elx+1)*nely+ely]
				elmat = np.zeros(24)
				cnt = 0
				for i in nodes:
					elmat[cnt] = i*3
					cnt+=1
					elmat[cnt] = i*3+1
					cnt+=1
					elmat[cnt] = i*3+2
					cnt+=1
				edofMat[el,:]=elmat
	# Construct the index pointers for the coo format
	iK = np.kron(edofMat,np.ones((24,1))).flatten()
	jK = np.kron(edofMat,np.ones((1,24))).flatten()
	# Filter: Build (and assemble) the index+data vectors for the coo matrix format
	nfilter=nele*((2*(np.ceil(rmin)-1)+1)**3)
	iH = np.zeros(nfilter)
	jH = np.zeros(nfilter)
	sH = np.zeros(nfilter)
	cc=0
	for k1 in xrange(nelz):
		for i1 in xrange(nelx):
			for j1 in xrange(nely):
				el1=k1*nelx*nely + i1*nely + j1
				ii1=int(np.maximum(i1-(np.ceil(rmin)-1),0))
				ii2=int(np.minimum(i1+np.ceil(rmin),nelx))
				jj1=int(np.maximum(j1-(np.ceil(rmin)-1),0))
				jj2=int(np.minimum(j1+np.ceil(rmin),nely))
				kk1=int(np.maximum(k1-(np.ceil(rmin)-1),0))
				kk2=int(np.minimum(k1+np.ceil(rmin),nelz))
				for k2 in xrange(kk1,kk2):
					for i2 in xrange(ii1,ii2):
						for j2 in xrange(jj1,jj2):
							el2=k2*nelx*nely + i2*nely + j2
							fac=rmin-np.sqrt((i1-i2)**2+(j1-j2)**2+(k1-k2)**2)
							iH[cc]=el1
							jH[cc]=el2
							sH[cc]=np.maximum(0.0,fac)
							cc=cc+1
	# Finalize assembly and convert to csc format
	#H=coo_matrix((sH,(iH,jH)),shape=(nelx*nely,nelx*nely)).tocsc()
	H=coo_matrix((sH,(iH,jH)),shape=(nele,nele)).tocsc()
	Hs=H.sum(1)

	# ------------------------------------------
	# boundary conditions and support
	# even numbers: freedom in x
	# odd numbers: freedom in y
	# ------------------------------------------
	dofs=np.arange(ndof)
	#p1 = dofs[0:2*(nely+1):2]
	#p1 = np.array([2*(nely-35),2*35*(nely+1)-1])
	#p2 = np.array([2*(nelx+1)*(nely+1)-1])
	#p2 = np.array([dofs.shape[0]-30,dofs.shape[0]-31])
	#fixed=np.union1d(p1,p2)

	left = []
	for ely in xrange(nely):
		for elz in xrange(nelz):
			el = elz*(nelx+1)*(nely+1) + ely
			left.append(3*el)
			left.append(3*el+1)
			left.append(3*el+2)
	#for elz in xrange(nelz+1):
	#	el = elz * (nelx+1)*(nely+1)
	#	left.append(3*el+1)
	fixed = np.array(left)
	free=np.setdiff1d(dofs,fixed) #Return the sorted, unique values in `ar1` that are not in `ar2`.

	#fixed = np.array([0,1,2,11085,11086,11087])
	#free=np.setdiff1d(dofs,fixed)

	# Solution and RHS (right-hand-side) vectors
	f=np.zeros((ndof,1)) # forces
	u=np.zeros((ndof,1)) # deformation matrix

	# ------------------------------------------
	# Set load
	# ------------------------------------------
	#f[1,0]=1 # MBB beam
	loadindex = 3*((nelz/2)*(nelx+1)*(nely+1)) +1 # +1 is for y-direction
	f[loadindex,0] = 1

	# Initialize plot and plot the initial design
	plt.ion() # Ensure that redrawing is possible
	fig,ax = plt.subplots()
	im = ax.imshow(-xPhys[:nelx*nely].reshape((nelx,nely)).T, cmap='gray',\
	interpolation='nearest',norm=colors.Normalize(vmin=-1,vmax=0))
	fig.show()
   	# Set loop counter and gradient vectors
	loop=0
	change=1
	dv = np.ones(nele)
	dc = np.ones(nele)
	ce = np.ones(nele)

	while change>tol and loop<maxloops:
		loop=loop+1
		# Setup and solve FE problem
		sK=((KE.flatten()[np.newaxis]).T*(Emin+(xPhys)**penal*(Emax-Emin))).flatten(order='F')
		K = coo_matrix((sK,(iK,jK)),shape=(ndof,ndof)).tocsc()
		# Remove constrained dofs from matrix
		K = K[free,:][:,free]
		# Solve system
		u[free,0]=spsolve(K,f[free,0]) # error: matrix is exactly singular
		# Objective and sensitivity
		ce[:] = (np.dot(u[edofMat].reshape(nele,24),KE) * u[edofMat].reshape(nele,24) ).sum(1)
		obj=( (Emin+xPhys**penal*(Emax-Emin))*ce ).sum()
		dc[:]=(-penal*xPhys**(penal-1)*(Emax-Emin))*ce
		dv[:] = np.ones(nele)
		# Sensitivity filtering:
		if ft==0:
			dc[:] = np.asarray((H*(x*dc))[np.newaxis].T/Hs)[:,0] / np.maximum(0.001,x)
		elif ft==1:
			dc[:] = np.asarray(H*(dc[np.newaxis].T/Hs))[:,0]
			dv[:] = np.asarray(H*(dv[np.newaxis].T/Hs))[:,0]
		# Optimality criteria
		xold[:]=x
		(x[:],g)=oc(nelx,nely,nelz,x,volfrac,dc,dv,g)
		# Filter design variables
		if ft==0:   xPhys[:]=x
		elif ft==1:	xPhys[:]=np.asarray(H*x[np.newaxis].T/Hs)[:,0]

		# Compute the change by the inf. norm
		change=np.linalg.norm(x.reshape(nele,1)-xold.reshape(nele,1),np.inf)

		# plot must be made new!
		# Plot to screen
		#im.set_array(-xPhys.reshape((nelx,nely)).T)
		fig.canvas.draw()
		# Write iteration history to screen (req. Python 2.6 or newer)
		print("it.: {0} , obj.: {1:.3f} Vol.: {2:.3f}, ch.: {3:.3f}".format(\
					loop,obj,(g+volfrac*nele)/(nele),change))
	# Make sure the plot stays and that the shell remains
	plt.show()
	raw_input("Press any key...")

#element stiffness matrix
def lk():
	E=1
	nu=0.3
	A = np.array([[32,6,-8,6,-6,4,3,-6,-10,3,-3,-3,-4,-8],
	              [-48,0,0,-24,24,0,0,0,12,-12,0,12,12,12]])
	#k = 1.0/144.0*A.T*np.array([1,nu]).T
	k = np.dot(1.0/144.0*A.T,np.array([1,nu]).T)
	K1 = np.array([[k[0],k[1],k[1],k[2],k[4],k[4]],
	      [k[1],k[0],k[1],k[3],k[5],k[6]],
	      [k[1],k[1],k[0],k[3],k[6],k[5]],
	      [k[2],k[3],k[3],k[0],k[7],k[7]],
	      [k[4],k[5],k[6],k[7],k[0],k[1]],
	      [k[4],k[6],k[5],k[7],k[1],k[0]]])
	K2 = np.array([[k[8], k[7], k[11],k[5], k[3], k[6]],
	      [k[7], k[8], k[11],k[4], k[2], k[4]],
	      [k[9],k[9],k[12],k[6], k[3], k[5]],
	      [k[5], k[4], k[10],k[8], k[1], k[9]],
	      [k[3], k[2], k[4], k[1], k[8], k[11]],
	      [k[10],k[3], k[5], k[11],k[9],k[12]]])
	K3 = np.array([[k[5], k[6], k[3], k[8], k[11],k[7]],
	      [k[6], k[5], k[3], k[9],k[12],k[9]],
	      [k[4], k[4], k[2], k[7], k[11],k[8]],
	      [k[8], k[9],k[1], k[5], k[10],k[4]],
	      [k[11],k[12],k[9],k[10],k[5], k[3]],
	      [k[1], k[11],k[8], k[3], k[4], k[2]]])
	K4 = np.array([[k[13],k[10],k[10],k[12],k[9],k[9]],
	      [k[10],k[13],k[10],k[11],k[8], k[7]],
	      [k[10],k[10],k[13],k[11],k[7], k[8]],
	      [k[12],k[11],k[11],k[13],k[6], k[6]],
	      [k[9],k[8], k[7], k[6], k[13],k[10]],
	      [k[9],k[7], k[8], k[6], k[10],k[13]]])
	K5 = np.array([[k[0],k[1], k[7], k[2],k[4], k[3]],
	      [k[1],k[0], k[7], k[3],k[5], k[10]],
	      [k[7],k[7], k[0], k[4],k[10],k[5]],
	      [k[2],k[3], k[4], k[0],k[7], k[1]],
	      [k[4],k[5], k[10],k[7],k[0], k[7]],
	      [k[3],k[10],k[5], k[1],k[7], k[0]]])
	K6 = np.array([[k[13],k[10],k[6], k[12],k[9],k[11]],
	      [k[10],k[13],k[6], k[11],k[8], k[1]],
	      [k[6], k[6], k[13],k[9],k[1], k[8]],
	      [k[12],k[11],k[9],k[13],k[6], k[10]],
	      [k[9],k[8], k[1], k[6], k[13],k[6]],
	      [k[11],k[1], k[8], k[10],k[6], k[13]]])
	KE = np.zeros((24,24))
	order = [K1, K2, K3, K4,
	 		 K2.T, K5, K6, K3.T,
	 	 	 K3.T, K6, K5.T, K2.T,
	 	 	 K4, K3, K2, K1.T]
	ind = 0
	for i in xrange(4):
		for j in xrange(4):
			KE[i*6:i*6+6,j*6:j*6+6] = order[ind]
			ind+=1
	KE = 1/((nu+1)*(1-2*nu))*KE
	'''KE = 1/((nu+1)*(1-2*nu))* \
		 np.array([[K1, K2, K3, K4],
		  [K2.T, K5, K6, K3.T],
		  [K3.T, K6, K5.T, K2.T],
		  [K4, K3, K2, K1.T]])'''
	return (KE)

# Optimality criterion
def oc(nelx,nely,nelz,x,volfrac,dc,dv,g):
	l1=0
	l2=1e9
	move=0.2
	# reshape to perform vector operations
	xnew=np.zeros(nelx*nely*nelz)
	cnt = 0
	while (l2-l1)/(l1+l2)>1e-3:
		lmid=0.5*(l2+l1) # Lagrange multiplier
		xnew[:]= np.maximum(0.0,np.maximum(x-move,np.minimum(1.0,np.minimum(x+move,x*np.sqrt(-dc/dv/lmid)))))
		gt=g+np.sum((dv*(xnew-x)))
		if gt>0 :
			l1=lmid
		else:
			l2=lmid
		if l1==l2:
			l2 = 1e9
	return (xnew,gt)

# The real main driver
if __name__ == "__main__":
	# Default input parameters
	nelx=20
	nely=15
	nelz=10
	volfrac=0.4
	rmin=1.4 #deeper: more branching (initial: 5.4, try 2.0!) proposal: 0.04 * nelx
	penal=3.0 # ensure black and white solution
	ft=0 # ft==0 -> sensitivity filtering, ft==1 -> density filtering
	import sys
	if len(sys.argv)>1: nelx   =int(sys.argv[1])
	if len(sys.argv)>2: nely   =int(sys.argv[2])
	if len(sys.argv)>3: nelz   =int(sys.argv[3])
	if len(sys.argv)>4: volfrac=float(sys.argv[4])
	if len(sys.argv)>5: rmin   =float(sys.argv[5])
	if len(sys.argv)>6: penal  =float(sys.argv[6])
	if len(sys.argv)>7: ft     =int(sys.argv[7])
	main(nelx,nely,nelz,volfrac,penal,rmin,ft)
