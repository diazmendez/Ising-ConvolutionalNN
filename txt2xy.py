"""
Reads the ising2d.c output files and writes them in the CNN input form
"""

import numpy as np
import sys

l=int(sys.argv[1])
path=sys.argv[2]
bname=sys.argv[3]


# upload and formating
#################################################

xx=np.loadtxt(path+'/snapshots.dat'); # inputs
xx=np.reshape(xx,(-1, 1, l, l))

yy=np.loadtxt(path+'/targets.dat'); # targets

ff=np.loadtxt(path+'/lasts.dat'); # final confs
ff=np.reshape(ff,(-1, l, l))

xx=(xx+1)/2  # normalizing the input between 0 and 1



# balancing the data
##################################

classes=yy.shape[1]

unbalance=np.zeros(classes)

for i in range(classes):
	unbalance[i]=int(yy[:,i].sum())
 
unbalance=unbalance-np.amin(unbalance)
unbalance=unbalance.astype(int)

#d0=int(yy[:,0].sum()) 
#d1=int(yy[:,1].sum())

for i in range (classes):
	for j in range(unbalance[i]):
		dice=int(np.random.rand()*yy.shape[0])
		while (yy[dice,i]==0):
			dice=int(np.random.rand()*yy.shape[0])	
		yy=np.delete(yy,dice,0)
		xx=np.delete(xx,dice,0)
		ff=np.delete(ff,dice,0)

"""
if (d0>d1): # delete elements from the first category
	for i in range(d0-d1):
		dice=int(np.random.rand()*yy.shape[0])
		while (yy[dice,0]==0):
			dice=int(np.random.rand()*yy.shape[0])	
		yy=np.delete(yy,dice,0)
		xx=np.delete(xx,dice,0)
		ff=np.delete(ff,dice,0)
if (d1>d0): # delete elements from the second category
	for i in range(d1-d0):
		dice=int(np.random.rand()*yy.shape[0])
		while (yy[dice,1]==0):
			dice=int(np.random.rand()*yy.shape[0])	
		yy=np.delete(yy,dice,0)
		xx=np.delete(xx,dice,0)
		ff=np.delete(ff,dice,0)
"""



# writing the output
###################################	

print ('xx shape: ',xx.shape) # checking all the sizes are ok
print ('yy shape: ',yy.shape)
print ('ff shape: ',ff.shape)


np.save(bname+'_X',xx)
np.save(bname+'_Y',yy)
np.save(bname+'_lasts',ff)


"""
# having a look
#################################
		
import matplotlib.pyplot as plt

for i in range(len(xx)):
	print(yy[i])

	# Plot the grid
	mx=np.matrix(xx[i])
	mf=np.matrix(ff[i])
	fig, axes = plt.subplots(nrows=1, ncols=2)
	axes.flat[0].imshow(mx)
	axes.flat[1].imshow(mf)
	plt.show()

"""



