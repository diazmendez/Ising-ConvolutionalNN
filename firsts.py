"""
Reads the ising2d.c output files and writes them in the CNN input form
"""

"""
Just to prepare firsts.dat as an input  X: _X0
"""

import numpy as np
import sys

l=int(sys.argv[1])
path=sys.argv[2]
bname=sys.argv[3]


# upload and formating
#################################################

xx=np.loadtxt(path+'/firsts.dat'); # inputs
xx=np.reshape(xx,(-1, 1, l, l))

yy=np.loadtxt(path+'/targets.dat'); # targets


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



# writing the output
###################################	

print ('x0 shape: ',xx.shape) # checking all the sizes are ok
print ('y0 shape: ',yy.shape)


np.save(bname+'_X0',xx)
np.save(bname+'_Y0',yy)




