import numpy as np

#data=np.loadtxt('L32s100000/magnetizations.dat')
data=np.loadtxt('../bising2d/bL32s100000/magnetizations.dat')



# balancing the data
##################################

classes=3

unbalance=np.zeros(classes)
ndata=len(data)	
yy=np.zeros((ndata,classes))

for i in range(len(data)):

	perc=data[i][0]
	mp=data[i][1]
	mf=data[i][2]


	if (mf==1.0): 
		yy[i,1]=1
	else: 
		if (mf==-1.0): 
			yy[i,2]=1
		else: 
			yy[i,0]=1
	







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
		data=np.delete(data,dice,0)







for k in range(50):

	th=.3*k/50

	acc=0

	for i in range(len(data)):

		perc=data[i][0]
		mp=data[i][1]
		mf=data[i][2]

		fstat=0
		if (mf==1.0): fstat=1
		if (mf==-1.0): fstat=-1


		pron=0
		if (mp>th): pron=1
		if (mp<-th): pron=-1
	
	
		if (pron==fstat): acc=acc+1

	acc=acc/len(data)

	print(th, acc)



	








