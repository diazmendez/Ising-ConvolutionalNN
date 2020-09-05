""" ***

This is my first Neural Network code, it diferentiates ones and zeroes in handwriten images
coded in gray-scale matrices from MNIST.

Rogelio Diaz-Mendez

Stockholm, Aug 2, 2017
 
"""

""" ***
Updated to work with critical percolation images from Ising2D.

"""

""" ***
Last modification: 2017/09/20
"""



import numpy as np

from keras.models import Sequential
from keras.optimizers import Adadelta
from keras.layers import (Activation, Dropout, Flatten, Dense, Conv2D, Convolution2D, MaxPooling2D)
from keras import backend as K

import sys


nb_epoch=int(float(sys.argv[1]))
filex=str(sys.argv[2])
filey=str(sys.argv[3])
filef=str(sys.argv[4])
directory=str(sys.argv[5])


K.set_image_dim_ordering('th')




# load the data (already prepared by a previous code)
###########################################################

xx=np.load(filex) 
yy=np.load(filey)
ff=np.load(filef)

print('xx shape:', xx.shape) # print the sizes
print('yy shape:', yy.shape)






## Some model and data processing constants
############################################################################3

nb_classes = yy.shape[1]
img_rows, img_cols = xx.shape[2], xx.shape[3]	# input image dimensions


nb_filters = 2	# number of convolutional filters to use
nb_pool = 2	# size of pooling area for max pooling
nb_conv = 3	# convolution kernel size


model = Sequential()

model.add(Conv2D(nb_filters, (nb_conv, nb_conv), input_shape=(1, img_rows, img_cols), padding="valid"))
model.add(Activation('relu'))

model.add(Conv2D(nb_filters, (nb_conv, nb_conv)))
model.add(Activation('relu'))

model.add(MaxPooling2D(pool_size=(nb_pool, nb_pool)))
model.add(Dropout(0.25))
model.add(Flatten())

model.add(Dense(150))
model.add(Activation('relu'))
model.add(Dropout(0.5))

model.add(Dense(nb_classes))
model.add(Activation('softmax'))

#model.compile(loss='categorical_crossentropy', optimizer='adadelta', metrics=['accuracy'])
model.compile(loss='categorical_crossentropy', optimizer=Adadelta(lr=1,epsilon=1e-8), metrics=['accuracy'])










# learn what it is there to be learned
##################################################################3
batch_size = 100000
#nb_epoch = 100
v_split = 0.33

# fit the model
training = model.fit(xx, yy, batch_size=batch_size, epochs=nb_epoch, verbose=1, validation_split=v_split)






##### Saving the result

import os

# making the path
os.system("touch "+directory)
os.system("rm -r "+directory)
os.system("mkdir "+directory)
os.system("mkdir "+directory+"/source")
# copying the sources
os.system("cp "+sys.argv[0]+" "+directory+"/source/")
os.system("cp "+filex+" "+directory+"/source/")
os.system("cp "+filey+" "+directory+"/source/")
os.system("cp "+filef+" "+directory+"/source/")

# saving the training history
fop = open(directory+'/history.dat', 'a')	
fop.write("#epoch\tloss....\tacc.....\tval_loss\tval_acc\n")
for i in range (0,len(training.history['loss'])):
	line = str(i)+"\t"+str("%.6f"%training.history['loss'][i])+"\t"+ \
	str("%.6f"%training.history['acc'][i])+"\t"+ \
	str("%.6f"%training.history['val_loss'][i])+"\t"+ \
	str("%.6f"%training.history['val_acc'][i])+"\n"
	fop.write(line)
fop.close()	



#for i in range (0,len(m2bh['loss'])):
#	print (i,"\t","%.6f"%m2bh['loss'][i],"\t","%.6f"%m2bh['acc'][i],"\t","%.6f"%m2bh['val_loss'][i],"\t","%.6f"%m2bh['val_acc'][i])


# Saving the model

import json

saved_model = model.to_json()
with open(directory+"/saved_model_architecture.json", 'w') as outfile:
        json.dump(saved_model, outfile)

model.save_weights(directory+"/saved_model_weights.h5")



# saving the model summary (not a very elegant trick)
orig_stdout = sys.stdout
f = open(directory+"/model_summary.txt", 'w')
sys.stdout = f
print(model.summary())
sys.stdout = orig_stdout
f.close()




# script for model uploader
"""
from keras.models import model_from_json
import json

import sys


filename_architecture=str(sys.argv[1])
filename_weights=str(sys.argv[2])

   
### Load architecture
with open(filename_architecture, 'r') as architecture_file:
        model_architecture = json.load(architecture_file)
     
     
loaded_model = model_from_json(model_architecture)
     
### Load weights
loaded_model.load_weights(filename_weights)
"""





# having a look if
#################################
                

if (0):

	import matplotlib.pyplot as plt
	# Create the plot
	#plt.plot(model.history['val_loss'], 'r', model_2_training.history['val_loss'], 'b')
	plt.plot(training.history['val_loss'])
	plt.xlabel('Epochs')
	plt.ylabel('Validation score')
	plt.show()





if (0):

	pred=model.predict(xx)

	import matplotlib.pyplot as plt
	for i in range(len(xx)):
        	print('target: ',yy[i])
	        print('prediction: ',pred[i])

        	# Plot the grid
        	mx=np.matrix(xx[i])
        	mf=np.matrix(ff[i])
        	fig, axes = plt.subplots(nrows=1, ncols=2)
        	axes.flat[0].imshow(mx)
        	axes.flat[1].imshow(mf, vmin=-1,vmax=1)
        	plt.show()













