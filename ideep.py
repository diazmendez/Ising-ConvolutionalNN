import numpy as np
import keras
import sys


nb_epoch=int(float(sys.argv[1]))
filex=str(sys.argv[2])
filey=str(sys.argv[3])
filef=str(sys.argv[4])
directory=str(sys.argv[5])



# load the data (already prepared by a previous code)
###########################################################

xx=np.load(filex)
yy=np.load(filey)
ff=np.load(filef)

#xx0=np.load(filex)
#yy0=np.load(filey)
#ff0=np.load(filef)
#xx=xx0[40000:49000]
#yy=yy0[40000:49000]
#ff=ff0[40000:49000]







print('xx shape:', xx.shape) # print the sizes
print('yy shape:', yy.shape)


## Some model and data processing constants
############################################################################3

nb_classes = yy.shape[1]
img_rows, img_cols = xx.shape[2], xx.shape[3]   # input image dimensions



model=keras.models.Sequential()
model.add(keras.layers.Dense(50,activation='relu',input_shape=(1,img_rows,img_cols)) )
#model.add(keras.layers.Dense(100,activation='relu'))
model.add(keras.layers.Dense(10,activation='relu'))
#model.add(keras.layers.Dense(50,activation='relu'))
model.add(keras.layers.core.Flatten())
model.add(keras.layers.Dense(nb_classes,activation='softmax'))

model.compile(optimizer='adam',loss='categorical_crossentropy',metrics=['accuracy'])





# learn what it is there to be learned
##################################################################3
batch_size=10000

training=model.fit(xx,yy,epochs=nb_epoch,batch_size=batch_size,validation_split=0.33)





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



