import numpy as np
import keras
import sys

# load the data (already prepared by a previous code)
xx=np.load(filex)
yy=np.load(filey)
ff=np.load(filef)

## Some model and data processing constants
nb_classes = yy.shape[1]
img_rows, img_cols = xx.shape[2], xx.shape[3]   # input image dimensions



model=keras.models.Sequential()
model.add(keras.layers.Dense(50,activation='relu',input_shape=(1,img_rows,img_cols)) )
model.add(keras.layers.Dense(10,activation='relu'))
model.add(keras.layers.core.Flatten())
model.add(keras.layers.Dense(nb_classes,activation='softmax'))

model.compile(optimizer='adam',loss='categorical_crossentropy',metrics=['accuracy'])


# learn what it is there to be learned
training=model.fit(xx,yy,epochs=30,batch_size=1000,validation_split=0.33)




