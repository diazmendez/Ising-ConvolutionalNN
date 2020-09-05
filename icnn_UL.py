

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




