import numpy as np
import tensorflow as tf
from tensorflow import keras
from keras.layers import TFSMLayer
from keras import Sequential
import sys


print("Imported: apply AI model script")



#model = keras.models.load_model("model_88", compile=False)
model = Sequential([TFSMLayer("/home/fecto19/Particle_Physics_FzF/Bachelor_Thesis/fit_N_histograms_with_Optuna/initial_files_needed/AI_scatt_params_to_corr_func/model_88", call_endpoint="serving_default")])


# Load training normalization constants

norms = np.load("/home/fecto19/Particle_Physics_FzF/Bachelor_Thesis/fit_N_histograms_with_Optuna/initial_files_needed/AI_scatt_params_to_corr_func/normalization_constants.npz")

#print("PRINTING NORMS.FILES:\n")
#print(norms.files)

#print("PRINting norms:\n", norms)

x_mean, x_std, y_mean, y_std = norms["x_mean"], norms["x_std"], norms["y_mean"], norms["y_std"]

#print("PRINTING NORMS:\n")
#print(norms.files)


#x_mean = np.array([X1_mean, X2_mean, X3_mean], dtype="float32")
#x_std  = np.array([X1_std,  X2_std,  X3_std], dtype="float32")
#y_mean = np.array([Y1_mean, Y2_mean], dtype="float32")
#y_std  = np.array([Y1_std,  Y2_std], dtype="float32")

def predict_values(x1, x2, x3):
    
    X = np.array([[x1, x2, x3]], dtype="float32")
    X = (X - x_mean) / x_std
    X = X.reshape(1, 3, 1)  

    y_pred_norm = model.predict(X, verbose=0)
    #print("y_pred_norm:\n", y_pred_norm)
    #print("type of y_pred_norm:\n", type(y_pred_norm))
    y_pred_norm_array = list(y_pred_norm.values())[0]

    # denormalize
    y_pred = y_pred_norm_array * y_std + y_mean
    return y_pred[0]

def pass_optuna_trial_guess(m, f, d):
    out1, out2 = predict_values(m, f, d)
    print("Predicted outputs:", out1, out2)
    return out1, out2

