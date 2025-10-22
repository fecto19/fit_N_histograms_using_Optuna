import numpy as np
import tensorflow as tf
from tensorflow import keras
import os

norm = np.load("normalization_constants.npz")
x_mean, x_std = norm["x_mean"], norm["x_std"]
y_mean, y_std = norm["y_mean"], norm["y_std"]


input_file = "output_SEED8.txt"
print(f"Loading input file: {input_file}")
data = np.loadtxt(input_file, skiprows=1).astype("float32")

X = data[:, [0, 3, 4]]
X_norm = (X - x_mean) / x_std
nevents, nsamples = X_norm.shape
X_norm = X_norm.reshape(nevents, nsamples, 1)


model_dir_prefix = "model_"
n_models = 120 

for i in range(n_models):
    model_name = f"{model_dir_prefix}{i}"

    print(f"Loading model {i}...")
    model = keras.models.load_model(model_name, compile=False)

    print(f"Predicting with model {i}...")
    Y_pred_norm = model.predict(X_norm, batch_size=512, verbose=0)
    Y_pred = Y_pred_norm * y_std + y_mean

    output_filename = f"predictions_model_{i}.txt"
    np.savetxt(output_filename, Y_pred, fmt="%.6e")
    print(f"Saved predictions to {output_filename}")

print("All predictions saved.")
