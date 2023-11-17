from matplotlib import pyplot as plt
import numpy as np
import os

L=10
PlotData=np.array([])
# complete...

Directory = "./L={L}"
FileExtension = '.txt'

Files = os.listdir(Directory)
for file in Files:
    # Check if File has the desired extension
    if file.endswith(FileExtension):
        # Extract number from file title (removing "beta=" and ".txt")
        Beta = file.replace("beta=", "").replace(".txt", "")
        try:
            # Convert Beta to float number
            Beta = float(Beta)
            Energy2,Energy4=np.loadtxt(file,unpack=True)
            Binder=np.mean( Energy4 )/( (np.mean(Energy2))**2 )
            PlotData.append([Beta, Binder])
            
        except ValueError:
            print(f"The value '{Beta}' in '{file}' can't be converted.")
