import numpy as np
import matplotlib.pyplot as plt

rawpoints = np.loadtxt("sigg1.csv", delimiter=",")

high=15
low=0
mins=np.min(rawpoints,axis=0)
maxs = np.max(rawpoints,axis=0)
rng=maxs-mins
scaled_points=high-(((high-low)*(max-rawpoints))/rng)
with open('sigg1_scaled.csv','w',newline="")as f:
    writer = csv.writer(f)
    writer.writerows(scaled_points)

plt.contourf(low,high,scaled_points)
plt.show()

