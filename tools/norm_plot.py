#!/usr/bin/python

import os
import subprocess
import math
import matplotlib.pyplot as plt

BINARY="/home/clapik/workspace/rds/bin/Debug/rds"
MESH_FOLDER="/home/clapik/workspace/meshes"
OUTPUT_FOLDER="/home/clapik/workspace/temp"

d = {}

for file in os.listdir(MESH_FOLDER):
    if file.endswith(".msh2"):
        numb = os.path.splitext(file)[0].split('_')[-1]
        a = int(numb)
        if a >= 100:
            b = float(str(subprocess.run([BINARY, os.path.join(MESH_FOLDER, file), os.path.join(OUTPUT_FOLDER, "plot" + numb + ".dat")], stdout=subprocess.PIPE).stdout).split("'")[1].split("\\")[0])
            d[a] = b

x = []
y = []
for k in sorted(d.keys()):
    x.append(math.sqrt(1 / k))
    y.append(d[k])
print(x)
print(y)
'''
plt.plot(x, y)
plt.xscale('log')
plt.yscale('log')
plt.show()
'''
