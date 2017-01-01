#!/usr/bin/python

import os
import subprocess
import math
import matplotlib.pyplot as plt

BINARY="/home/mic/workspace/przejsciowka/rds/cmake-build-release/rds"
MESH_FOLDER="/home/mic/workspace/przejsciowka/meshes"
OUTPUT_FOLDER="/home/mic/workspace/przejsciowka/temp"

d = {}

for method in range(6):
    for file in os.listdir(MESH_FOLDER):
        if file.endswith(".msh2"):
            numb = os.path.splitext(file)[0].split('_')[-1]
            a = int(numb)
            if a >= 500:
                print(numb)
                b = float(str(subprocess.run([BINARY, "-i", os.path.join(MESH_FOLDER, file), "-o", os.path.join(OUTPUT_FOLDER, "plot" + numb + ".dat"), "-m", str(method)], stdout=subprocess.PIPE).stdout).split("'")[1].split("\\")[0])
                d[a] = b
                print("done")

    x = []
    y = []
    for k in sorted(d.keys()):
        x.append(math.sqrt(2 / k))
        y.append(d[k])
    print(x)
    print(y)
'''
plt.plot(x, y)
plt.xscale('log')
plt.yscale('log')
plt.show()
'''
