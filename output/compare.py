data = {}
with open("lena.txt") as f:
        while (True):
                line = f.readline()
                if (line == ""):
                        break
                algo_name, iteration, time, fro_dist = line.split(" ")
                iteration = int(iteration)
                time = float(time)
                fro_dist = float(fro_dist)
                if algo_name not in data:
                        data[algo_name] = []
                data[algo_name].append((iteration, time, fro_dist))
        
import matplotlib.pyplot as plt
import random
import math
for_legend = []
for algo_name, iterations in data.items(): 
        r = random.random()
        b = random.random()
        g = random.random()
        color = (r, g, b)
        for_legend.append(plt.plot([i[1] for i in iterations], [math.log(i[2]) for i in iterations], label=algo_name, c=color)[0])
plt.legend(handles=for_legend)
plt.ylabel("Frobenius distance (log)")
plt.xlabel("Time in seconds")
plt.title("Log Frobenius Error by Time on Image of Lena")
plt.show()