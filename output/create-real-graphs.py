files = ["movielens2.txt", "lena.txt"]
names = ["Movie Lens", "Lena"]
for i, file in enumerate(files):
        data = {}
        with open(file) as f:
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
        import math
        for_legend = []
        for algo_name, iterations in data.items(): 
                for_legend.append(plt.plot([i[1] for i in iterations], [math.log(i[2]) for i in iterations], label=algo_name)[0])
        plt.legend(handles=for_legend)
        plt.ylabel("Frobenius distance (log)")
        plt.xlabel("Time in seconds")
        plt.title(f"Log Frobenius Error by Time on {names[i]}")
        s = names[i].replace(" ", "_")
        plt.savefig(f"Matrix_{s}.png", bbox_inches='tight')
        plt.close()