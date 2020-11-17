import matplotlib.pyplot as plt
import math

with open("test.txt") as f:
        line = f.readline()
        while (True):
                if (line == ""):
                        break
                if (line.startswith("Parameters")):
                        _, matrix_size, rank, percent_shown, cond_num = line.split(" ")
                        line = f.readline()
                        data = {}
                        while (line != "" and not line.startswith("Parameters")):
                                algo_name, iteration, time, fro_dist = line.split(" ")
                                iteration = int(iteration)
                                time = float(time)
                                fro_dist = float(fro_dist)
                                if algo_name not in data:
                                        data[algo_name] = []
                                data[algo_name].append((iteration, time, fro_dist))
                                line = f.readline()
                        
                        for algo_name in data.keys():
                                last = data[algo_name][-1]
                                data[algo_name].append((last[0] + 1, 60, last[2]))

                        for_legend = []
                        for algo_name, iterations in data.items(): 
                                for_legend.append(plt.plot([i[1] for i in iterations], [math.log10(i[2]) for i in iterations], label=algo_name)[0])
                        plt.legend(handles=for_legend)
                        plt.ylabel("Frobenius distance (log 10)")
                        plt.xlabel("Time in seconds")
                        plt.title(f"Size = {matrix_size}, Rank = {rank}, Shown = {percent_shown.strip()}%, Condition Number = {cond_num}")
                        plt.savefig(f"Matrix_{matrix_size}_{rank}_{percent_shown.strip()}_{cond_num.strip()}.png", bbox_inches='tight')
                        plt.close()
                else:
                        line = f.readline()