import matplotlib.pyplot as plt
import math


all_data = {}

files = ["results/randoms.txt"]
for file_name in files:
	with open(file_name) as f:
		line = f.readline()
		while (True):
			if (line == ""):
					break
			if (line.startswith("Parameters")):
				split = line.split(" ")
				key = (int(split[1]), int(split[2]), int(split[3]), int(split[4]))
				if key not in all_data:
						all_data[key] = {}
				data = all_data[key]
				line = f.readline()
				while (line != "" and not line.startswith("Parameters")):
					algo_name, iteration, time, fro_dist = line.split(" ")
					if float(time) < 60:
						iteration = int(iteration)
						time = float(time)
						fro_dist = float(fro_dist)
						if algo_name not in data:
								data[algo_name] = []
						data[algo_name].append((iteration, time, fro_dist))
					line = f.readline()
			else:
					line = f.readline()

import math
result = {}
print(len(all_data))
for params, data in all_data.items():
	if len(data) == 0:
		continue
	best = float("inf")
	for iterations in data.values():
		for i in iterations:
			val = math.log10(i[2])
			if val < best:
				best = val
	best = math.floor(best)
	if best not in result:
		result[best] = 0
	result[best] += 1
s = list(result.items())
s.sort()
print(s)
lower = -10
upper = 0
lower_count = 0
upper_count = 0
total_count = 0
for val, count in s:
	total_count += count
	if val <= lower:
		lower_count += count
	if val >= upper:
		upper_count += count

middle_count = total_count - upper_count - lower_count
print(lower_count, middle_count, upper_count, lower_count / total_count, middle_count / total_count, upper_count / total_count)

for params, data in all_data.items():
	
	matrix_size, rank, percent_shown, cond_num = params
	for algo_name in data.keys():
		last = data[algo_name][-1]
		data[algo_name].append((last[0] + 1, 60, last[2]))
		first = data[algo_name][0]
		data[algo_name].insert(0, (first[0] - 1, 0, first[2]))
		
	for_legend = []
	for algo_name, iterations in data.items(): 
		for_legend.append(plt.plot([i[1] for i in iterations], [math.log10(i[2]) for i in iterations], label=algo_name)[0])
	plt.legend(handles=for_legend)
	plt.ylabel("Frobenius distance (log 10)")
	plt.xlabel("Time in seconds")
	plt.title(f"Size = {matrix_size}, Rank = {rank}, Shown = {percent_shown}%, Condition Number = {cond_num}")
	plt.savefig(f"Matrix_{matrix_size}_{rank}_{percent_shown}_{cond_num}.png", bbox_inches='tight')
	plt.close()
