import matplotlib.pyplot as plt

# Data provided by the user
thread_counts = [1, 2, 4, 8, 10, 12, 16, 20, 24, 28, 32, 36, 40]
real_times = [23.082, 11.633, 5.864, 2.993, 2.433, 2.041, 1.595, 1.314, 1.251, 1.215, 1.142, 1.105, 1.089]
best_sequential_time = 24.000  # seconds

# Calculating speedup
speedup = [best_sequential_time / time for time in real_times]

# Plotting the graph for real times
plt.figure(figsize=(14, 6))

plt.subplot(1, 2, 1)
plt.plot(thread_counts, real_times, marker='o', color='blue')
plt.title('Execution Time vs Number of Threads')
plt.xlabel('Number of Threads')
plt.ylabel('Real Time (seconds)')
plt.grid(True)

# Plotting the graph for speedup
plt.subplot(1, 2, 2)
plt.plot(thread_counts, speedup, marker='o', color='red')
plt.title('Speedup vs Number of Threads')
plt.xlabel('Number of Threads')
plt.ylabel('Speedup')
plt.grid(True)

plt.tight_layout()
# save plot
plt.savefig('graph.png')
#plt.show()