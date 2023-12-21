import matplotlib.pyplot as plt

# Updated sequential time
best_sequential_time = 21.515  # seconds

# Updated data with new values added
thread_counts = [1, 2, 4, 8, 10, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80]
real_times = [22.115, 11.175, 5.665, 2.888, 2.376, 1.982, 1.533, 1.265, 1.225, 1.187, 1.144, 1.118, 1.089, 1.278, 1.315, 1.339, 1.384, 1.421, 1.471, 1.506, 1.555, 1.576, 1.617]

# Calculating speedup
speedup = [best_sequential_time / time for time in real_times]

# Ideal speedup calculation: Linear up to 20 threads, then constant
ideal_speedup = [min(tc, 20) for tc in thread_counts]
ideal_speedup = [ideal_speedup[i] if thread_counts[i] <= 20 else 20 for i in range(len(thread_counts))]

# Plotting the graph for real times
plt.figure(figsize=(10, 8))
plt.plot(thread_counts, real_times, marker='o', color='blue')
plt.title('Execution Time vs Number of Threads', fontsize=30)
plt.xlabel('Number of Threads', fontsize=36)
plt.ylabel('Real Time (seconds)', fontsize=36)
plt.xticks(range(0, max(thread_counts) + 1, 4))
plt.tick_params(axis='both', which='major', labelsize=20)
plt.grid(True)
plt.tight_layout()
plt.savefig('../images/time.png')

# Clear plot
plt.clf()

# Plotting the graph for speedup
plt.figure(figsize=(10, 8))
plt.plot(thread_counts, speedup, marker='o', color='red')
plt.plot(thread_counts, ideal_speedup, linestyle='--', color='green', label='Ideal Speedup')  # Corrected ideal speedup line
plt.title('Speedup vs Number of Threads', fontsize=30)
plt.xlabel('Number of Threads', fontsize=36)
plt.ylabel('Speedup', fontsize=36)
plt.xticks(range(0, max(thread_counts) + 1, 4))
plt.tick_params(axis='both', which='major', labelsize=20)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('../images/speedup.png')
