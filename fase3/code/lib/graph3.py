import matplotlib.pyplot as plt

# Updated sequential time
best_sequential_time = 21.515  # seconds

# Updated data with new values added
thread_counts = [1, 2, 4, 8, 10, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80]
real_times = [2.107, 2.053, 1.72, 1.548, 1.564, 1.783, 2.114, 2.641, 2.794, 26.665, 31.325, 35.71, 42.85, 47.277, 57.683]

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
