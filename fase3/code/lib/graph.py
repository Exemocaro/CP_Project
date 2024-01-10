import matplotlib.pyplot as plt

# Updated sequential time
best_sequential_time = 21.515  # seconds

# Updated data with new values added
thread_counts = [1, 2, 4, 8, 10, 12, 16, 20, 24, 28, 32, 36, 40]
mpi_times = [28.436, 22.928, 13.440, 7.820, 6.496, 5.528, 4.322, 3.657, 3.149, 4.017, 3.917, 3.627, 3.785]
openmp_times = [22.115, 11.175, 5.665, 2.888, 2.376, 1.982, 1.533, 1.265, 1.225, 1.187, 1.144, 1.118, 1.089]
joint_approach_times = [2.107, 2.053, 1.720, 1.548, 1.564, 1.783, 2.114, 2.641, 2.794, 26.665, 31.325, 35.710, 42.850]

# Calculating speedup
mpi_speedup = [best_sequential_time / time for time in mpi_times]
openmp_speedup = [best_sequential_time / time for time in openmp_times]
joint_approach_speedup = [best_sequential_time / time for time in joint_approach_times]

# print speed up with values rounded for 3 cases
for i in range(len(thread_counts)):
    print("MPI Speedup for {} threads: {:.3f}".format(thread_counts[i], mpi_speedup[i]))
    print("OpenMP Speedup for {} threads: {:.3f}".format(thread_counts[i], openmp_speedup[i]))
    print("Joint Approach Speedup for {} threads: {:.3f}".format(thread_counts[i], joint_approach_speedup[i]))

# Ideal speedup calculation: Linear up to 20 threads, then constant
ideal_speedup = [min(tc, 20) for tc in thread_counts]
ideal_speedup = [ideal_speedup[i] if thread_counts[i] <= 20 else 20 for i in range(len(thread_counts))]

# Plotting the graph for real times
plt.figure(figsize=(10, 8))
plt.plot(thread_counts, mpi_times, marker='o', color='blue', label='MPI Time')
plt.plot(thread_counts, openmp_times, marker='o', color='orange', label='OpenMP Time')
plt.plot(thread_counts, joint_approach_times, marker='o', color='purple', label='Joint Approach Time')
plt.title('Execution Time', fontsize=30)
plt.xlabel('Number of Threads/Processes', fontsize=36)
plt.ylabel('Real Time (seconds)', fontsize=36)
plt.xticks(range(0, max(thread_counts) + 1, 4))
plt.yticks(range(0, int(max(mpi_times + openmp_times + joint_approach_times)) + 1, 5))  # Set more ticks on the vertical axis
plt.tick_params(axis='both', which='major', labelsize=20)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('../images/time.png')

# Clear plot
plt.clf()

# Plotting the graph for speedup
plt.figure(figsize=(10, 8))
plt.plot(thread_counts, mpi_speedup, marker='o', color='red', label='MPI Speedup')
plt.plot(thread_counts, openmp_speedup, marker='o', color='green', label='OpenMP Speedup')
plt.plot(thread_counts, joint_approach_speedup, marker='o', color='purple', label='Joint Approach Speedup')
plt.plot(thread_counts, ideal_speedup, linestyle='--', color='gray', label='Ideal Speedup')  # Corrected ideal speedup line
plt.title('Speedup', fontsize=30)
plt.xlabel('Number of Threads/Processes', fontsize=36)
plt.ylabel('Speedup', fontsize=36)
plt.xticks(range(0, max(thread_counts) + 1, 4))
plt.yticks(range(0, int(max(mpi_speedup + openmp_speedup + joint_approach_speedup)) + 1, 2))  # Set more ticks on the vertical axis
plt.tick_params(axis='both', which='major', labelsize=20)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('../images/speedup.png')
