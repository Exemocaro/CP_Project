import matplotlib.pyplot as plt

# Assumed parallelization percentage
P = 0.989

# Number of threads (up to 20)
thread_counts = list(range(1, 21))

# Calculating maximum speedup using Amdahl's Law
amdahls_law_speedup = [1 / ((1 - P) + P / n) for n in thread_counts]

# Preparing the figure and axis
fig, ax = plt.subplots(figsize=(10, 5))

# Calculating the non-parallelizable part and parallelizable part
non_parallel_part = [(1 - P) * 100] * len(thread_counts)  # Non-parallelizable part does not change
parallel_part = [P * 100 / s for s in amdahls_law_speedup]  # Parallel part is reduced by speedup

# Stacking the non-parallelizable and parallelizable parts
ax.bar(thread_counts, non_parallel_part, label='Non-parallelizable Work')
ax.bar(thread_counts, parallel_part, bottom=non_parallel_part, label='Parallelizable Work')

# Adding the assumed parallelization percentage text at the center right
text_x = max(thread_counts) * 0.95  # Placing text near the right edge
text_y = max(parallel_part) / 2  # Centering text vertically
ax.text(text_x, text_y, f'f = {P}', fontsize=12, ha='center', va='center', color='black')

# Setting the labels and title
ax.set_xlabel('Number of Threads', fontsize=18)
ax.set_ylabel('Percentage of Work (%)', fontsize=18)
ax.set_title('Amdahl\'s Law - Work Distribution with Speedup', fontsize=18)
ax.set_xticks(thread_counts)  # Set the ticks to be at the thread counts
ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticklabels(thread_counts)

# Adding a legend
ax.legend()

# Displaying the graph
plt.tight_layout()
plt.savefig('../images/amdahl.png')
