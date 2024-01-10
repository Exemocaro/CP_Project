import matplotlib.pyplot as plt

# Data for the graph
processes = [1, 2, 4, 8, 20, 40]
instructions = [149, 164, 174, 181, 201, 227]
cycles = [89, 134, 162, 172, 215, 356]

# Plotting the graph
plt.figure(figsize=(10, 6))
plt.plot(processes, instructions, marker='o', label='Instructions (Billions)', color='blue')
plt.plot(processes, cycles, marker='o', label='Cycles (Billions)', color='orange')

# Adding labels and title with larger font sizes
plt.title('Instructions and Cycles Analysis', fontsize=36)
plt.xlabel('Processes', fontsize=36)
plt.ylabel('Billions', fontsize=36)
plt.tick_params(axis='both', which='major', labelsize=20)

# Adding legend with larger font size
plt.legend()

# Display the graph
plt.grid(True)
plt.tight_layout()
plt.savefig('../images/instructions_cycles_graph.png')
#plt.show()
