import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import glob

# Choose the time step you want to visualize
time_step = 35  # Change this to visualize different steps

# Get all CSV files for the chosen time step
csv_files = glob.glob(f"velocity_step_{time_step}_rank_*.csv")

# Combine data from all ranks
df_list = [pd.read_csv(file) for file in csv_files]
df = pd.concat(df_list, ignore_index=True)

# Create a quiver plot for velocity vectors
plt.figure(figsize=(10, 8))
plt.quiver(df['x'], df['y'], df['u'], df['v'], scale=20)
plt.title(f"Velocity Field at Step {time_step}")
plt.xlabel("X")
plt.ylabel("Y")

# Save and show plot
plt.savefig(f"velocity_plot_step_{time_step}.png")
matplotlib.use("TkAgg")
plt.show()
