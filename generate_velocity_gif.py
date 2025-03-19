import pandas as pd
import matplotlib.pyplot as plt
import glob
import imageio

# Define the number of time steps
num_steps = 100  # Make sure this matches your C++ simulation

# Create a list to store image filenames
image_filenames = []

# Generate velocity plots for each time step
for t in range(num_steps):
    # Get all CSV files for the current time step (across all MPI ranks)
    csv_files = glob.glob(f"velocity_step_{t}_rank_*.csv")
    
    if not csv_files:
        print(f"No data found for time step {t}, skipping...")
        continue

    # Combine data from all ranks
    df_list = [pd.read_csv(file) for file in csv_files]
    df = pd.concat(df_list, ignore_index=True)

    # Create a quiver plot for velocity vectors
    plt.figure(figsize=(10, 8))
    plt.quiver(df['x'], df['y'], df['u'], df['v'], scale=20)
    plt.title(f"Velocity Field at Step {t}")
    plt.xlabel("X")
    plt.ylabel("Y")

    # Save the image
    image_filename = f"velocity_plot_step_{t}.png"
    plt.savefig(image_filename, dpi=300)
    image_filenames.append(image_filename)

    # Close the figure to save memory
    plt.close()

# Create a GIF from the saved images
gif_filename = "velocity_simulation.gif"
with imageio.get_writer(gif_filename, mode="I", duration=0.2) as writer:
    for image_filename in image_filenames:
        writer.append_data(imageio.imread(image_filename))

print(f"✅ GIF created: {gif_filename}")
