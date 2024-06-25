import matplotlib.pyplot as plt
import numpy as np
import os

def save_plot(dist1, dist2, dist3, area, output_folder):
    # Squeeze data for plotting
    df1 = np.squeeze(dist1)
    df2 = np.squeeze(dist2)
    df3 = np.squeeze(dist3)
    df4 = np.squeeze(area)

    # Plot distances and area
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111)
    ax.set_title('Distances between the Center of the Mass of Nucleobases')
    ax.plot(df1, color='red', label="d1")
    ax.plot(df2, color='blue', label="d2")
    ax.plot(df3, color='green', label="d3")
    ax.plot(df4, color='black', label="area")

    ax.set_xlabel('Frame')
    ax.set_ylabel('Distance (Å)')
    ax.legend(loc='best')

    # Create the directory if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Save the plot using an absolute path
    plot_file = os.path.join(output_folder, 'distances_and_area_plot.png')
    plt.savefig(plot_file)
    plt.close(fig)  # Close the figure to free up resources

    print(f"Plot saved: {plot_file}")
