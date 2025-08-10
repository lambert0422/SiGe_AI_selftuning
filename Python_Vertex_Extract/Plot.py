import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Read the data files
df = pd.read_csv('bandedges.dat', sep='\s+', comment='#', header=0)
df_density = pd.read_csv('density_hole.dat', sep='\s+', comment='#', header=0)

# Create figure and primary axis
fig, ax1 = plt.subplots(figsize=(12, 6))

# Plot band edges on primary y-axis (left)
if 'HH[eV]' in df.columns:
    ax1.plot(df.iloc[:, 0], df['HH[eV]'],
             label='Heavy Hole Band',
             color='red',
             linewidth=2)

    ax1.set_xlabel(df.columns[0], fontsize=12)
    ax1.set_ylabel('Energy (eV)', fontsize=12)
    ax1.set_ylim(-0.2, 0.05)  # Energy range
    ax1.set_xlim(-30, 15)  # X-axis range
    ax1.tick_params(axis='y', labelcolor='red')
    ax1.grid(True, alpha=0.3)
else:
    print(f"Error: 'HH[eV]' column not found. Available columns: {list(df.columns)}")

# Create secondary y-axis (right)
ax2 = ax1.twinx()

# Plot hole density on secondary y-axis
if len(df_density.columns) >= 2:  # Check if density file has at least 2 columns
    ax2.plot(df_density.iloc[:, 0], df_density.iloc[:, 1],
             label='Hole Density',
             color='blue',
             linewidth=2,
             linestyle='--')

    ax2.set_ylabel('Density (1e18 cm^-3)', fontsize=12)
    ax2.tick_params(axis='y', labelcolor='blue')
else:
    print("Error: Density file doesn't have enough columns")
ax2.set_ylim(0, 0.05)  # Energy range
# Add legends
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

plt.title('Band Edge and Hole Density', fontsize=14)
plt.tight_layout()
plt.savefig('band_and_density.png', dpi=300)
plt.show()