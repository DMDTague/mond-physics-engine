import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import glob
import os

# 1. Setup the 3D Plot
fig = plt.figure(figsize=(10, 10), facecolor='black')
ax = fig.add_subplot(111, projection='3d', facecolor='black')

# Hide the grid/axes for that "Space" look
ax.set_axis_off()
ax.grid(False)

# 2. Load Data
snapshot_files = sorted(glob.glob('snapshot_*.dat'), key=lambda x: int(x.split('_')[1].split('.')[0]))

if not snapshot_files:
    print("❌ No data found! You must run './mond_sim' and type 'y' first.")
    exit()

print(f"🚀 Found {len(snapshot_files)} frames. Loading engine...")

# 3. Animation Update Function
scatter = ax.scatter([], [], [], s=2, c='cyan', alpha=0.6)

def update(frame_idx):
    # Clear previous frame data
    ax.clear()
    ax.set_axis_off()
    
    # Set camera angle (rotate slowly)
    ax.view_init(elev=30, azim=frame_idx * 2)
    
    # Set fixed limits so the universe doesn't jump around
    limit = 30 # kpc
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_zlim(-limit, limit)
    
    # Load data for this frame
    try:
        data = np.loadtxt(snapshot_files[frame_idx])
        # Assuming columns are: x, y, z, vx, vy, vz, mass, type
        x = data[:, 0]
        y = data[:, 1]
        z = data[:, 2]
        
        # Plot stars
        ax.scatter(x, y, z, s=3, c='#00ffff', marker='o', alpha=0.8)
        
        # Add Title
        ax.text2D(0.05, 0.95, f"MOND Physics Engine\nTime: {frame_idx*100} Myr", 
                  transform=ax.transAxes, color='white', fontsize=14)
        
    except Exception as e:
        print(f"Frame {frame_idx} error: {e}")

# 4. Run Animation
anim = FuncAnimation(fig, update, frames=len(snapshot_files), interval=50)

print("🎥 Rendering 3D Physics...")
plt.show()
