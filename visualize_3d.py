import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import glob

# Load data
snapshot_files = sorted(glob.glob("snapshot_*.dat"))
if not snapshot_files:
    print("❌ No snapshots found!")
    exit()

fig = plt.figure(figsize=(10, 8), facecolor='black')
ax = fig.add_subplot(111, projection='3d')
ax.set_facecolor('black')

def update(frame_idx):
    ax.clear()
    ax.set_facecolor('black')
    # Hide axes for a cinematic look
    ax.grid(False)
    ax.set_axis_off()
    
    try:
        data = np.loadtxt(snapshot_files[frame_idx])
        # x, y, z are columns 0, 1, 2
        ax.scatter(data[:, 0], data[:, 1], data[:, 2], s=1, c='#00ffff', alpha=0.7)
        ax.set_title(f"MOND 3D Evolution: {frame_idx*100} Myr", color='white', fontsize=12)
        # Set consistent limits so it doesn't "jump"
        ax.set_xlim([-50, 50]); ax.set_ylim([-50, 50]); ax.set_zlim([-20, 20])
    except Exception as e:
        print(f"Error: {e}")

anim = FuncAnimation(fig, update, frames=len(snapshot_files), interval=100)

print("🎥 Rendering 3D Animation to galaxy_3d.gif (this may take a minute)...")
# We use 'pillow' because it doesn't require extra software like ffmpeg
anim.save('galaxy_3d.gif', writer='pillow', fps=10)
print("✅ Done! galaxy_3d.gif created.")
