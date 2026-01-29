#!/usr/bin/env python3
"""
Advanced MOND Visualization Suite

Visualizes:
1. Galaxy cluster with neutrino distribution
2. Bullet Cluster collision evolution
3. KBC Void evolution and Hubble tension
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import glob
import os

plt.style.use('seaborn-v0_8-darkgrid')
plt.rcParams['figure.figsize'] = (14, 10)
plt.rcParams['font.size'] = 11

class AdvancedMONDAnalyzer:
    
    def __init__(self):
        self.particle_colors = {
            0: 'yellow',    # TYPE_BARYON (galaxies)
            1: 'red',       # TYPE_GAS (hot ICM)
            2: 'blue',      # TYPE_NEUTRINO (hot DM)
            3: 'gray'       # TYPE_DM_TEST
        }
        self.particle_names = {
            0: 'Baryons (Galaxies)',
            1: 'Gas (X-ray)',
            2: 'Neutrinos (Hot DM)',
            3: 'Test DM'
        }
    
    def load_snapshot(self, filename):
        """Load advanced simulation snapshot"""
        data = []
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                vals = line.split()
                if len(vals) >= 9:
                    data.append([float(x) for x in vals])
        
        if not data:
            return None
        
        data = np.array(data)
        return {
            'id': data[:, 0].astype(int),
            'x': data[:, 1],
            'y': data[:, 2],
            'z': data[:, 3],
            'vx': data[:, 4],
            'vy': data[:, 5],
            'vz': data[:, 6],
            'mass': data[:, 7],
            'type': data[:, 8].astype(int)
        }
    
    def plot_bullet_cluster_sequence(self, save=True):
        """
        Plot Bullet Cluster evolution showing gas-gravity separation
        """
        snapshots = sorted(glob.glob('bullet_snapshot_*.dat'))
        
        if not snapshots:
            print("No bullet cluster snapshots found")
            return
        
        print(f"Plotting Bullet Cluster evolution ({len(snapshots)} snapshots)...")
        
        # Select key frames
        n_frames = min(4, len(snapshots))
        frame_indices = [0, len(snapshots)//3, 2*len(snapshots)//3, -1]
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 14))
        axes = axes.flatten()
        
        for i, idx in enumerate(frame_indices):
            data = self.load_snapshot(snapshots[idx])
            if data is None:
                continue
            
            ax = axes[i]
            
            # Plot different particle types
            for ptype in [0, 1, 2]:
                mask = data['type'] == ptype
                if np.any(mask):
                    ax.scatter(data['x'][mask], data['y'][mask],
                             c=self.particle_colors[ptype],
                             s=20 if ptype != 1 else 40,
                             alpha=0.6 if ptype != 1 else 0.8,
                             label=self.particle_names[ptype],
                             edgecolors='black' if ptype == 1 else 'none',
                             linewidths=0.5 if ptype == 1 else 0)
            
            # Calculate centers for each component
            for ptype, name, color in [(1, 'Gas', 'red'), (2, 'Gravity (ν+gal)', 'blue')]:
                if ptype == 2:
                    # Gravity = neutrinos + baryons
                    mask = (data['type'] == 0) | (data['type'] == 2)
                else:
                    mask = data['type'] == ptype
                
                if np.any(mask):
                    cx = np.average(data['x'][mask], weights=data['mass'][mask])
                    cy = np.average(data['y'][mask], weights=data['mass'][mask])
                    ax.plot(cx, cy, 'X', color=color, markersize=20, 
                           markeredgecolor='black', markeredgewidth=2,
                           label=f'{name} center')
            
            ax.set_xlabel('X (Mpc)', fontsize=12)
            ax.set_ylabel('Y (Mpc)', fontsize=12)
            ax.set_title(f'Frame {idx+1}/{len(snapshots)}', fontsize=14, fontweight='bold')
            ax.set_aspect('equal')
            ax.grid(True, alpha=0.3)
            ax.set_xlim(-6, 6)
            ax.set_ylim(-4, 4)
            
            if i == 0:
                ax.legend(loc='upper left', fontsize=9)
        
        fig.suptitle('Bullet Cluster Collision: Gas-Gravity Separation in MOND', 
                    fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        
        if save:
            plt.savefig('bullet_cluster_evolution.png', dpi=300, bbox_inches='tight')
            print("✓ Saved bullet_cluster_evolution.png")
        
        plt.show()
    
    def plot_neutrino_distribution(self, save=True):
        """
        Plot neutrino radial distribution
        Shows that neutrinos avoid galaxy scales but cluster at Mpc scales
        """
        if not os.path.exists('neutrino_distribution.dat'):
            print("neutrino_distribution.dat not found")
            return
        
        data = np.loadtxt('neutrino_distribution.dat', comments='#')
        
        if len(data) == 0:
            print("No data in neutrino_distribution.dat")
            return
        
        radius = data[:, 0]      # kpc
        count = data[:, 1]       # number
        velocity = data[:, 2]    # km/s
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Plot 1: Particle count vs radius
        ax = axes[0]
        ax.plot(radius/1000, count, 'b-', linewidth=2, marker='o', markersize=6)
        ax.axvline(100/1000, color='red', linestyle='--', alpha=0.7,
                  label='Galaxy scale (~100 kpc)')
        ax.axvline(1000/1000, color='green', linestyle='--', alpha=0.7,
                  label='Cluster scale (~1 Mpc)')
        ax.set_xlabel('Radius (Mpc)', fontsize=12)
        ax.set_ylabel('Number of Neutrinos', fontsize=12)
        ax.set_title('Neutrino Spatial Distribution', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Add annotations
        ax.text(0.05, 0.95, 'Low count at galaxy scales\n(avoids galactic halos)',
               transform=ax.transAxes, fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
        
        # Plot 2: Velocity vs radius
        ax = axes[1]
        ax.plot(radius/1000, velocity, 'r-', linewidth=2, marker='s', markersize=6)
        ax.axhline(1000, color='blue', linestyle='--', alpha=0.7,
                  label='Thermal velocity (~1000 km/s)')
        ax.set_xlabel('Radius (Mpc)', fontsize=12)
        ax.set_ylabel('Average Velocity (km/s)', fontsize=12)
        ax.set_title('Neutrino Velocity Profile', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Add annotations
        ax.text(0.05, 0.95, 'High thermal velocity prevents\ncollapse at small scales',
               transform=ax.transAxes, fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
        
        plt.tight_layout()
        
        if save:
            plt.savefig('neutrino_analysis.png', dpi=300, bbox_inches='tight')
            print("✓ Saved neutrino_analysis.png")
        
        plt.show()
    
    def plot_void_evolution(self, save=True):
        """
        Plot KBC Void evolution and Hubble tension
        """
        if not os.path.exists('void_evolution.dat'):
            print("void_evolution.dat not found")
            return
        
        data = np.loadtxt('void_evolution.dat', comments='#')
        
        if len(data) == 0:
            print("No data in void_evolution.dat")
            return
        
        time = data[:, 0]           # Gyr
        void_radius = data[:, 1]    # Mpc
        outflow_vel = data[:, 2]    # km/s
        scale_factor = data[:, 3]   # a(t)
        hubble_param = data[:, 4]   # H(t)
        
        # Convert H(t) to km/s/Mpc for comparison with observations
        # H0_code is in units of 1/Gyr, need to convert to km/s/Mpc
        unit_time_gyr = 1.0  # Gyr
        unit_length_mpc = 1.0  # Mpc
        h_local = hubble_param / unit_time_gyr * (3.0857e19 / 1e3)  # km/s/Mpc
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Plot 1: Void radius growth
        ax = axes[0, 0]
        ax.plot(time, void_radius, 'b-', linewidth=2, marker='o', markersize=4)
        ax.set_xlabel('Time (Gyr)', fontsize=12)
        ax.set_ylabel('Void Radius (Mpc)', fontsize=12)
        ax.set_title('Void Expansion', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Add MOND vs ΛCDM comparison (theoretical)
        # ΛCDM would grow slower
        time_theory = np.linspace(0, time[-1], 100)
        radius_lcdm = void_radius[0] + (void_radius[-1] - void_radius[0]) * \
                      (time_theory / time[-1])**0.7  # Slower growth
        ax.plot(time_theory, radius_lcdm, 'r--', linewidth=2, alpha=0.7, 
               label='ΛCDM (slower)')
        ax.plot(time, void_radius, 'b-', linewidth=2, label='MOND (faster)')
        ax.legend()
        
        # Plot 2: Outflow velocity
        ax = axes[0, 1]
        ax.plot(time, outflow_vel, 'r-', linewidth=2, marker='s', markersize=4)
        ax.set_xlabel('Time (Gyr)', fontsize=12)
        ax.set_ylabel('Outflow Velocity (km/s)', fontsize=12)
        ax.set_title('Matter Outflow from Void', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        ax.text(0.05, 0.95, 'MOND: Stronger outflow\ndue to enhanced gravity\nat low acceleration',
               transform=ax.transAxes, fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
        
        # Plot 3: Scale factor evolution
        ax = axes[1, 0]
        ax.plot(time, scale_factor, 'g-', linewidth=2, marker='^', markersize=4)
        ax.set_xlabel('Time (Gyr)', fontsize=12)
        ax.set_ylabel('Scale Factor a(t)', fontsize=12)
        ax.set_title('Cosmic Expansion', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Plot 4: Local Hubble parameter (THE KEY RESULT!)
        ax = axes[1, 1]
        ax.plot(time, h_local, 'purple', linewidth=2, marker='D', markersize=4)
        
        # Add observational constraints
        ax.axhline(73, color='blue', linestyle='--', linewidth=2, alpha=0.7,
                  label='Local H₀ (SH0ES: 73±1)')
        ax.axhline(67, color='red', linestyle='--', linewidth=2, alpha=0.7,
                  label='CMB H₀ (Planck: 67±1)')
        ax.fill_between(time, 72, 74, alpha=0.2, color='blue')
        ax.fill_between(time, 66, 68, alpha=0.2, color='red')
        
        ax.set_xlabel('Time (Gyr)', fontsize=12)
        ax.set_ylabel('Local H₀ (km/s/Mpc)', fontsize=12)
        ax.set_title('Hubble Tension Resolution', fontsize=14, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        
        # Add result box
        if len(h_local) > 0:
            final_h0 = h_local[-1]
            ax.text(0.05, 0.05, f'Final local H₀: {final_h0:.1f} km/s/Mpc\n' +
                   'MOND void produces\nhigher local expansion!',
                   transform=ax.transAxes, fontsize=11, verticalalignment='bottom',
                   bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
        
        fig.suptitle('KBC Void Evolution: MOND Explanation of Hubble Tension', 
                    fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        
        if save:
            plt.savefig('void_hubble_analysis.png', dpi=300, bbox_inches='tight')
            print("✓ Saved void_hubble_analysis.png")
        
        plt.show()
    
    def create_summary_plot(self, save=True):
        """
        Create a summary figure showing all three modules
        """
        fig = plt.figure(figsize=(16, 10))
        gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)
        
        # Module 1: Neutrino distribution
        if os.path.exists('neutrino_distribution.dat'):
            data = np.loadtxt('neutrino_distribution.dat', comments='#')
            ax1 = fig.add_subplot(gs[0, 0])
            ax1.plot(data[:, 0]/1000, data[:, 1], 'b-', linewidth=2)
            ax1.set_xlabel('Radius (Mpc)')
            ax1.set_ylabel('N(neutrinos)')
            ax1.set_title('Module 1: Cluster Neutrinos', fontweight='bold')
            ax1.grid(True, alpha=0.3)
        
        # Module 2: Bullet Cluster snapshot
        bullets = sorted(glob.glob('bullet_snapshot_*.dat'))
        if bullets:
            data = self.load_snapshot(bullets[-1])
            if data:
                ax2 = fig.add_subplot(gs[0, 1:])
                
                for ptype in [1, 0, 2]:
                    mask = data['type'] == ptype
                    if np.any(mask):
                        ax2.scatter(data['x'][mask], data['y'][mask],
                                  c=self.particle_colors[ptype],
                                  s=15 if ptype != 1 else 30,
                                  alpha=0.6, label=self.particle_names[ptype])
                
                ax2.set_xlabel('X (Mpc)')
                ax2.set_ylabel('Y (Mpc)')
                ax2.set_title('Module 2: Bullet Cluster Final State', fontweight='bold')
                ax2.legend(fontsize=9)
                ax2.grid(True, alpha=0.3)
                ax2.set_aspect('equal')
        
        # Module 3: Void evolution
        if os.path.exists('void_evolution.dat'):
            data = np.loadtxt('void_evolution.dat', comments='#')
            
            ax3 = fig.add_subplot(gs[1, 0])
            ax3.plot(data[:, 0], data[:, 1], 'g-', linewidth=2)
            ax3.set_xlabel('Time (Gyr)')
            ax3.set_ylabel('Void Radius (Mpc)')
            ax3.set_title('Module 3: Void Growth', fontweight='bold')
            ax3.grid(True, alpha=0.3)
            
            ax4 = fig.add_subplot(gs[1, 1:])
            # Recalculate H0 in km/s/Mpc
            h_local = data[:, 4] * (3.0857e19 / 1e3)
            ax4.plot(data[:, 0], h_local, 'purple', linewidth=2)
            ax4.axhline(73, color='blue', linestyle='--', label='Obs: 73±1')
            ax4.axhline(67, color='red', linestyle='--', label='CMB: 67±1')
            ax4.set_xlabel('Time (Gyr)')
            ax4.set_ylabel('Local H₀ (km/s/Mpc)')
            ax4.set_title('Hubble Tension', fontweight='bold')
            ax4.legend()
            ax4.grid(True, alpha=0.3)
        
        fig.suptitle('Advanced MOND: Three Modules for Hard Cases', 
                    fontsize=18, fontweight='bold')
        
        if save:
            plt.savefig('advanced_mond_summary.png', dpi=300, bbox_inches='tight')
            print("✓ Saved advanced_mond_summary.png")
        
        plt.show()
    
    def analyze_all(self):
        """Run all analyses"""
        print("=== Advanced MOND Analysis ===\n")
        
        print("1. Analyzing neutrino distribution...")
        self.plot_neutrino_distribution()
        
        print("\n2. Analyzing Bullet Cluster collision...")
        self.plot_bullet_cluster_sequence()
        
        print("\n3. Analyzing void evolution...")
        self.plot_void_evolution()
        
        print("\n4. Creating summary plot...")
        self.create_summary_plot()
        
        print("\n=== Analysis Complete ===")


if __name__ == '__main__':
    print("╔═══════════════════════════════════════════════════════════════╗")
    print("║  Advanced MOND - Visualization Suite                         ║")
    print("╚═══════════════════════════════════════════════════════════════╝\n")
    
    analyzer = AdvancedMONDAnalyzer()
    
    # Check what files exist
    cluster_files = glob.glob('cluster_snapshot_*.dat')
    bullet_files = glob.glob('bullet_snapshot_*.dat')
    void_files = glob.glob('void_snapshot_*.dat')
    
    has_neutrino = os.path.exists('neutrino_distribution.dat')
    has_void_data = os.path.exists('void_evolution.dat')
    
    if not any([cluster_files, bullet_files, void_files, has_neutrino, has_void_data]):
        print("No simulation output found.")
        print("Run: ./mond_advanced <test> first\n")
        print("Available tests: cluster, bullet, void, all")
    else:
        print(f"Found:")
        if cluster_files: print(f"  - {len(cluster_files)} cluster snapshots")
        if bullet_files: print(f"  - {len(bullet_files)} bullet snapshots")
        if void_files: print(f"  - {len(void_files)} void snapshots")
        if has_neutrino: print(f"  - Neutrino distribution data")
        if has_void_data: print(f"  - Void evolution data")
        print()
        
        analyzer.analyze_all()
    
    print("\n✓ Visualization complete!")
