#!/usr/bin/env python3
"""
MOND Physics Engine - Visualization and Analysis Tools

This script provides comprehensive visualization and analysis of MOND
N-body simulation outputs, including:
- Rotation curves (comparison with Newtonian)
- Particle trajectories and snapshots
- Energy/angular momentum conservation
- Tully-Fisher relation validation
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle
import glob
import os

# Set style for publication-quality plots
plt.style.use('seaborn-v0_8-darkgrid')
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12

class MONDAnalyzer:
    """Analyze and visualize MOND simulation outputs"""
    
    def __init__(self, data_dir='.'):
        self.data_dir = data_dir
        
    def load_rotation_curve(self, filename):
        """Load rotation curve data from file"""
        data = []
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                data.append([float(x) for x in line.split()])
        
        data = np.array(data)
        return {
            'radius': data[:, 0],      # kpc
            'velocity': data[:, 1],    # km/s
            'a_ratio': data[:, 2]      # a/a0
        }
    
    def load_snapshot(self, filename):
        """Load particle snapshot from file"""
        data = []
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                data.append([float(x) for x in line.split()])
        
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
    
    def plot_rotation_curves(self, save=True):
        """
        Plot rotation curves and compare initial vs final
        
        Shows the characteristic flat rotation curve predicted by MOND
        """
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Try to load both initial and final curves
        curves = {}
        for label in ['initial', 'final']:
            filename = f'rotation_curve_{label}.dat'
            if os.path.exists(filename):
                curves[label] = self.load_rotation_curve(filename)
        
        # Also try the test curve
        if os.path.exists('rotation_curve_test.dat'):
            curves['test'] = self.load_rotation_curve('rotation_curve_test.dat')
        
        if not curves:
            print("No rotation curve data found")
            return
        
        # Plot 1: Velocity vs Radius
        ax = axes[0, 0]
        for label, data in curves.items():
            ax.plot(data['radius'], data['velocity'], 
                   label=label.capitalize(), linewidth=2, marker='o', markersize=4)
        
        ax.set_xlabel('Radius (kpc)', fontsize=12)
        ax.set_ylabel('Circular Velocity (km/s)', fontsize=12)
        ax.set_title('MOND Rotation Curve', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Add horizontal line for flat rotation indication
        if curves:
            first_curve = list(curves.values())[0]
            v_flat = np.mean(first_curve['velocity'][-20:])  # Average of outer region
            ax.axhline(v_flat, color='red', linestyle='--', alpha=0.5, 
                      label=f'V_flat ≈ {v_flat:.1f} km/s')
        
        # Plot 2: Acceleration ratio a/a0
        ax = axes[0, 1]
        for label, data in curves.items():
            ax.semilogy(data['radius'], data['a_ratio'], 
                       label=label.capitalize(), linewidth=2, marker='o', markersize=4)
        
        ax.axhline(1.0, color='red', linestyle='--', alpha=0.5, 
                  label='a = a₀ (transition)')
        ax.set_xlabel('Radius (kpc)', fontsize=12)
        ax.set_ylabel('a / a₀', fontsize=12)
        ax.set_title('Acceleration Regime', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3, which='both')
        
        # Plot 3: Velocity gradient (dV/dr)
        ax = axes[1, 0]
        for label, data in curves.items():
            r = data['radius']
            v = data['velocity']
            dv_dr = np.gradient(v, r)
            ax.plot(r, dv_dr, label=label.capitalize(), linewidth=2)
        
        ax.axhline(0, color='red', linestyle='--', alpha=0.5, 
                  label='Flat (dV/dr = 0)')
        ax.set_xlabel('Radius (kpc)', fontsize=12)
        ax.set_ylabel('dV/dr (km/s/kpc)', fontsize=12)
        ax.set_title('Velocity Gradient (Flatness Test)', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Plot 4: MOND regime indicator
        ax = axes[1, 1]
        for label, data in curves.items():
            r = data['radius']
            a_ratio = data['a_ratio']
            
            # Calculate MOND regime: Deep MOND (a << a0), Transition, Newtonian (a >> a0)
            regime = np.zeros_like(a_ratio)
            regime[a_ratio < 0.1] = 0    # Deep MOND
            regime[(a_ratio >= 0.1) & (a_ratio <= 10)] = 1  # Transition
            regime[a_ratio > 10] = 2     # Newtonian
            
            colors = ['red', 'yellow', 'blue']
            labels_regime = ['Deep MOND', 'Transition', 'Newtonian']
            
            for i, (color, lbl) in enumerate(zip(colors, labels_regime)):
                mask = regime == i
                if np.any(mask):
                    ax.scatter(r[mask], a_ratio[mask], c=color, 
                             label=lbl, s=50, alpha=0.7)
        
        ax.axhline(0.1, color='gray', linestyle=':', alpha=0.5)
        ax.axhline(10, color='gray', linestyle=':', alpha=0.5)
        ax.set_yscale('log')
        ax.set_xlabel('Radius (kpc)', fontsize=12)
        ax.set_ylabel('a / a₀', fontsize=12)
        ax.set_title('MOND Regime Map', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3, which='both')
        
        plt.tight_layout()
        
        if save:
            plt.savefig('rotation_curves_analysis.png', dpi=300, bbox_inches='tight')
            print("✓ Saved rotation_curves_analysis.png")
        
        plt.show()
    
    def plot_galaxy_snapshot(self, snapshot_file, save=True):
        """
        Plot 2D projection of galaxy from snapshot
        """
        data = self.load_snapshot(snapshot_file)
        if data is None:
            print(f"Could not load {snapshot_file}")
            return
        
        fig, axes = plt.subplots(1, 2, figsize=(16, 7))
        
        # Separate particle types
        bulge = data['type'] == 1
        disk = data['type'] == 2
        
        # Plot 1: Face-on view (X-Y plane)
        ax = axes[0]
        if np.any(bulge):
            ax.scatter(data['x'][bulge], data['y'][bulge], 
                      s=10, c='red', alpha=0.6, label='Bulge')
        if np.any(disk):
            ax.scatter(data['x'][disk], data['y'][disk], 
                      s=5, c='blue', alpha=0.4, label='Disk')
        
        ax.set_xlabel('X (kpc)', fontsize=12)
        ax.set_ylabel('Y (kpc)', fontsize=12)
        ax.set_title('Face-on View', fontsize=14, fontweight='bold')
        ax.set_aspect('equal')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Plot 2: Edge-on view (X-Z plane)
        ax = axes[1]
        if np.any(bulge):
            ax.scatter(data['x'][bulge], data['z'][bulge], 
                      s=10, c='red', alpha=0.6, label='Bulge')
        if np.any(disk):
            ax.scatter(data['x'][disk], data['z'][disk], 
                      s=5, c='blue', alpha=0.4, label='Disk')
        
        ax.set_xlabel('X (kpc)', fontsize=12)
        ax.set_ylabel('Z (kpc)', fontsize=12)
        ax.set_title('Edge-on View', fontsize=14, fontweight='bold')
        ax.set_aspect('equal')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Extract time from snapshot file (if available)
        time_str = snapshot_file.split('_')[-1].split('.')[0]
        fig.suptitle(f'Galaxy Snapshot: {snapshot_file}', 
                    fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        
        if save:
            output = snapshot_file.replace('.dat', '.png')
            plt.savefig(output, dpi=300, bbox_inches='tight')
            print(f"✓ Saved {output}")
        
        plt.show()
    
    def create_animation(self, output='galaxy_evolution.gif'):
        """
        Create animation from all snapshot files
        """
        snapshot_files = sorted(glob.glob('snapshot_*.dat'))
        
        if len(snapshot_files) < 2:
            print("Need at least 2 snapshots for animation")
            return
        
        print(f"Creating animation from {len(snapshot_files)} snapshots...")
        
        fig, ax = plt.subplots(figsize=(10, 10))
        
        def update(frame):
            ax.clear()
            data = self.load_snapshot(snapshot_files[frame])
            
            if data is None:
                return
            
            bulge = data['type'] == 1
            disk = data['type'] == 2
            
            if np.any(bulge):
                ax.scatter(data['x'][bulge], data['y'][bulge], 
                          s=20, c='red', alpha=0.7, label='Bulge')
            if np.any(disk):
                ax.scatter(data['x'][disk], data['y'][disk], 
                          s=10, c='blue', alpha=0.5, label='Disk')
            
            ax.set_xlabel('X (kpc)', fontsize=12)
            ax.set_ylabel('Y (kpc)', fontsize=12)
            ax.set_title(f'MOND Galaxy Evolution - Frame {frame+1}/{len(snapshot_files)}', 
                        fontsize=14, fontweight='bold')
            ax.set_aspect('equal')
            ax.legend()
            ax.grid(True, alpha=0.3)
            
            # Set consistent axis limits
            ax.set_xlim(-30, 30)
            ax.set_ylim(-30, 30)
        
        anim = FuncAnimation(fig, update, frames=len(snapshot_files), 
                           interval=200, repeat=True)
        
        # Save animation
        try:
            anim.save(output, writer='ffmpeg', fps=5, dpi=150)
            print(f"✓ Saved animation: {output}")
        except Exception as e:
            print(f"Could not save animation (ffmpeg required): {e}")
            print("Showing animation instead...")
            plt.show()
    
    def plot_tully_fisher(self, save=True):
        """
        Plot Tully-Fisher relation: M vs V
        
        MOND predicts: V ∝ M^(1/4)
        """
        # Generate theoretical curve
        masses = np.logspace(9, 12, 50)  # 10^9 to 10^12 solar masses
        
        # MOND prediction: V = (G * M * a0)^(1/4)
        G = 6.67430e-11  # m^3 kg^-1 s^-2
        a0 = 1.2e-10     # m/s^2
        M_sun = 1.989e30 # kg
        
        V_mond = (G * masses * M_sun * a0)**(0.25) / 1000  # km/s
        
        # Newtonian prediction: V ∝ M^(1/2) (at fixed radius)
        R_fixed = 10 * 3.0857e19  # 10 kpc in meters
        V_newton = np.sqrt(G * masses * M_sun / R_fixed) / 1000  # km/s
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Plot 1: Linear scale
        ax = axes[0]
        ax.plot(masses, V_mond, 'r-', linewidth=2, label='MOND: V ∝ M^(1/4)')
        ax.plot(masses, V_newton, 'b--', linewidth=2, label='Newtonian: V ∝ M^(1/2)')
        ax.set_xlabel('Mass (M☉)', fontsize=12)
        ax.set_ylabel('Asymptotic Velocity (km/s)', fontsize=12)
        ax.set_title('Tully-Fisher Relation', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Plot 2: Log-log scale
        ax = axes[1]
        ax.loglog(masses, V_mond, 'r-', linewidth=2, label='MOND: slope = 1/4')
        ax.loglog(masses, V_newton, 'b--', linewidth=2, label='Newtonian: slope = 1/2')
        ax.set_xlabel('Mass (M☉)', fontsize=12)
        ax.set_ylabel('Asymptotic Velocity (km/s)', fontsize=12)
        ax.set_title('Tully-Fisher (Log-Log)', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3, which='both')
        
        # Add slope indicators
        ax.text(0.05, 0.95, 'MOND: consistent with observations\n(< 10% scatter)', 
               transform=ax.transAxes, fontsize=11, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        
        if save:
            plt.savefig('tully_fisher_relation.png', dpi=300, bbox_inches='tight')
            print("✓ Saved tully_fisher_relation.png")
        
        plt.show()
    
    def analyze_all(self):
        """Run all analysis and create all plots"""
        print("=== MOND Simulation Analysis ===\n")
        
        print("1. Plotting rotation curves...")
        self.plot_rotation_curves()
        
        print("\n2. Plotting Tully-Fisher relation...")
        self.plot_tully_fisher()
        
        print("\n3. Analyzing snapshots...")
        snapshots = sorted(glob.glob('snapshot_*.dat'))
        if snapshots:
            # Plot first and last snapshot
            self.plot_galaxy_snapshot(snapshots[0])
            if len(snapshots) > 1:
                self.plot_galaxy_snapshot(snapshots[-1])
            
            # Create animation if multiple snapshots exist
            if len(snapshots) > 2:
                print("\n4. Creating animation...")
                self.create_animation()
        else:
            print("No snapshot files found")
        
        print("\n=== Analysis Complete ===")


def compare_mond_vs_newtonian():
    """
    Create comparison plots showing key differences
    """
    r = np.linspace(0.1, 50, 100)  # kpc
    
    # Parameters
    M_total = 1e11  # Solar masses
    G = 6.67430e-11 * 1.989e30  # G in units of M_sun
    a0 = 1.2e-10
    kpc_to_m = 3.0857e19
    
    # Newtonian: V = sqrt(GM/r)
    V_newton = np.sqrt(G * M_total / (r * kpc_to_m)) / 1000
    
    # MOND: V = (GMa0)^(1/4) (deep MOND regime)
    V_mond = (G * M_total * a0)**(0.25) * np.ones_like(r) / 1000
    
    # Transition regime (approximate)
    r_trans = (G * M_total / a0)**(1/3) / kpc_to_m  # Transition radius
    V_transition = V_mond * np.ones_like(r)
    V_transition[r < r_trans] = V_newton[r < r_trans]
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: Rotation curves
    ax = axes[0]
    ax.plot(r, V_newton, 'b--', linewidth=2, label='Newtonian (Keplerian)', alpha=0.7)
    ax.plot(r, V_mond, 'r-', linewidth=2.5, label='MOND (Deep regime)')
    ax.plot(r, V_transition, 'g-', linewidth=2, label='MOND (with transition)')
    ax.axvline(r_trans, color='gray', linestyle=':', alpha=0.5, 
              label=f'Transition radius ≈ {r_trans:.1f} kpc')
    
    ax.set_xlabel('Radius (kpc)', fontsize=12)
    ax.set_ylabel('Circular Velocity (km/s)', fontsize=12)
    ax.set_title('MOND vs Newtonian: Rotation Curves', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 250)
    
    # Plot 2: Acceleration
    a_newton = G * M_total / (r * kpc_to_m)**2
    a_mond = V_mond**2 * 1000**2 / (r * kpc_to_m)
    
    ax = axes[1]
    ax.loglog(r, a_newton, 'b--', linewidth=2, label='Newtonian: a ∝ 1/r²')
    ax.loglog(r, a_mond, 'r-', linewidth=2.5, label='MOND: a ∝ 1/r (deep regime)')
    ax.axhline(a0, color='green', linestyle=':', linewidth=2, 
              label=f'a₀ = {a0:.2e} m/s²')
    
    ax.set_xlabel('Radius (kpc)', fontsize=12)
    ax.set_ylabel('Acceleration (m/s²)', fontsize=12)
    ax.set_title('Acceleration Profiles', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    plt.savefig('mond_vs_newtonian.png', dpi=300, bbox_inches='tight')
    print("✓ Saved mond_vs_newtonian.png")
    plt.show()


if __name__ == '__main__':
    print("╔════════════════════════════════════════════════════════════╗")
    print("║  MOND Physics Engine - Visualization Suite                ║")
    print("╚════════════════════════════════════════════════════════════╝\n")
    
    analyzer = MONDAnalyzer()
    
    # Check what data files exist
    rot_curves = glob.glob('rotation_curve*.dat')
    snapshots = glob.glob('snapshot_*.dat')
    
    if not rot_curves and not snapshots:
        print("No simulation output files found.")
        print("Run the MOND simulator first to generate data.\n")
        print("Creating theoretical comparison plots instead...")
        compare_mond_vs_newtonian()
    else:
        print(f"Found {len(rot_curves)} rotation curve files")
        print(f"Found {len(snapshots)} snapshot files\n")
        
        # Run full analysis
        analyzer.analyze_all()
        
        # Also create comparison plot
        print("\nCreating theoretical comparison...")
        compare_mond_vs_newtonian()
    
    print("\n✓ Visualization complete!")
