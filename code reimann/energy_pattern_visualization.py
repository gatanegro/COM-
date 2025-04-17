#!/usr/bin/env python3
"""
Energy Pattern Visualization for Quantum Measurement in COM Framework

This script implements the visualization of energy patterns as described in the
Continuous Oscillatory Model (COM) framework, specifically for addressing the
quantum measurement problem.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from matplotlib.colors import hsv_to_rgb
import os

# Constants from the COM framework
LZ = 1.23498  # Fundamental scaling constant
HQS = 0.235 * LZ  # Harmonic Quantum Scalar
E_MIN = 0.01  # Minimum energy state (replacing zero)

def octave_reduce(value):
    """
    Reduce a value to its fundamental oscillatory nature through octave reduction.
    
    Parameters:
    - value: The value to reduce
    
    Returns:
    - Reduced value between 1 and 9
    """
    return (value - 1) % 9 + 1

def map_to_octave(value, layer):
    """
    Map a reduced value to an octave structure.
    
    Parameters:
    - value: The reduced value (1-9)
    - layer: The octave layer
    
    Returns:
    - (x, y) coordinates in the octave structure
    """
    angle = (value / 9) * 2 * np.pi
    x = np.cos(angle) * (layer + 1)
    y = np.sin(angle) * (layer + 1)
    return x, y

def create_energy_pattern(modes, amplitudes, phases):
    """
    Create a visual representation of an energy pattern.
    
    Parameters:
    - modes: List of oscillatory modes
    - amplitudes: Energy amplitude in each mode
    - phases: Phase of oscillation in each mode
    
    Returns:
    - 3D coordinates and visual properties for rendering
    """
    coordinates = []
    colors = []
    sizes = []
    
    for i, (mode, amplitude, phase) in enumerate(zip(modes, amplitudes, phases)):
        # Reduce mode to single digit using octave reduction
        reduced_mode = octave_reduce(mode)
        
        # Map to circular octave
        x, y = map_to_octave(reduced_mode, i)
        z = i  # Each mode gets its own layer for clarity
        
        # Scale point size by amplitude
        size = amplitude * 50
        
        # Color based on phase
        hue = phase / (2 * np.pi)
        saturation = 0.8
        value = min(1.0, amplitude / max(amplitudes) if max(amplitudes) > 0 else 0)
        
        coordinates.append((x, y, z))
        colors.append((hue, saturation, value))
        sizes.append(size)
    
    return np.array(coordinates), np.array(colors), np.array(sizes)

def visualize_energy_pattern(coordinates, colors, sizes, title="Energy Pattern Visualization"):
    """
    Visualize an energy pattern in 3D.
    
    Parameters:
    - coordinates: 3D coordinates of energy modes
    - colors: HSV colors for each mode
    - sizes: Sizes of points representing energy amplitude
    - title: Plot title
    
    Returns:
    - Figure and axes objects
    """
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Convert HSV colors to RGB for matplotlib
    rgb_colors = np.array([hsv_to_rgb(color) for color in colors])
    
    # Plot points representing energy modes
    ax.scatter(
        coordinates[:, 0],
        coordinates[:, 1],
        coordinates[:, 2],
        s=sizes,
        c=rgb_colors,
        alpha=0.8
    )
    
    # Connect points with lines to show the pattern structure
    ax.plot(
        coordinates[:, 0],
        coordinates[:, 1],
        coordinates[:, 2],
        'k-',
        alpha=0.3
    )
    
    # Add labels and adjust the view
    ax.set_title(title)
    ax.set_xlabel("X (Horizontal Oscillation)")
    ax.set_ylabel("Y (Vertical Oscillation)")
    ax.set_zlabel("Z (Octave Layer)")
    
    # Set equal aspect ratio
    max_range = np.array([
        coordinates[:, 0].max() - coordinates[:, 0].min(),
        coordinates[:, 1].max() - coordinates[:, 1].min(),
        coordinates[:, 2].max() - coordinates[:, 2].min()
    ]).max() / 2.0
    
    mid_x = (coordinates[:, 0].max() + coordinates[:, 0].min()) * 0.5
    mid_y = (coordinates[:, 1].max() + coordinates[:, 1].min()) * 0.5
    mid_z = (coordinates[:, 2].max() + coordinates[:, 2].min()) * 0.5
    
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    
    return fig, ax

def apply_hamiltonian(energy_pattern, hamiltonian):
    """
    Apply energy transformation operator (Hamiltonian) to an energy pattern.
    
    Parameters:
    - energy_pattern: Dictionary with 'modes', 'amplitudes', and 'phases'
    - hamiltonian: Matrix representing energy transformation
    
    Returns:
    - Updated energy pattern
    """
    new_pattern = energy_pattern.copy()
    
    # Apply Hamiltonian to amplitudes
    new_amplitudes = np.dot(hamiltonian, energy_pattern['amplitudes'])
    
    # Ensure no amplitude goes below minimum energy
    new_amplitudes = np.maximum(new_amplitudes, E_MIN)
    
    # Update phases based on energy transfer
    phase_shifts = np.angle(np.exp(1j * np.array(energy_pattern['phases'])) * 
                           np.exp(1j * np.dot(hamiltonian, np.array(energy_pattern['phases']))))
    
    new_pattern['amplitudes'] = new_amplitudes
    new_pattern['phases'] = [(p + shift) % (2 * np.pi) for p, shift in 
                            zip(energy_pattern['phases'], phase_shifts)]
    
    return new_pattern

def evolve_energy_pattern(energy_pattern, hamiltonian, phase_steps):
    """
    Evolve an energy pattern according to the COM evolution equation.
    
    Parameters:
    - energy_pattern: Initial energy pattern (dict with modes, amplitudes, phases)
    - hamiltonian: Energy transformation operator
    - phase_steps: Number of phase steps to simulate
    
    Returns:
    - Time series of evolved energy patterns
    """
    evolution = [energy_pattern.copy()]
    current_pattern = energy_pattern.copy()
    
    for step in range(phase_steps):
        # Apply energy transformation operator
        new_pattern = apply_hamiltonian(current_pattern, hamiltonian)
        
        evolution.append(new_pattern.copy())
        current_pattern = new_pattern
    
    return evolution

def calculate_synchronization(system_pattern, measurement_pattern):
    """
    Calculate phase synchronization between two energy patterns.
    
    Parameters:
    - system_pattern: Energy pattern of the quantum system
    - measurement_pattern: Energy pattern of the measurement device
    
    Returns:
    - Synchronization matrix
    """
    sync_matrix = np.zeros((len(system_pattern['modes']), len(measurement_pattern['modes'])))
    
    for i, (s_amp, s_phase) in enumerate(zip(system_pattern['amplitudes'], system_pattern['phases'])):
        for j, (m_amp, m_phase) in enumerate(zip(measurement_pattern['amplitudes'], measurement_pattern['phases'])):
            # Calculate phase difference
            phase_diff = (s_phase - m_phase) % (2 * np.pi)
            
            # Synchronization is strongest when phases align
            sync_strength = s_amp * m_amp * np.cos(phase_diff)
            sync_matrix[i, j] = sync_strength
    
    return sync_matrix

def apply_interaction(system_pattern, measurement_pattern, sync_matrix, interaction_strength):
    """
    Update both patterns based on their interaction.
    
    Parameters:
    - system_pattern: Energy pattern of the quantum system
    - measurement_pattern: Energy pattern of the measurement device
    - sync_matrix: Synchronization matrix between patterns
    - interaction_strength: Strength of coupling between patterns
    
    Returns:
    - Updated system and measurement patterns
    """
    new_system = system_pattern.copy()
    new_measurement = measurement_pattern.copy()
    
    # Calculate energy transfer based on synchronization
    system_transfer = np.zeros(len(system_pattern['amplitudes']))
    measurement_transfer = np.zeros(len(measurement_pattern['amplitudes']))
    
    for i in range(len(system_pattern['amplitudes'])):
        for j in range(len(measurement_pattern['amplitudes'])):
            # Energy flows from system to measurement and vice versa based on synchronization
            transfer = interaction_strength * sync_matrix[i, j]
            
            # Energy conservation: what flows out of one flows into the other
            system_transfer[i] -= transfer
            measurement_transfer[j] += transfer
    
    # Apply energy transfers
    new_system['amplitudes'] = np.maximum(
        np.array(system_pattern['amplitudes']) + system_transfer, 
        E_MIN
    )
    
    new_measurement['amplitudes'] = np.maximum(
        np.array(measurement_pattern['amplitudes']) + measurement_transfer,
        E_MIN
    )
    
    # Update phases based on energy transfer
    for i in range(len(system_pattern['amplitudes'])):
        # Phase shifts toward measurement modes with strongest synchronization
        strongest_sync_idx = np.argmax(np.abs(sync_matrix[i, :]))
        phase_target = measurement_pattern['phases'][strongest_sync_idx]
        
        # Phase shift is proportional to synchronization strength and interaction strength
        phase_shift = interaction_strength * np.max(np.abs(sync_matrix[i, :])) * 0.1
        
        # Move phase toward target
        current_phase = system_pattern['phases'][i]
        phase_diff = (phase_target - current_phase) % (2 * np.pi)
        if phase_diff > np.pi:
            phase_diff -= 2 * np.pi
        
        new_system['phases'][i] = (current_phase + phase_shift * np.sign(phase_diff)) % (2 * np.pi)
    
    # Similar phase updates for measurement pattern
    for j in range(len(measurement_pattern['amplitudes'])):
        strongest_sync_idx = np.argmax(np.abs(sync_matrix[:, j]))
        phase_target = system_pattern['phases'][strongest_sync_idx]
        
        phase_shift = interaction_strength * np.max(np.abs(sync_matrix[:, j])) * 0.1
        
        current_phase = measurement_pattern['phases'][j]
        phase_diff = (phase_target - current_phase) % (2 * np.pi)
        if phase_diff > np.pi:
            phase_diff -= 2 * np.pi
        
        new_measurement['phases'][j] = (current_phase + phase_shift * np.sign(phase_diff)) % (2 * np.pi)
    
    return new_system, new_measurement

def simulate_measurement_interaction(system_pattern, measurement_pattern, interaction_strength, phase_steps):
    """
    Simulate the interaction between a quantum system and measurement device.
    
    Parameters:
    - system_pattern: Energy pattern of the quantum system
    - measurement_pattern: Energy pattern of the measurement device
    - interaction_strength: Strength of coupling between patterns
    - phase_steps: Number of phase steps to simulate
    
    Returns:
    - Time series of both patterns during interaction
    """
    system_evolution = [system_pattern.copy()]
    measurement_evolution = [measurement_pattern.copy()]
    
    current_system = system_pattern.copy()
    current_measurement = measurement_pattern.copy()
    
    for step in range(phase_steps):
        # Calculate phase synchronization between patterns
        sync_matrix = calculate_synchronization(current_system, current_measurement)
        
        # Update both patterns based on interaction
        new_system, new_measurement = apply_interaction(
            current_system, 
            current_measurement,
            sync_matrix,
            interaction_strength
        )
        
        system_evolution.append(new_system.copy())
        measurement_evolution.append(new_measurement.copy())
        
        current_system = new_system
        current_measurement = new_measurement
    
    return system_evolution, measurement_evolution

def visualize_collapse_process(system_evolution):
    """
    Visualize the apparent collapse process during measurement.
    
    Parameters:
    - system_evolution: Time series of system energy patterns during measurement
    
    Returns:
    - Figure showing energy concentration over time
    """
    # Track energy in each mode over time
    mode_energies = []
    
    for pattern in system_evolution:
        mode_energies.append(pattern['amplitudes'])
    
    # Convert to numpy array for easier analysis
    mode_energies = np.array(mode_energies)
    
    # Calculate entropy of energy distribution over time
    entropy = []
    for distribution in mode_energies:
        normalized = distribution / np.sum(distribution)
        ent = -np.sum([p * np.log(p) if p > 0 else 0 for p in normalized])
        entropy.append(ent)
    
    # Create visualization
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Plot energy in each mode over time
    for i in range(mode_energies.shape[1]):
        ax1.plot(mode_energies[:, i], label=f"Mode {i+1}")
    
    ax1.set_xlabel("Phase Steps")
    ax1.set_ylabel("Energy Amplitude")
    ax1.set_title("Energy Redistribution During Measurement")
    ax1.legend()
    
    # Plot entropy over time
    ax2.plot(entropy)
    ax2.set_xlabel("Phase Steps")
    ax2.set_ylabel("Entropy")
    ax2.set_title("Entropy of Energy Distribution")
    
    plt.tight_layout()
    return fig

def animate_measurement_process(system_evolution, measurement_evolution, interval=200):
    """
    Create an animation of the measurement process.
    
    Parameters:
    - system_evolution: Time series of system energy patterns
    - measurement_evolution: Time series of measurement device patterns
    - interval: Time between frames in milliseconds
    
    Returns:
    - Animation object
    """
    fig = plt.figure(figsize=(15, 8))
    
    # Create two subplots for system and measurement device
    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122, projection='3d')
    
    # Function to update the plot for each animation frame
    def update(frame):
        ax1.clear()
        ax2.clear()
        
        # Get current patterns
        system = system_evolution[frame]
        measurement = measurement_evolution[frame]
        
        # Create visual representations
        sys_coords, sys_colors, sys_sizes = create_energy_pattern(
            system['modes'], system['amplitudes'], system['phases']
        )
        
        meas_coords, meas_colors, meas_sizes = create_energy_pattern(
            measurement['modes'], measurement['amplitudes'], measurement['phases']
        )
        
        # Convert HSV colors to RGB
        sys_rgb = np.array([hsv_to_rgb(color) for color in sys_colors])
        meas_rgb = np.array([hsv_to_rgb(color) for color in meas_colors])
        
        # Plot system
        ax1.scatter(
            sys_coords[:, 0], sys_coords[:, 1], sys_coords[:, 2],
            s=sys_sizes, c=sys_rgb, alpha=0.8
        )
        ax1.plot(
            sys_coords[:, 0], sys_coords[:, 1], sys_coords[:, 2],
            'k-', alpha=0.3
        )
        
        # Plot measurement device
        ax2.scatter(
            meas_coords[:, 0], meas_coords[:, 1], meas_coords[:, 2],
            s=meas_sizes, c=meas_rgb, alpha=0.8
        )
        ax2.plot(
            meas_coords[:, 0], meas_coords[:, 1], meas_coords[:, 2],
            'k-', alpha=0.3
        )
        
        # Set titles and labels
        ax1.set_title(f"Quantum System (Step {frame})")
        ax2.set_title(f"Measurement Device (Step {frame})")
        
        for ax in [ax1, ax2]:
            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            ax.set_zlabel("Z")
            
            # Set consistent view limits
            ax.set_xlim(-10, 10)
            ax.set_ylim(-10, 10)
            ax.set_zlim(0, len(system['modes']))
    
    # Create animation
    anim = animation.FuncAnimation(
        fig, update, frames=len(system_evolution),
        interval=interval, blit=False
    )
    
    return anim

def is_even_resonant(amplitude):
    """
    Determine if an energy amplitude is even-resonant.
    
    Parameters:
 
(Content truncated due to size limit. Use line ranges to read in chunks)