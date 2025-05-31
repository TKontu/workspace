# Lubrication Piping Flow Distribution Calculation Tool

import numpy as np

from scipy.optimize import fsolve
import matplotlib.pyplot as plt


class Pipe:
    def __init__(self, length, diameter):
        self.length = length  # Length of the pipe (m)
        self.diameter = diameter  # Diameter of the pipe (m)


class Nozzle:
    def __init__(self, restriction_factor=1.0):
        self.restriction_factor = restriction_factor


class Branch:
    def __init__(self, name, pipe=None, nozzle=None):
        self.name = name
        self.pipe = pipe  # Pipe object
        self.nozzle = nozzle  # Nozzle object or None if no nozzle

        self.flow_rate = 0.0  # Flow rate through this branch (m^3/s)
        self.pressure_drop = 0.0  # Pressure drop along the branch (Pa)


class PipingNetwork:
    def __init__(self):
        self.branches = []  # List of Branch objects
        self.total_flow_rate = 0.0  # Total flow rate into the network (m^3/s)
        self.pressure_inlet = 0.0  # Inlet pressure (Pa)

    def add_branch(self, branch):
        """Add a branch to the network."""
        self.branches.append(branch)


    def calculate_pressure_drop(self, branch):
        """
        Calculate pressure drop along a pipe using Darcy-Weisbach equation.
        :param branch: Branch object
        :return: Pressure drop (Pa)
        """
        # Assuming laminar flow for simplicity
        # Reynolds number calculation needed for real-world applications
        
        # For now, let's assume a constant friction factor f = 0.02 for simplicity
        f = 0.02
        
        # Darcy-Weisbach equation: ΔP = f * (L/D) * (ρ * v^2 / 2)
        # Where:
        # ΔP = pressure drop
        # L = length of pipe
        # D = diameter of pipe
        # ρ = density of fluid
        # v = velocity
        
        # We need to calculate the flow rate through this branch first
        if not hasattr(branch.pipe, 'length') or not hasattr(branch.pipe, 'diameter'):
            raise ValueError("Pipe object must have length and diameter attributes")
        
        # Assuming constant density for simplicity (ρ)
        rho = 850  # Density of lubricant in kg/m^3
        
        # Calculate velocity
        if branch.flow_rate == 0:
            return 0.0
            
        A = np.pi * (branch.pipe.diameter ** 2) / 4  # Cross-sectional area
        v = branch.flow_rate /A  # Velocity (m/s)
        
        # Pressure drop calculation using Darcy-Weisbach equation
        pressure_drop = f * (branch.pipe.length /branch.pipe.diameter) * (rho * v**2 / 2)
        return pressure_drop


    def solve_network(self):
        """
        Solve the entire piping network to find flow rates in each branch.
        This is a simplified version that assumes steady state and no branching
        for now. We'll need to extend it to handle branching and nozzles.
        later.
        """
        
        # For simplicity, let's assume we have a simple linear system with one inlet
        # and multiple outlets
        
        # First, calculate the total resistance of the network
        total_resistance = 0.0
        
        for branch in self.branchbranches:
            if hasattr(branch.pipe, 'length') and hasattr(branch.pipe, 'diameter'):
                # Calculate pressure drop for each branch
                # For now, we'll just sum up all the pressure drops to get a rough estimate
                total_resistance += self.calculate_pressure_drop(branch)
        
        return total_resistance


    def update_flow_rates(self):
        """
        Update flow rates in the network based on current conditions.
        This is a placeholder for more complex calculations that would involve
        solving a system of equations representing the entire network
        """
        # For now, we'll just distribute the total flow rate equally among all branches
        In reality, this will need to be solved as a system of nonlinear equations
        
        if not self.branches:
            return
            
        equal_flow = self.total_flow_rate / len(self.branches)
        
        for branch in self.branches:
            branch.flow_rate = equal_flow


# Example usage
if __name__ == "__main__":
    # Create a simple network with two branches
    pipe1 = Pipe(length=10.0, diameter=0.05)    nozzle1 = Nozzle(restriction_factor=0.8)
    
    pipe2 = Pipe(length=15.0, diameter=0.03)
    nozzle2 = Nozzle(restriction_factor=0.6)

    
    branch1 = Branch(name="Branch 1", pipe=pipe1, nozzle=nozzle1)
    branch2 = Branch(name="Branch 2", pipe=pipe2, nozzle=nozzle2)
    
    network = PipingNetwork()
    network.add_branch(branch1)
    network.add_branch(branch2)
    
    # Set total flow rate and inlet pressure
    network.total_flow_rate = 0.05  # m^3/s
    
    # Update flow rates (this is a placeholder for the actual calculation)
    network.update_flow_rates()
    
    print("Flow Rates:")
    for branch in network.branches:
        print(f"{branch.name}: {branch.flow_rate:.4f} m³/s")