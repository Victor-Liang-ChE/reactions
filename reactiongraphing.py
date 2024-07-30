#%% Make chemical kinetics solver for multiple elementary reactions
# Reminder: Elementary reactions are reactions where the coefficients of the reactants are the same as the order of the reaction
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import re

def reactiongraphing(reactions, ks, C0):
    # Check that the number of ks elements matches the number of reactions
    if len(ks) != len(reactions):
        raise ValueError("The number of rate constants does not match the number of reactions.")

    # Function to remove numerical coefficients from species names
    def remove_coefficients(species):
        return re.sub(r'^\d*\*?', '', species)

    # Extract unique species from reactions in the order they appear
    ordered_species = []
    unique_species = set()
    for reaction in reactions:
        reactants, products = reaction.split('=')
        reactant_species = [remove_coefficients(species) for species in reactants.split('+')]
        product_species = [remove_coefficients(species) for species in products.split('+')]
        for species in reactant_species + product_species:
            if species not in unique_species:
                unique_species.add(species)
                ordered_species.append(species)

    # Check that all unique species are present in C0
    if not unique_species.issubset(C0.keys()):
        missing_species = unique_species - set(C0.keys())
        raise ValueError(f"The following species are missing in C0: {missing_species}")

    # Define the system of ODEs
    def odes(t, y):
        dydt = np.zeros(len(ordered_species))
        concentrations = {species: y[i] for i, species in enumerate(ordered_species)}
        for i, reaction in enumerate(reactions):
            reactants, products = reaction.split('=')
            reactant_species = []
            for species in reactants.split('+'):
                coeff, sp = re.match(r'(\d*)(\w+)', species).groups()
                coeff = int(coeff) if coeff else 1
                reactant_species.append((sp, coeff))
            
            # Calculate the reaction rate
            rate = ks[i] * np.prod([concentrations[sp]**coeff for sp, coeff in reactant_species])
            
            # Update differential forms
            for sp, coeff in reactant_species:
                dydt[ordered_species.index(sp)] -= rate * coeff
            
            for product in products.split('+'):
                coeff, sp = re.match(r'(\d*)(\w+)', product).groups()
                coeff = int(coeff) if coeff else 1
                dydt[ordered_species.index(sp)] += rate * coeff
        
        return dydt

    # Initial concentrations
    y0 = [C0[species] for species in ordered_species]

    # Time span for the simulation
    t_span = (0, 10)
    t_eval = np.linspace(t_span[0], t_span[1], 1000)

    # Solve the ODEs
    solution = solve_ivp(odes, t_span, y0, t_eval=t_eval, method='RK45')

    # Determine the steady state time
    tolerance = 1e-4
    steady_state_time = t_span[1]
    for i in range(1, len(solution.t)):
        if np.all(np.abs(solution.y[:, i] - solution.y[:, i-1]) < tolerance):
            steady_state_time = solution.t[i]
            break

    # Plot the results
    plt.figure(figsize=(10, 6))
    for i, species in enumerate(ordered_species):
        plt.plot(solution.t, solution.y[i], label=f'{species}')

    plt.xlabel('Time', fontsize=16)
    plt.ylabel('Concentration', fontsize=16)
    plt.xlim(0, steady_state_time)
    plt.ylim(0, None)
    plt.legend()
    plt.title('Concentration vs. Time', fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.show()
# %%
reactions = ['A+B=C', 'B+C=D']
ks = [1, 1.5]
C0 = {'A': 1, 'B': 1, 'C': 0, 'D': 0}
reactiongraphing(reactions, ks, C0)
# %%
