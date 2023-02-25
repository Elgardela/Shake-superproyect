import numpy as np
import matplotlib.pyplot as plt

def plot_energies(steps_array: np.ndarray, ekin_array, epot_array, etot_array) -> None:

    fig, ax = plt.subplots()

    ax.set_xlabel('Simulation step')
    ax.set_ylabel('Energy')

    ax.plot(steps_array, ekin_array, label='Ekin')
    ax.plot(steps_array, epot_array, label='Epot')
    ax.plot(steps_array, etot_array, label='Etot')

    plt.legend()

    plt.show()
    plt.close()

    return None

def plot_temperature(steps_array, temp_array) -> None:

    fig, ax = plt.subplots()

    ax.set_xlabel('Simulation step')
    ax.set_ylabel('Temperature')

    ax.plot(steps_array, temp_array, label='Ekin')

    plt.legend()

    plt.show()

    return None

def plot_lambda(steps_array, lambda_array) -> None:

    fig, ax = plt.subplots()

    ax.set_xlabel('Simulation step')
    ax.set_ylabel('Temperature')

    ax.plot(steps_array, lambda_array, label='Ekin')

    plt.legend()

    plt.show()

    return None

if __name__ == '__main__':

    # Loading data from the output file
    data = np.loadtxt('thermdata.out', dtype=np.float64, skiprows=1)

    # Separating observables
    steps: np.ndarray = data[:, 0]

    ekinetik: np.ndarray = data[:, 1]
    epotential: np.ndarray = data[:, 2]
    etotal: np.ndarray = data[:, 3]
    temperature: np.ndarray = data[:, 4]
    lambda_val: np.ndarray = data[:, 5]

    plot_energies(steps, ekinetik, epotential, etotal)
    plot_temperature(steps, temperature)
    plot_lambda(steps, lambda_val)
