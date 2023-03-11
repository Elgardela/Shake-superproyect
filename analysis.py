import numpy as np
import matplotlib.pyplot as plt

# Parameters
SIGMA: float = 3.405  # [A]
EPSILON: float = 119.8 # [K]
KB: float = 1.380649e-20  # [kJ/K]
R_G: float = 8.314472673 #J/(mol*K)
MASS: float = 39.948
utemps = SIGMA * np.sqrt(MASS/EPSILON) * np.sqrt(10.0/R_G)

EPSILON *= 1.380649e-20  # [kJ]

RED_UNITS: bool = False  # Compute in reduced units

def plot_energies(steps_array: np.ndarray, ekin_array, epot_array, etot_array) -> None:

    fig, ax = plt.subplots()

    if not RED_UNITS:
        ax.set_xlabel(r'Simulation time [ps]')
        ax.set_ylabel(r'Energy [kJ]')

        ax.plot(steps_array*utemps, ekin_array*EPSILON, label=r'$E_{kin}$', color='red')
        ax.plot(steps_array*utemps, epot_array*EPSILON, label=r'$E_{pot}$', color='green')
        ax.plot(steps_array*utemps, etot_array*EPSILON, label=r'$E_{tot}$', color='black')
    else:
        ax.set_xlabel(r'Simulation time [-]')
        ax.set_ylabel(r'Energy [-]')

        ax.plot(steps_array, ekin_array, label='Ekin', color='red')
        ax.plot(steps_array, epot_array, label='Epot', color='green')
        ax.plot(steps_array, etot_array, label='Etot', color='black')

    fig.legend()
    ax.grid(alpha=0.5, linestyle=':')

    fig.tight_layout()

    fig.savefig('figs/energy.eps')
    plt.close()

    return None

def plot_temperature(steps_array, temp_array) -> None:

    N: int = 110

    fig, ax = plt.subplots()

    if not RED_UNITS:
        ax.set_xlabel(r'Simulation time [ps]')
        ax.set_ylabel(r'Temperature [K]')
        ax.plot(steps_array*utemps, temp_array*EPSILON/KB, label="Sampled data")
        conv_array = np.convolve(temp_array*EPSILON/KB, np.ones(N)/N, mode='same')
        time_arr = steps_array*utemps
        ax.plot(time_arr[conv_array > 220.0], conv_array[conv_array > 220.0], label=f'Running average. {N=}')
    else:
        ax.set_xlabel(r'Simulation time [-]')
        ax.set_ylabel(r'Temperature [-]')
        ax.plot(steps_array, temp_array)

    ax.set_ylim([210., 240.])
    ax.grid(alpha=0.5, linestyle=':')
    fig.legend()
    fig.tight_layout()

    fig.savefig('figs/temperature.eps')
    plt.close()

    return None

def plot_lambda(steps_array, lambda_array) -> None:

    fig, ax = plt.subplots()

    if not RED_UNITS:
        ax.set_xlabel(r'Simulation time [ps]')
        ax.plot(steps_array, lambda_array)
    else:
        ax.set_xlabel(r'Simulation time [ps]')
        ax.plot(steps_array*utemps, lambda_array)
    
    ax.set_ylabel('Berendsen Lambda [-]')

    ax.grid(alpha=0.5, linestyle=':')
    fig.tight_layout()
    
    fig.savefig('figs/b_lambda_t.eps')
    plt.close()

    return None

def plot_gdr(radi_array, gdr_array) -> None:

    fig, ax = plt.subplots()

    if not RED_UNITS:
        ax.set_xlabel(r'Distance $r$[$\AA$]')
        ax.set_ylabel(r'Radial distribution function, $g(r)$')
        
        ax.plot(radi_array*SIGMA, gdr_array)
    else: 
        ax.set_xlabel(r'Distance $r$ [-]')
        ax.set_ylabel(r'Radial distribution function, $g(r)$')
        
        ax.plot(radi_array, gdr_array)

    ax.grid(alpha=0.5, linestyle=':')
    ax.axhline(y=1.0, linestyle='--', linewidth=0.5, color='black')

    print(f"RDF peak = {radi_array[gdr_array == np.max(gdr_array)][0]*SIGMA}")

    fig.tight_layout()

    fig.savefig('figs/rdf.eps')

    plt.close()

    return None

if __name__ == '__main__':

    # Loading data from the output file
    data = np.loadtxt('thermdata.out', dtype=np.float64, skiprows=1)
    data_g_r = np.loadtxt('radial_func.out', dtype=np.float64, skiprows=1)

    SHAKE_iterations = np.loadtxt('SHAKE_iters.out', dtype=np.uint16)

    # Separating observables
    steps: np.ndarray = data[:, 0]

    ekinetik: np.ndarray = data[:, 1]
    epotential: np.ndarray = data[:, 2]
    etotal: np.ndarray = data[:, 3]
    temperature: np.ndarray = data[:, 4]
    lambda_val: np.ndarray = data[:, 5]
    
    radi: np.ndarray = data_g_r[:, 0]
    gdr: np.ndarray = data_g_r[:, 1]

    # plot_energies(steps, ekinetik, epotential, etotal)
    # plot_temperature(steps, temperature)
    # plot_lambda(steps, lambda_val)
    # plot_gdr(radi, gdr)

    print(f"Mean iterations: {np.mean(SHAKE_iterations)} +- {np.std(SHAKE_iterations)}")
