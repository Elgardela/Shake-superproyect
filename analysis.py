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

def plot_atom_distance(mol_array, dist_array, std_array):
    
    fig, ax = plt.subplots()
    
    if not RED_UNITS:
        ax.set_xlabel('Molecule')
        ax.set_ylabel(r'Distance $r$[$\AA$]')
        
        
        ax.errorbar(mol_array, dist_array*SIGMA, yerr=std_array*SIGMA, fmt=' ', alpha=0.7, capsize=2)
        ax.plot(mol_array, dist_array*SIGMA, '.')
        ax.axhline(y=3.0, linestyle='--', linewidth=0.5, color='black', label=r'$3\AA$')
        
    else: 
         ax.set_xlabel('Molecule')
         ax.set_ylabel('Distance')
        
         ax.errorbar(mol_array, dist_array, yerr=std_array, fmt=' ', alpha=0.7, capsize=2)
         ax.plot(mol_array, dist_array, '.')
         ax.axhline(y=3.0/SIGMA, linestyle='--', linewidth=0.5, color='black', label=r'$3\AA$')

         
    ax.grid(alpha=0.5, linestyle=':')
    fig.tight_layout()
    ax.legend(title=r'Tolerance= $10^{10} \AA$', loc='best')
    
    fig.savefig('figs/dist_atoms.eps')
    plt.close()
    return None

def plot_torque(torque_array):
    my_cmap = plt.cm.inferno
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    nrow = len(torque_array[:,0])
    yticks = np.arange(nrow)

    if not RED_UNITS: torque_array *= EPSILON*1e3

    xbins = np.linspace(torque_array.min(), torque_array.max(), 50)

    xcenter = np.convolve(xbins, np.ones(2), "valid")/2
    xwidth = np.diff(xbins)

    for i, ytick in enumerate(yticks):

    #extract the current column from your df by its number
        row =  torque_array[ytick,:]

    #determine the histogram values, here you have to adapt it to your needs
        histvals, _ = np.histogram(row, bins=xbins)

    #plot the histogram as a bar for each bin
    #now with continuous color mapping and edgecolor, but thinner lines, so we can better see all bars
        ax.bar(left=xcenter, height=histvals, width=xwidth, zs=ytick, zdir="y", color=my_cmap(i/nrow), alpha=0.666, edgecolor="grey", linewidth=0.1)


    ax.set_xlabel(r'Torque [N$\cdot$m]')
    ax.set_ylabel('Time steps')
    ax.set_zlabel('Frequency')
    fig.savefig('figs/hist_torque.eps')
    plt.close()

    return None


def plot_angular_mom(ang_mom_array):
    my_cmap = plt.cm.inferno
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    nrow = len(ang_mom_array[:,0])
    yticks = np.arange(nrow)

    xbins = np.linspace(ang_mom_array.min(), ang_mom_array.max(), 50)

    xcenter = np.convolve(xbins, np.ones(2), "valid")/2
    xwidth = np.diff(xbins)

    for i, ytick in enumerate(yticks):

    #extract the current column from your df by its number
        row =  ang_mom_array[ytick,:]

    #determine the histogram values, here you have to adapt it to your needs
        histvals, _ = np.histogram(row, bins=xbins)

    #plot the histogram as a bar for each bin
    #now with continuous color mapping and edgecolor, but thinner lines, so we can better see all bars
        ax.bar(left=xcenter, height=histvals, width=xwidth, zs=ytick, zdir="y", color=my_cmap(i/nrow), alpha=0.666, edgecolor="grey", linewidth=0.1)


    ax.set_xlabel(r'Angular momentum [kg m$^2$/s]')
    ax.set_ylabel('Time steps')
    ax.set_zlabel('Frequency')
    fig.savefig('figs/hist_ang_mom.eps')
    plt.close()

    return None

def plot_ang_mom_torque(torque_array, ang_mom_array, steps_array):

    colors = ['#a6cee3','#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#49006a', '#dd3497']
    fig, ax1 = plt.subplots()

    ax2 = ax1.twinx()

    for i in range(6,7):
        ax1.plot(steps_array*utemps, torque_array[:,i]*EPSILON*1e3, color=colors[i-3], label='Torque')
        ax2.plot(steps_array[1:]*utemps, ang_mom_array[:,i], color=colors[i], linestyle= 'dashed', label='Angular momentum')
    
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel(r'Torque [N$\cdot$m]')
    ax2.set_ylabel('Angular momentum')

    fig.legend()
    fig.savefig('figs/angmom_torque.eps')
    plt.close()
    
    return None

if __name__ == '__main__':

    # Loading data from the output file
    data = np.loadtxt('thermdata.out', dtype=np.float64, skiprows=1)
    data_g_r = np.loadtxt('radial_func.out', dtype=np.float64, skiprows=1)

    data_dist_atoms = np.loadtxt('distance_atoms.data', dtype=np.float64, skiprows=1)

    #SHAKE_iterations = np.loadtxt('SHAKE_iters.out', dtype=np.uint16)

    torque = np.loadtxt('torque.data', dtype=np.float64, skiprows=1)

    ang_mom = np.loadtxt('angular_mom.data', dtype=np.float64, skiprows=1)    

    # Separating observables
    steps: np.ndarray = data[:, 0]

    ekinetik: np.ndarray = data[:, 1]
    epotential: np.ndarray = data[:, 2]
    etotal: np.ndarray = data[:, 3]
    temperature: np.ndarray = data[:, 4]
    lambda_val: np.ndarray = data[:, 5]
    
    radi: np.ndarray = data_g_r[:, 0]
    gdr: np.ndarray = data_g_r[:, 1]

    mol: np.ndarray = data_dist_atoms[:, 0]
    dist: np.ndarray = data_dist_atoms[:, 1]
    std_dist: np.ndarray = data_dist_atoms[:, 2]

    # plot_energies(steps, ekinetik, epotential, etotal)
    # plot_temperature(steps, temperature)
    # plot_lambda(steps, lambda_val)
    # plot_gdr(radi, gdr)
    # plot_atom_distance(mol, dist, std_dist)
    # plot_torque(torque)
    # plot_angular_mom(ang_mom)
    plot_ang_mom_torque(torque, ang_mom, steps)

    #print(f"Mean iterations: {np.mean(SHAKE_iterations)} +- {np.std(SHAKE_iterations)}")
