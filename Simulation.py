import random
import matplotlib.pyplot as plt
import numpy as np
import math

# Input parameters ####################

# time steps
dt = 0.01
timesteps = math.floor(1000 / dt)

# Simulation interval
sim_interval = 10000

# Simulation start and end times
sim_start = timesteps - sim_interval
sim_end = timesteps


########## Functions ############
def get_input_molecules(X, alpha, beta):
    dt_ = 0.01
    if X == 0:
        beta_dt = 0
        alpha_dt = alpha
    elif X > 0:
        beta_dt = beta * X * dt_
        alpha_dt = alpha * dt_
    p = round(random.uniform(0, 1), 2)
    if p <= alpha_dt:
        print("Birth")
        X += 1
    elif alpha_dt < p <= (alpha_dt + beta_dt):
        print("Death")
        X -= 1
    else:
        print("Do nothing: Stay")
    return X


def plot_state_diagram(state_data):
    plt.figure(figsize=[5, 1.5])
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False
    plt.xlabel('Time (t)')
    plt.ylabel('State')
    y = [0, 1]
    y_labels = ['off', 'on']
    plt.yticks(y, y_labels)
    plt.xticks([])
    plt.step(state_data[:, 0], state_data[:, 1])
    plt.show()


def plot_fig1b(Kon_arr, p_tau, P_tau_analytical):
    plt.style.use('seaborn-darkgrid')
    plt.xlabel('t')
    plt.ylabel(r'$P(\tau > t)$')
    plt.yscale('log')
    plt.ylim(10e-10, 10e0)
    palette = plt.get_cmap('Set1')
    plt.scatter(range(len(p_tau[Kon_arr[0]])), p_tau[Kon_arr[0]], label='kon = koff = 0.05', marker='o',
                color=palette(0), facecolors='none')
    plt.scatter(range(len(p_tau[Kon_arr[1]])), p_tau[Kon_arr[1]], label='kon = koff = 0.1', marker='^',
                color=palette(1))
    plt.scatter(range(len(p_tau[Kon_arr[2]])), p_tau[Kon_arr[2]], label='kon = koff = 0.2', marker='s',
                color=palette(2))
    plt.plot(range(len(P_tau_analytical[Kon_arr[0]])), P_tau_analytical[Kon_arr[0]], label='kon = koff = 0.05',
             color=palette(0))
    plt.plot(range(len(P_tau_analytical[Kon_arr[1]])), P_tau_analytical[Kon_arr[1]], label='kon = koff = 0.1',
             color=palette(1))
    plt.plot(range(len(P_tau_analytical[Kon_arr[2]])), P_tau_analytical[Kon_arr[2]], label='kon = koff = 0.2',
             color=palette(2))
    plt.legend()
    plt.show()


def plot_fig2b(Kon_arr, p_tau, P_analytical, P_analytical_without_noise):
    plt.xlabel('t')
    plt.ylabel(r'$P(\tau > t)$')
    plt.yscale('log')
    plt.ylim(10e-10, 10e0)
    palette = plt.get_cmap('Set1')
    plt.scatter(range(len(p_tau[Kon_arr[0]])), p_tau[Kon_arr[0]], label=r'Simulation', marker='o',
                color=palette(0), facecolors='none')
    plt.plot(range(len(P_analytical[Kon_arr[0]])), P_analytical[Kon_arr[0]], label=r'$exp(-\~k_{on}\alpha t)$',
             color=palette(2))
    plt.plot(range(len(P_analytical_without_noise[Kon_arr[0]])), P_analytical_without_noise[Kon_arr[0]], label=r'$exp(-k_{on}\alpha t)$', color=palette(1))
    plt.legend()
    plt.show()

#########################################


# Main loop
def run_simulation(alpha, beta, Kon_arr, Koff_arr, steps, noise=False, constant_input=False, ):
    # Initial number of input molecules
    X = 0

    # Initial state
    state = 0

    # Tuple of tau frequency and duration of tau intervals
    tau_data = {}

    # Probability of tau > t
    p_tau = {}

    # Analytical
    P_tau_analytical = {}

    # Run the simulation for all values of kon and koff
    for idx in range(len(Kon_arr)):
        tau_arr = []
        kon = Kon_arr[idx]
        koff = Koff_arr[idx]
        state_data = []
        tau_start = sim_start
        tau_end = 0
        tau_count = 0
        if noise:
            zs = 1 / (kon + 1)
            kon_noise = zs * kon
            kon_eff = kon_noise
        else:
            kon_eff = kon

        for t in range(timesteps):
            if noise and constant_input:
                # Constant input rate with X(t) = alpha
                X = alpha
            else:
                X = get_input_molecules(X, alpha, beta)

            # print('X(t): {}'.format(X))
            if sim_start <= t < sim_end:
                if X == 0:
                    Kondt = Koffdt = 0
                else:

                    Kondt = kon_eff * X * dt
                    Koffdt = koff * dt

                p = random.uniform(0, dt)

                state_data.append((t, state))
                if state == 0 and p <= Kondt:
                    tau_end = t
                    tau_count += 1
                    tau_arr.append(tau_end - tau_start)
                    state = 1
                elif state == 1 and p <= Koffdt:
                    tau_start = t
                    state = 0
                print(state)
        tau_data[kon] = (tau_count, tau_arr)

    print('tau count: {}'.format(tau_count))
    print('tau intervals: {}'.format(tau_data))
    state_data = np.array(state_data)
    #plot_state_diagram(state_data)

    # Calculate experimental P(tau > t)
    for key in tau_data:
        tau_count = tau_data[key][0]
        p_tau_arr = []
        for i in range(steps):
            tau_freq = sum(j > i for j in tau_data[key][1])
            if tau_freq > 0:
                p_tau_arr.append(tau_freq / tau_count)
        p_tau[key] = p_tau_arr

    # Calculate the analytical value of P(tau > t)
    for i in range(len(Kon_arr)):
        if noise:
            zs = 1 / (Kon_arr[i] + 1)
            kon_noise = zs * Kon_arr[i]
            kon_eff = kon_noise
        else:
            kon_eff = Kon_arr[i]

        P_arr = []
        # Calculate the probability for all the time steps
        for t in range(steps):
            P = np.exp(-kon_eff * alpha * t)
            P_arr.append(P)

        # Add the P values to a dictionary
        P_tau_analytical[Kon_arr[i]] = P_arr
    return p_tau, P_tau_analytical


def main():
    ############### Run 1 ######################
    # float; on rates of the switch
    Kon_arr = [0.05, 0.1, 0.2]

    # float; off rates of the switch
    Koff_arr = [0.05, 0.1, 0.2]

    # Birth and death rate of input molecules
    alpha = 1
    beta = 1
    steps = 120

    # Without noise
    # run_simulation(noise=False)

    # With noise
    p_tau, P_analytical = run_simulation(alpha, beta, Kon_arr, Koff_arr, steps, noise=True, constant_input=True)
    plot_fig1b(Kon_arr, p_tau, P_analytical)

    ##################### Run 2 ####################
    Kon_arr = [1]
    Koff_arr = [1]

    alpha = 1
    beta = 0.2
    steps = 15
    p_tau_without_noise, P_analytical_without_noise = run_simulation(alpha, beta, Kon_arr, Koff_arr, steps)
    p_tau, P_analytical = run_simulation(alpha, beta, Kon_arr, Koff_arr, steps, noise=True, constant_input=True)
    plot_fig2b(Kon_arr, p_tau, P_analytical, P_analytical_without_noise)


if __name__ == "__main__":
    main()
