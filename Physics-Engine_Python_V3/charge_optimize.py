import csv 
import numpy as np 
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import minimize
import sys


input_file = sys.argv[1]
output_file = sys.argv[2]
verbose = sys.argv[3].lower() == 'true' if len(sys.argv) > 3 else False

def read_data_from_csv(filename):
    t_log = []
    energy_log = []
    charged_nodes_positions_log = []
    with open(filename, mode='r') as file:
        reader = csv.reader(file)
        headings = next(reader)  
        for row in reader:
            t_log.append(float(row[0]))
            energy_log.append(float(row[1]))
            num_charges = (len(row) - 2) // 2
            charged_positions = []
            for i in range(num_charges):
                x = float(row[2 + 2*i])
                y = float(row[3 + 2*i])
                charged_positions.append((x, y))
            charged_nodes_positions_log.append(charged_positions)
    return t_log, energy_log, charged_nodes_positions_log



t_log, energy_log, charged_nodes_positions_log = read_data_from_csv(input_file)
t_log = np.array(t_log)
energy_log = np.array(energy_log)
charged_nodes_positions_log = np.array(charged_nodes_positions_log)
if verbose:
    print("Read data from CSV: {}".format(input_file))


minima_indices, _ = find_peaks([-e for e in energy_log])
max_idx = minima_indices[np.argmax(t_log[minima_indices])]
t_log_filtered = t_log[:max_idx+1]
energy_log_filtered = energy_log[:max_idx+1]

maxima_indices, _ = find_peaks(energy_log_filtered)


α = int((max_idx)*1.1)
β = t_log[int(maxima_indices[0]*0.7)]
σ = 0.5*np.abs(α-β)
a = np.max(energy_log_filtered) * 0.001

def gaussian(t, μ, σ, a):
    return a * np.exp(-((t - μ) / σ) ** 2)



def resulting_energy_function_non_interacting(t, β, σ, a):
    val =  -gaussian(t, μ=β, σ=σ, a=a) + gaussian(0, μ=β, σ=σ, a=a)
    return val
final_energy_values = resulting_energy_function_non_interacting(t_log_filtered, β=β, σ=σ, a=a)

required_charge_energy_values = final_energy_values-energy_log_filtered
if verbose:
    print("Computed required charge energy values based on non-interacting model.")


num_charges = charged_nodes_positions_log.shape[1]
def compute_A_all(positions: np.ndarray):
    """
    positions shape: (T, N, dim)
    """
    T, N, _ = positions.shape
    A_all = np.zeros((T, N, N))

    for t in range(T):
        pos = positions[t]

        diff = pos[:, None, :] - pos[None, :, :]
        dist = np.linalg.norm(diff, axis=2)

        with np.errstate(divide='ignore', invalid='ignore'):
            A = np.where(dist > 0, 1.0 / dist, 0.0)
        np.fill_diagonal(A, 0.0)
        A_all[t] = A
    return A_all



def energy(q, A):
    return 0.5 * q @ A @ q
def energy_time_series(q, A_all):
    return np.array([energy(q, A) for A in A_all])
def objective(q, A_all, U_target):
    U_pred = energy_time_series(q, A_all)
    diff = U_pred - U_target
    return np.sum(diff**2)
def gradient(q, A_all, U_target):
    grad = np.zeros_like(q)

    for t, A in enumerate(A_all):
        U = 0.5 * q @ A @ q
        diff = U - U_target[t]

        grad += 2 * diff * (A @ q)

    return grad

def solve_charges(positions, U_target):

    A_all = compute_A_all(positions)

    N = positions.shape[1]

    q0 = np.random.randn(N)

    result = minimize(
        fun=objective,
        x0=q0,
        args=(A_all, U_target),
        jac=gradient,
        method='L-BFGS-B'
    )

    return result.x, result

if verbose:
    print("Starting charge optimization...")
charges, result = solve_charges(
    positions=charged_nodes_positions_log[:len(required_charge_energy_values)],
    U_target=required_charge_energy_values
)
if verbose:
    print("Charge optimization completed.")


def write_charges_to_csv(charges: list[float], filename="optimized_charges.csv"):
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Charge"])
        for charge in charges:
            writer.writerow([charge])

if verbose:
    print("Writing optimized charges to CSV: {}".format(output_file))
write_charges_to_csv(charges, filename=output_file)
if verbose:
    print("Optimized charges written to CSV successfully.")