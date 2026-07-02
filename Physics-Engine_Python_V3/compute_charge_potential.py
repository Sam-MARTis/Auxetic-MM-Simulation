import pickle
import csv
import numpy as np
import scipy.sparse as sp
from scipy import sparse
from Objects import Solid
from Constants import CHARGE_CONSTANT, OPTIMIZED_CHARGE_STORAGE_FILE, OBJECT_SAVE_FILE

def load_solid(filename: str) -> Solid:
    with open(filename, 'rb') as f:
        solid = pickle.load(f)
    return solid

def load_optimized_charges(filename: str) -> np.ndarray:
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        charges = [float(row[0]) for row in reader]
    return np.array(charges, dtype=np.float64)

def main():
    solid = load_solid(OBJECT_SAVE_FILE)

    # load optimized charges if present and assign (scaled by CHARGE_CONSTANT)
    try:
        opt_q = load_optimized_charges(OPTIMIZED_CHARGE_STORAGE_FILE)
        solid.charges = opt_q * CHARGE_CONSTANT
        print(f"Loaded {len(opt_q)} optimized charges and scaled by CHARGE_CONSTANT={CHARGE_CONSTANT}")
    except Exception as e:
        print(f"Could not load optimized charges ({e}), using existing solid.charges")

    # Compute pairwise A matrix for charged positions and the scalar energy
    charged_pos = solid.node_positions.reshape(solid.num_nodes, 2)[solid.charged_indices]
    q = solid.charges

    A = solid.compute_A(charged_pos)
    energy = 0.5 * q @ A @ q

    print(f"Charge energy (0.5 q^T A q) = {energy:.12e}")

    # Compute gradient and hessian from the Solid routine
    g, H = solid.assemble_electrostatic_gradient_and_hessian()

    # Save outputs
    np.save('charge_gradient.npy', g)
    sp.save_npz('charge_hessian.npz', H)
    with open('charge_energy.txt', 'w') as f:
        f.write(str(energy) + '\n')

    print('Saved: charge_energy.txt, charge_gradient.npy, charge_hessian.npz')

if __name__ == '__main__':
    main()
