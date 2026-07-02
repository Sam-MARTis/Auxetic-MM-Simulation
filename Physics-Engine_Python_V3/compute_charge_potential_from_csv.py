import csv
import numpy as np
import scipy.sparse as sp
from Constants import CHARGE_CONSTANT

def read_first_charge_positions(filename='charge_distance_data.csv'):
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        row = next(reader)
    # find columns for charge_*_x and charge_*_y
    xs = []
    ys = []
    for i, col in enumerate(header):
        if col.endswith('_x'):
            xs.append(float(row[i]))
    for i, col in enumerate(header):
        if col.endswith('_y'):
            ys.append(float(row[i]))
    positions = np.vstack([xs, ys]).T
    return positions

def read_optimized_charges(filename='optimized_charges.csv'):
    qs = []
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        try:
            next(reader)
        except StopIteration:
            return np.array([], dtype=float)
        for row in reader:
            if not row: continue
            qs.append(float(row[0]))
    return np.array(qs, dtype=float) * CHARGE_CONSTANT

def compute_energy_gradient_hessian(positions: np.ndarray, q: np.ndarray):
    m = len(q)
    assert positions.shape == (m, 2)
    # all pairs i<j
    idx_i, idx_j = np.triu_indices(m, k=1)

    dx = positions[idx_i, 0] - positions[idx_j, 0]
    dy = positions[idx_i, 1] - positions[idx_j, 1]
    d = np.sqrt(dx*dx + dy*dy)

    live = d > 1e-12
    idx_i, idx_j, dx, dy, d = idx_i[live], idx_j[live], dx[live], dy[live], d[live]
    qprod = q[idx_i] * q[idx_j]

    # Energy
    energy = 0.5 * np.sum(qprod / d)

    # Gradient (2m vector): contributions to i and j
    inv = 1.0 / d
    inv2 = inv * inv
    inv3 = inv2 * inv
    ex = dx * inv
    ey = dy * inv

    gx = - qprod * inv2 * ex
    gy = - qprod * inv2 * ey

    g = np.zeros(2*m, dtype=float)
    np.add.at(g, 2*idx_i,   gx);  np.add.at(g, 2*idx_i+1, gy)
    np.add.at(g, 2*idx_j,  -gx);  np.add.at(g, 2*idx_j+1,-gy)

    # Hessian as sparse matrix
    scale = qprod * inv3
    c00 = scale * (1.0 - 3.0 * ex * ex)
    c01 = scale * (    - 3.0 * ex * ey)
    c11 = scale * (1.0 - 3.0 * ey * ey)

    r1x, r1y = 2*idx_i,   2*idx_i+1
    r2x, r2y = 2*idx_j,   2*idx_j+1

    rows = np.concatenate([r1x,r1x,r1y,r1y, r2x,r2x,r2y,r2y,
                           r1x,r1x,r1y,r1y, r2x,r2x,r2y,r2y])
    cols = np.concatenate([r1x,r1y,r1x,r1y, r2x,r2y,r2x,r2y,
                           r2x,r2y,r2x,r2y, r1x,r1y,r1x,r1y])
    vals = np.concatenate([ c00, c01, c01, c11,
                            c00, c01, c01, c11,
                           -c00,-c01,-c01,-c11,
                           -c00,-c01,-c01,-c11])

    H = sp.csc_matrix((vals, (rows, cols)), shape=(2*m, 2*m))

    return energy, g, H

def main():
    positions = read_first_charge_positions()
    q = read_optimized_charges()
    if len(q) == 0:
        raise RuntimeError('No optimized charges found')
    if positions.shape[0] != len(q):
        raise RuntimeError(f'Mismatch: {positions.shape[0]} positions vs {len(q)} charges')

    energy, g, H = compute_energy_gradient_hessian(positions, q)

    print(f'Charge energy = {energy:.12e}')
    np.save('charge_gradient.npy', g)
    sp.save_npz('charge_hessian.npz', H)
    with open('charge_energy.txt', 'w') as f:
        f.write(str(energy) + '\n')
    print('Saved: charge_energy.txt, charge_gradient.npy, charge_hessian.npz')

if __name__ == "__main__":
    main()
