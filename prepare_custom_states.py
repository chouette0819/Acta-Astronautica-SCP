import sys
import numpy as np

def main():
    if len(sys.argv) < 3:
        print("usage: prepare_custom_states.py <npy_path> <out_csv>")
        sys.exit(1)
    npy_path = sys.argv[1]
    out_csv = sys.argv[2]
    arr = np.load(npy_path, allow_pickle=True)
    r = arr.reshape(arr.shape[0], -1)
    # Expect RTN positions in km: columns [R, T, N]
    r0_km = r[0, :3]
    rf_km = r[-1, :3]
    # Convert to meters and build [x, vx, y, vy, z, vz]
    r0 = (r0_km * 1e3).astype(float)
    rf = (rf_km * 1e3).astype(float)
    x0 = np.array([r0[0], 0.0, r0[1], 0.0, r0[2], 0.0], dtype=float)
    xf = np.array([rf[0], 0.0, rf[1], 0.0, rf[2], 0.0], dtype=float)
    data = np.vstack([x0, xf])
    np.savetxt(out_csv, data, delimiter=",", fmt="%.9f")
    # Also emit suggested N based on samples
    print(f"N_suggest={r.shape[0]}")

if __name__ == "__main__":
    main()

