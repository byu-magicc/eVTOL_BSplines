import numpy as np

# Set seed for reproducibility
np.random.seed(0)

# Create known matrices
A = np.random.rand(5, 3)       # 5x3
Z_true = np.random.rand(3, 4)  # 3x4
B = A @ Z_true                 # B is 5x4

# Solve for Z in AZ = B using lstsq
# Each column: A @ Z[:, i] = B[:, i]

Z_est, residuals, rank, s = np.linalg.lstsq(A, B, rcond=None)

# Compare result
print("True Z:\n", Z_true)
print("\nEstimated Z:\n", Z_est)
print("\nDifference (Z_true - Z_est):\n", Z_true - Z_est)


potato = 0