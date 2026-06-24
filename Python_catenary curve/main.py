#
# Calculation of the shape of a catenary curve and the tensile stresses developed in it.
#
# Dmitriy Makhnovskiy, 24.06.2026
#

import numpy as np
import sys
import matplotlib.pyplot as plt

# Enter the initial parameters
vd = 2710  # wire volume density (kg/m^3)
S = 1.5  # wire cross-section (mm^2)
p1 = (0.0, 20.0)  # coordinates (m, m) of the left fixed end
p2 = (100.0, 30.0)  # coordinates (m, m) of the right fixed end

# Newton iteration parameters (recommended)
s0 = 1.0
tol=1e-12
max_iter=100

# Interactive length and recalculated S (do not change)
g = 9.81  # free fall acceleration (m / s^2)
S = S * 1.0e-6  # wire cross-section (m^2)
x1 = p1[0]
y1 = p1[1]
x2 = p2[0]
y2 = p2[1]

if x2 <= x1:
    print("Error: the coordinate x2 must be larger than x1.")
    print("Program stopped.")
    sys.exit()

Lmin = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)  # minimum wire length

while True:
    print(f"Minimum wire length Lmin = {Lmin:.6g} m")
    try:
        L = float(input("Choose a wire length (m), larger than Lmin: "))
    except ValueError:
        print("Invalid input. Please enter a number.")
        continue
    if L > Lmin:
        break
    print("Wire length must be larger than Lmin.")

rho = vd * S  # linear mass density (kg/m)

def f(s):
    C = np.sqrt(L**2 - (y1 - y2)**2) / (x2 - x1)
    if abs(s) < 1e-8:
        # limit of (exp(s)-1)/(s*exp(s/2)) as s -> 0
        phi = 1 + s**2 / 24
    else:
        phi = (np.exp(s) - 1) / (s * np.exp(s / 2))
    return phi - C

def df(s):
    if abs(s) < 1e-8:
        # derivative near s = 0
        return s / 12
    else:
        numerator = (
            (s / 2) * np.exp(3 * s / 2)
            - np.exp(3 * s / 2)
            + np.exp(s / 2)
            + (s / 2) * np.exp(s / 2)
        )
        denominator = s**2 * np.exp(s)
        return numerator / denominator

def newton_s():
    s = s0
    for _ in range(max_iter):
        fs = f(s)
        dfs = df(s)
        s_new = s - fs / dfs
        error = abs(s_new - s)
        if error < tol:
            return s_new
        s = s_new
    raise RuntimeError("Newton method did not converge")

s = newton_s()
dy = y1 - y2
r = 0.5 * np.log(np.exp(s) * (L - dy) / (L + dy))
q = r - s
a = (x2 - x1) / s
b = x1 - a * q
lam = (a * rho * g / 2) * (np.exp((x1 - b) / a) + np.exp(-(x1 - b) / a)) - y1 * rho * g

print('')
print('Parameters of the catenary curve:')
print('a = ', a)
print('b = ', b)
print('lambda = ', lam)
print('rho = ', rho, ' kg / m')

print('')
print('Tension in the catenary curve:')
u1 = (x1 - b) / a
u2 = (x2 - b) / a

T1 = g * rho * a * np.sqrt(1 + 0.25 * (np.exp(u1) - np.exp(-u1))**2)
T2 = g * rho * a * np.sqrt(1 + 0.25 * (np.exp(u2) - np.exp(-u2))**2)
sigma1 = (g * rho * a / S) * np.sqrt(1 + 0.25 * (np.exp(u1) - np.exp(-u1))**2)
sigma2 = (g * rho * a / S) * np.sqrt(1 + 0.25 * (np.exp(u2) - np.exp(-u2))**2)


if x1 <= b <= x2:
    Tmin = g * rho * a
    sigma_min = g * rho * a / S
else:
    Tmin = min(T1, T2)
    sigma_min = min(sigma1, sigma2)

Tmax = max(T1, T2)
sigma_max = max(sigma1, sigma2)

print("T(x1) =", T1, "N")
print("T(x2) =", T2, "N")
print("Tmin =", Tmin, "N")
print("Tmax =", Tmax, "N")
print('')
print(f"sigma(x1) = {sigma1:.6e} Pa")
print(f"sigma(x2) = {sigma2:.6e} Pa")
print(f"sigma_min = {sigma_min:.6e} Pa")
print(f"sigma_max = {sigma_max:.6e} Pa")

# Plotting the calculated catenary curve and checking the fixed end points
# Coordinate array
x = np.linspace(x1, x2, 500)
u = (x - b) / a

# Catenary curve y(x)
y = (a / 2.0) * (np.exp(u) + np.exp(-u)) - lam / (rho * g)

# Tension T(x)
T = g * rho * a * np.sqrt(
    1.0 + 0.25 * (np.exp(u) - np.exp(-u))**2
)

# Tensile stress sigma(x)
sigma = T / S

# First figure: catenary curve and fixed end points
plt.figure(figsize=(9, 5))

plt.plot(x, y, linewidth=2, label="Calculated catenary curve")
plt.scatter([x1, x2], [y1, y2], s=120, linewidths=2.5, edgecolors="black", label="Fixed end points")

plt.text(x1, y1, "  p1", fontsize=12, fontweight="bold")
plt.text(x2, y2, "  p2", fontsize=12, fontweight="bold")

plt.xlabel("x, m")
plt.ylabel("y, m")
plt.title("Calculated catenary curve")
plt.grid(True)
plt.axis("equal")
plt.legend()
plt.show()

# Second figure: T(x) and sigma(x)
fig, ax1 = plt.subplots(figsize=(9, 5))

ax1.plot(x, T, linewidth=2)
ax1.set_xlabel("x, m")
ax1.set_ylabel("Tension T(x), N")
ax1.grid(True)

ax2 = ax1.twinx()
ax2.plot(x, sigma, linewidth=2, linestyle="--")
ax2.set_ylabel("Tensile stress sigma(x), Pa")
ax2.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

plt.show()
