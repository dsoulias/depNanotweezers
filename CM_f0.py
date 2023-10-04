import numpy as np
import matplotlib.pyplot as plt

# Parameters

F = 96485.332 # [C/mol]; Faraday constant
R = 8.314 # [J/(K*mol)]; gas constant
k_B = 1.380649e-23 # [J/K]; Boltzmann's constant
q = 1.6e-19 # [C]; elementary charge
T = 25 + 273.15 # [K]; Room temperature
r_p = 1e-6 # [m]; radius of particle (latex polystyrene bead)
e_0 = 8.8541878e-12 # [F/m]; vacuum permittivity
e_m = 79 # [1]; relative permittivity of suspending medium (PBS pH 7.4)
e_p = 2.56 # [1]; relative permittivity of particle (latex polystyrene bead)
R_m = 18.2e4 # [Ohm*m]; resistivity of suspending medium (double-distilled water)
s_m = 1/R_m # [S/m]; conductivity of suspending medium (double-distilled water)
# s_p = 6e-3 # [S/m]; conductivity of particle (latex polystyrene bead)
visPBS = 0.99e-3 # [Pa*s=N*s/m^2]; viscosity of 10x PBS (pH 7.4) {CRC Handbook for 30 mM}
visH2O = 0.888e-3 # viscosity of water
c_b = 34 # [mol/m^3]; suspending medium concentration (PBS pH 7.4)
D_Na = 1.334e-9 # [m^2/s]; diffusion coefficient of Na+
mob_Na = D_Na * (F / (R * T))
z = -2.4e-3 # [V]; mean zeta potential
dz = 0.5e-3 # [V]; zeta potential standard error of the mean

# Double layer length

DebyeLength = np.sqrt((e_0 * e_m * R * T) / (2 * c_b * F**2)) # [m]
k = 1 / DebyeLength # [m^-1]

# Surface conductance of Stern inner layer

K_s = - e_0 * e_m * z * mob_Na * ((1 + k * r_p) / r_p)
K_s_er = e_0 * e_m * dz * mob_Na * ((1 + k * r_p) / r_p)

# Dukhin number

m = ((R * T / F)**2) * 2 * e_0 * e_m / (3 * visPBS * D_Na)

# Surface conductance of the diffusion layer

K_d = np.cosh((q * z / (2 * k_B * T)) - 1) * 4 * (F**2) * c_b * D_Na * (1 + (3 * m)) / (R * T * k)
K_d_er = -dz*(q / (2 * k_B * T) )* (4 * (F**2) * c_b * D_Na * (1 + (3 * m)) / (R * T * k)) * np.sinh(1 - (q * z / (2 * k_B * T)))

# Particle conductivity based on the contrivution from the diffusion layer

K = K_s + K_d
K_er = np.sqrt(K_s_er**2 + K_d_er**2)
s_p = (2 * K / r_p)
s_p_er = 2 * K_er / r_p

# Calculation of crossover frequency (Re[CM] = 0)

f_CO = (1 / (2 * np.pi *e_0)) * np.sqrt(((s_m - s_p) * (s_p + (2 * s_m))) / ((e_p - e_m) * (e_p + (2 * e_m))))
f_CO_er = -s_p_er * (1 / (2 * np.pi *e_0)) * (s_m + 2 * s_p) / ((2 * (e_p - e_m) * (e_p + (2 * e_m))) * np.sqrt((2 * (s_m **2) - (s_m * s_p) - s_p**2)/(2 * (e_p - e_m) * (e_p + (2 * e_m)))))

print (round(f_CO*1e-6,2), "MHz")

#%%
# Plotting Re[CM] over applied frequency

f = np.arange(0,1e9,step=1000) # [Hz]
f_rad = 2*np.pi*f # [Hz]
t = (e_p + (2*e_m))/(s_p+(2*s_m))

Re_CM = (e_0**2*f_rad**2*(e_p-e_m)*(e_p+2*e_m)+(s_p+2*s_m)*(s_p-s_m))/(e_0**2*f_rad**2*(e_p+2*e_m)**2+2*(s_p+2*s_m)**2)

plt.figure(figsize=(4.45, 2.95))
plt.plot(f,Re_CM, linewidth=3)
plt.xlabel('Frequency [Hz]', fontname="Arial", fontsize=14)
plt.ylabel('Re[f$_{CM}$]', fontname="Arial", fontsize=14)
plt.yticks(np.arange(-0.5, 0.6, 0.25), fontname="Arial", fontsize=11)
plt.ylim([-0.6,0.6])
plt.xscale("log")
plt.xticks(fontname="Arial", fontsize=10)
plt.tick_params(axis="x", which='major', direction="in", length=6, top=True)
plt.tick_params(axis="x", which='minor', direction="in", length=3, top=True)
plt.tick_params(axis="y", which='major', direction="in", length=6, right=True)
plt.tick_params(axis="y", which='minor', direction="in", length=3, right=True)
plt.axvline(x=980e3, ymin=0, ymax=0.5, color='black', linestyle='--', linewidth=2)
plt.tight_layout()
plt.show()
plt.savefig(r"C:\Users\User\Desktop\Re_CMvsf", dpi=300)