import numpy as np
import matplotlib.pyplot as plt

# Parameters

R = 1e-6 #[m], radius of particle
s = np.arange(1e-9, 1e-4, 5e-10) # [m], distance between glass wall and particle circumference
h = s + R # [m], distance between glass wall and particle's COM

#%% Faxen's correction factors

# component parallel to xy plane (valid for s >> 0)

parFaxen = (1 - ((9/16)*(R/h)) + ((1/8)*(R/h)**3) - ((45/256)*(R/h)**4) - ((1/16)*(R/h)**5))**(-1)

# component perpendicular to xy plane

perFaxen = (1 - ((9/8)*(R/h)) + ((1/2)*(R/h)**3))**(-1)

#%% Brenner's correction factors

# component parallel to xy plane (valid for s ~= 0)

parBrenner = 0.9588 - ((8/15) * np.log(s/R))

#%% Plotting

plt.figure(figsize=(2.95, 2.95))
plt.plot(s, parFaxen, '--', color='black', label=r"$\rm \gamma_{//}^F$")
plt.plot(s, parBrenner, '-.', color='grey', label=r"$\rm \gamma_{//}^B$")
plt.axhline(y = 1, color = 'red', linestyle = ':', label="Bulk")
plt.xscale("log")
plt.ylim(bottom=0.0)
plt.xlabel("Distance s (m)", fontname="Arial" , fontsize=14)
plt.ylabel(r" Correction factor $\rm \gamma_{//}$", fontname="Arial" , fontsize=14)
plt.xticks(fontname="Arial", fontsize=12) 
plt.yticks(fontname="Arial", fontsize=12)
plt.tick_params(axis="x", which='major', direction="in", length=6, top=True)
plt.tick_params(axis="x", which='minor', direction="in", length=3, top=True)
plt.tick_params(axis="y", which='major', direction="in", length=6, right=True)
plt.tick_params(axis="y", which='minor', direction="in", length=3, right=True) 
plt.legend(prop={'size': 8, 'family': 'Arial'}, title="Component", title_fontproperties={'family': 'Arial', 'size': 10, 'weight':'bold'}, loc='upper right')
plt.tight_layout()
plt.show()

instr = input("Save figure? (yes/no): ")
if instr == "yes":
    figFolder = input("Folder path to save figure: ")
    plt.savefig(figFolder + "/FaxenBrenner_300.png", dpi=300)
    print("Figure saved!")
else:
    print("Figure was not saved!")

#%% Distance at which Faxen and Brenner parallel components are equal

difList = [abs(parFaxen - parBrenner)]
difMin = min(difList[0])
eqPoint = np.where(difList[0] == difMin)
print("The parallel components are equal for a correction factor =", round(float(parFaxen[eqPoint[0]]), 3))
print("at s =", float(s[eqPoint[0]]))