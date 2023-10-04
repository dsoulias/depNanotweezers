import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42 # Nature's guidance for python graphs
matplotlib.rcParams['font.family'] = 'Arial'

difCoef_220808 = [0.18692875, 0.18853941, 0.20073641, 0.18174713, 0.16913403,
                  0.21302487, 0.1929759, 0.18791494, 0.14651821, 0.2105694,
                  0.16482973, 0.09595498, 0.18597475, 0.06882037, 0.1662355,
                  0.1571227, 0.1363569, 0.16199517, 0.1366351, 0.10553301]

difCoef_221010 = [0.18883896, 0.10694386, 0.168782, 0.12989387, 0.13757786,
                  0.19378314, 0.15446586, 0.21585646, 0.16292237, 0.1614728, 
                  0.14803347, 0.18346576, 0.30263322, 0.22091672, 0.12234933,
                  0.31867698, 0.16340911, 0.21693684, 0.19626537, 0.14125344]

difCoef_221011 = [0.20600543, 0.16113976, 0.1546856, 0.1811758, 0.17295129,
                  0.15624694, 0.20921312, 0.17411476, 0.17864975, 0.20344049,
                  0.21578583, 0.14527122, 0.1839055, 0.13855913, 0.20808559,
                  0.166775, 0.14895506, 0.19452901, 0.2795419, 0.17753803]

difCoef_221012 = [0.10396726, 0.16611844, 0.15381223, 0.13490262, 0.19232924,
                  0.219087, 0.15990424, 0.14816554, 0.15458216, 0.14262635,
                  0.17823451, 0.15354489, 0.15255916, 0.20971119, 0.16102918,
                  0.17884777, 0.14616543, 0.13414079, 0.1571568, 0.05121192]

difCoef_221013 = [0.16576054, 0.15574929, 0.20949764, 0.1813968, 0.15103316,
                  0.19827558, 0.13395518, 0.21312437, 0.14824534, 0.14316445,
                  0.20327308, 0.17372289, 0.17353279, 0.17479939, 0.16371922,
                  0.20904334, 0.14676926, 0.18380771, 0.20070533, 0.16222564]

difCoef_221014 = [0.1512116, 0.17545919, 0.18009803, 0.19801951, 0.15612859,
                  0.14068639, 0.19761084, 0.12423937, 0.14609144, 0.12861151,
                  0.13145877, 0.18225295, 0.20183657, 0.17651817, 0.11746556,
                  0.21122913, 0.17736272, 0.34477545, 0.16725828, 0.21158472]

difCoefDF = pd.DataFrame({'1': difCoef_220808,
                          '2': difCoef_221010,
                          '3': difCoef_221011,
                          '4': difCoef_221012,
                          '5': difCoef_221013,
                          '6': difCoef_221014})

fig, ax = plt.subplots(figsize=(2.95, 2.95))
difCoefDF.plot.kde(ax=ax)
plt.xlabel(r"$\rm D_T\ (\mu m^2/s)$", fontsize=14)
plt.ylabel("PDF", fontsize=14)
plt.xticks(fontsize=12) 
plt.yticks(fontsize=12)
plt.tick_params(axis="x", direction="in", length=6, top=True)
plt.tick_params(axis="y", direction="in", length=6, right=True)
plt.ylim(0)
plt.legend(prop={'size': 8}, title="Set", title_fontproperties={'size': 10, 'weight':'bold'}, loc='upper right')
plt.tight_layout()
plt.show()

instr = input("Save figure? (yes/no): ")
if instr == "yes":
    figFolder = input("Folder path to save figure: ")
    plt.savefig(figFolder + "/meanDiffusion_300.png", dpi=300)
    print("Figure saved!")
else:
    print("Figure was not saved!")

#%% Temperature influence at the beginning of illumination and 30 minutes after

avgFirstList = []
avgLastList = []
stdFirstList = []
stdLastList = []
measurement = []

for column in difCoefDF:
    avgFirst = np.mean(difCoefDF[column][:10])
    stdFirst = np.std(difCoefDF[column][:10])
    avgLast = np.mean(difCoefDF[column][10:-1])
    stdLast = np.std(difCoefDF[column][10:-1])
    avgFirstList.append(avgFirst)
    avgLastList.append(avgLast)
    stdFirstList.append(stdFirst)
    stdLastList.append(stdLast)
    measurement.append(column)
    
difCoefTempDF = pd.DataFrame({"Set": measurement, "Mean first": avgFirstList,
                              "Std first": stdFirstList, "Mean last": avgLastList,
                              "Std last": stdLastList})

#%% KDE plot of calculated medium viscosity

R = 1e-6 # particle's radius, [m]
k_B = 1.3806485279e-23 # Boltzmann constant, [J/K]
T = 273 + 25 # room temperature, [K]
visc0 = 0.88e-3 # water viscosity (25 degC), [Pa*s] (measured)

columns = list(difCoefDF)
fullViscList = []
viscDict= {}

for clm in columns:
    viscList = []
    for D_T in difCoefDF[clm]:
        visc = ((k_B * T) / (6 * np.pi * D_T * 1e-12 * R)) / visc0 # Stokes-Einstein equation
        viscList.append(visc)
    fullViscList.extend(viscList)
    viscDict["{0}".format(clm)] = viscList

viscDF =  pd.DataFrame(viscDict)
viscListDF = pd.DataFrame(fullViscList)
meanRatio = np.mean(fullViscList)
stdRatio = np.std(fullViscList)
print(round(meanRatio, 2), round(stdRatio, 2))

fig, ax = plt.subplots(figsize=(2.95, 2.95))
viscDF.plot.kde(ax=ax)
plt.axvline(x = visc0 / visc0, color = 'gray', linestyle = '--', label = "Bulk")
plt.xlabel(r"$\rm \gamma_{//}$", fontsize=14)
plt.ylabel("PDF", fontname="Arial" , fontsize=14)
plt.xticks(fontsize=12) 
plt.yticks(fontsize=12)
plt.tick_params(axis="x", direction="in", length=6, top=True)
plt.tick_params(axis="y", direction="in", length=6, right=True)
plt.ylim(0)
plt.legend(prop={'size': 8}, title="Set", title_fontproperties={'size': 10, 'weight':'bold'}, loc='upper right')
plt.tight_layout()
plt.show()

instr = input("Save figure? (yes/no): ")
if instr == "yes":
    figFolder = input("Folder path to save figure: ")
    plt.savefig(figFolder + "/correctionFactor_300.png", dpi=300)
    print("Figure saved!")
else:
    print("Figure was not saved!")
    
fig, ax = plt.subplots(figsize=(2.95, 2.95))
viscListDF.plot.kde(ax=ax, legend=False)
plt.axvline(x = visc0 / visc0, color = 'gray', linestyle = '--', label = "Bulk")
plt.xlabel(r"$\rm \gamma_{//}$", fontsize=14)
plt.ylabel("PDF", fontsize=14)
plt.xticks(fontsize=12) 
plt.yticks(fontsize=12)
plt.tick_params(axis="x", direction="in", length=6, top=True)
plt.tick_params(axis="y", direction="in", length=6, right=True)
plt.ylim(0)
plt.tight_layout()
plt.show()

instr = input("Save figure? (yes/no): ")
if instr == "yes":
    figFolder = input("Folder path to save figure: ")
    plt.savefig(figFolder + "/meanCorrectionFactor_300.png", dpi=300)
    print("Figure saved!")
else:
    print("Figure was not saved!")
        
#%% KDE plot of calculated distance between glass wall and particle circumference

eqPoint = 3.45e-08 # distance at which the parallel components are equal
eqFactor = 2.755 # factor for which the parallel components are equal

distanceDict ={}
fullList = []

for key in viscDict.keys():
    distanceList = []
    for visc in viscDict[key]:
        if visc >= eqFactor:
            s = R * np.exp((15/8) * (0.9588 - visc))
        else:
            p = np.poly1d([1/16, 45/256, -1/8, 0, 9/16, (1/visc) - 1])
            pSol = p.roots
            for w in pSol:
                if np.isreal(w) == True:
                    if w > 0:    
                        s = (R / w.real) - R
        distanceList.append(s*1e6)
    distanceDict["{0}".format(key)] = distanceList
    fullList.extend(distanceList)

distanceDF = pd.DataFrame(distanceDict)
fullDF = pd.DataFrame(fullList)
meanDistance = np.mean(fullList)
stdDistance = np.std(fullList)
print(round(meanDistance, 2), round(stdDistance, 2))

fig, ax = plt.subplots(figsize=(2.95, 2.95))
distanceDF.plot.kde(ax=ax)
plt.xlabel(r"Distance s ($\rm \mu m$)", fontsize=14)
plt.ylabel("PDF", fontsize=14)
plt.xticks(fontname="Arial", fontsize=12) 
plt.yticks(fontname="Arial", fontsize=12)
plt.tick_params(axis="x", direction="in", length=6, top=True)
plt.tick_params(axis="y", direction="in", length=6, right=True)
plt.ylim(0)
plt.legend(prop={'size': 8}, title="Set", title_fontproperties={'size': 10, 'weight':'bold'}, loc='upper right')
plt.tight_layout()
plt.show()

fig, ax = plt.subplots(figsize=(2.95, 2.95))
fullDF.plot.kde(ax=ax, legend=False)
plt.xlabel(r"Distance s ($\rm \mu m$)", fontsize=14)
plt.ylabel("PDF", fontsize=14)
plt.xticks(fontsize=12) 
plt.yticks(fontsize=12)
plt.tick_params(axis="x", direction="in", length=6, top=True)
plt.tick_params(axis="y", direction="in", length=6, right=True)
plt.ylim(0)
plt.tight_layout()
plt.show()