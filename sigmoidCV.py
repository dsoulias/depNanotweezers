import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42 # Nature's guidance for python graphs
matplotlib.rcParams['font.family'] = 'Arial'

# Loading data file

data = pd.read_excel(input("CSV filepath: "), usecols=(7, 8, 9)) # delete quotes when pasting

cyclesRange = range(1, 1 + data["cycle number"].iloc[-1])
splitData = {}

for cycle in cyclesRange:
    splitData["data{0}".format(cycle)] = data.loc[data["cycle number"] == cycle]
    splitData["data{0}".format(cycle)].reset_index(drop=True, inplace=True) # resets index of each [dfs] Dataframe
    
ivR = {}
ivO = {}

for key in splitData.keys():
    splitIndex = splitData[key].index[splitData[key]["Ewe/V"] == min(splitData[key]["Ewe/V"])][0]
    ivR[key] = splitData[key][:splitIndex]
    ivO[key] = splitData[key][(splitIndex + 1):-1]
    ivO[key].reset_index(drop=True, inplace=True)

peakCurrentList = []

for key in splitData.keys():
    if key != "data1":
        basePointR1 = np.abs(ivR[key]["Ewe/V"] + 0.010).argmin()
        basePointR2 = np.abs(ivR[key]["Ewe/V"] + 0.030).argmin()
        baseVoltagesR = [ivR[key]["Ewe/V"][basePointR1], ivR[key]["Ewe/V"][basePointR2]]
        baseCurrentsR = [ivR[key]["<I>/mA"][basePointR1], ivR[key]["<I>/mA"][basePointR2]]
        pBaseR = np.polyfit(baseVoltagesR, baseCurrentsR, 1)
        platPointR1 = np.abs(ivR[key]["Ewe/V"] + 0.300).argmin()
        platPointR2 = np.abs(ivR[key]["Ewe/V"] + 0.350).argmin()
        platVoltagesR = [ivR[key]["Ewe/V"][platPointR1], ivR[key]["Ewe/V"][platPointR2]]
        platCurrentsR = [ivR[key]["<I>/mA"][platPointR1], ivR[key]["<I>/mA"][platPointR2]]
        pPlatR = np.polyfit(platVoltagesR, platCurrentsR, 1)
        midPointR = (ivR[key]["Ewe/V"].iloc[-1] - ivR[key]["Ewe/V"][0]) / 2
        peakCurrentR = (pPlatR[0] * midPointR + pPlatR[1]) - (pBaseR[0] * midPointR + pBaseR[1])
        peakCurrentList.append(peakCurrentR)
        
for key in splitData.keys():
    basePointO1 = np.abs(ivO[key]["Ewe/V"] + 0.010).argmin()
    basePointO2 = np.abs(ivO[key]["Ewe/V"] + 0.030).argmin()
    baseVoltagesO = [ivO[key]["Ewe/V"][basePointO1], ivO[key]["Ewe/V"][basePointO2]]
    baseCurrentsO = [ivO[key]["<I>/mA"][basePointO1], ivO[key]["<I>/mA"][basePointO2]]
    pBaseO = np.polyfit(baseVoltagesO, baseCurrentsO, 1)
    platPointO1 = np.abs(ivO[key]["Ewe/V"] + 0.300).argmin()
    platPointO2 = np.abs(ivO[key]["Ewe/V"] + 0.350).argmin()
    platVoltagesO = [ivO[key]["Ewe/V"][platPointO1], ivO[key]["Ewe/V"][platPointO2]]
    platCurrentsO = [ivO[key]["<I>/mA"][platPointO1], ivO[key]["<I>/mA"][platPointO2]]
    pPlatO = np.polyfit(platVoltagesO, platCurrentsO, 1)
    midPointO = (ivO[key]["Ewe/V"][0] - ivO[key]["Ewe/V"].iloc[-1]) / 2
    peakCurrentO = (pPlatO[0] * midPointO + pPlatO[1]) - (pBaseO[0] * midPointO + pBaseO[1])
    peakCurrentList.append(peakCurrentO)

meanCurrent = np.mean(peakCurrentList)
stdCurrent = np.std(peakCurrentList)

print("The mean steady-state current is =", round(meanCurrent*1e6, 1), u"\u00B1", round(stdCurrent*1e6, 1), "nA")

F_const = 96485.3312 # [C/mol], Faraday constant
cRedox = 10 # [mM = mol/m**3] concentration of Ruhex in 0.1 M KCl
DifCoeff = 9.1*1e-10 # [m^2 / s], [Ref. CRC Handbook 2017]
electrodeRadius = 1e6 * (abs(meanCurrent) * 1e-3 / (4 * F_const * DifCoeff * cRedox))
electrodeRadiusStd = 1e6 * (abs(stdCurrent) * 1e-3 / (4 * F_const * DifCoeff * cRedox))
print("Electrode radius =", round(electrodeRadius, 2), u"\u00B1", round(electrodeRadiusStd, 2), "um")

plt.figure(figsize=(6.68, 5.65))
plt.plot(ivR[key]["Ewe/V"]*1e3, ivR[key]["<I>/mA"]*1e6)
plt.plot(ivR[key]["Ewe/V"]*1e3, 1e6*(pBaseR[0] * ivR[key]["Ewe/V"] + pBaseR[1]), '--', color='red')
plt.plot(ivR[key]["Ewe/V"]*1e3, 1e6*(pPlatR[0] * ivR[key]["Ewe/V"] + pPlatR[1]), '--', color='red')
plt.axvline(x = midPointR*1e3, color = 'black', linestyle = ':')
plt.xlabel(r"$\rm E_{we}\ (mV)$", fontsize=18)
plt.ylabel("Current (nA)", fontsize=18)
plt.xticks(np.arange(-400, 100, 100), fontsize=16) 
plt.yticks(fontsize=16)
plt.tick_params(axis="x", direction="in", length=6, top=True)
plt.tick_params(axis="y", direction="in", length=6, right=True)
plt.tight_layout()
plt.show()