import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42 # Nature's guidance for python graphs
matplotlib.rcParams['font.family'] = 'Arial'

# # Load XML data file...!delete quotes when pasting the filepath in Console!

# data = pd.read_csv(input("XML filepath: "), sep="\t", header=None) # for tab separated file without any headers
# data.columns = ["t", "x", "y"] # label each file column to represented value

# Load CSV data file...

data = pd.read_csv(input("CSV filepath: "), usecols=(4, 5, 7))
data.columns = ["x", "y", "t"] # label each file column to represented value

#%% Trajectory splitting

# Define function to split dataframes based on desired position (Ref. https://datagy.io/split-pandas-dataframe/ )

def split_dataframe_by_position(df, splits):
    """
    Takes a dataframe and an integer of the number of splits to create.
    Returns a list of dataframes.
    """
    dataframes = []
    index_to_split = len(df) // splits
    start = 0
    end = index_to_split
    for split in range(splits):
        temporary_df = df.iloc[start:end, :]
        dataframes.append(temporary_df)
        start += index_to_split
        end += index_to_split
    return dataframes

# Split trajectory into shorter segments

Nseg = [1, 2, 3, 4, 5, 10, 20, 40] # number of tracks per dataframe
splitData = {} # empty dictionary to store all splitted new data sets

for i in Nseg:
    splitData["tracks{0}".format(i)] = split_dataframe_by_position(data, i)
    
#%% (optional) Plotting full trajectory

plt.figure(figsize=(2.95, 2.95))
plt.plot(splitData["tracks1"][0]["x"], splitData["tracks1"][0]["y"], color='black')
plt.locator_params(axis='both', nbins=4)
plt.xlabel(r"x ($\rm \mu m$)", fontsize=14)
plt.ylabel(r"y ($\rm \mu m$)", fontsize=14)
plt.xticks(fontsize=12) 
plt.yticks(fontsize=12)
plt.tick_params(axis="x", which='major', direction="in", length=6, top=True)
plt.tick_params(axis="x", which='minor', direction="in", length=3, top=True)
plt.tick_params(axis="y", which='major', direction="in", length=6, right=True)
plt.tick_params(axis="y", which='minor', direction="in", length=3, right=True) 
plt.tight_layout()
plt.show()

instr = input("Save figure? (yes/no): ")
if instr == "yes":
    figFolder = input("Folder path to save figure: ")
    plt.savefig(figFolder + "/figure1a_300.png", dpi=300)
    print("Figure saved!")
else:
    print("Figure was not saved!")

#%% (optional) Plotting segmented trajectories

cmap = plt.get_cmap("viridis") # viridis is an example

fig, ax = plt.subplots(1, 1, figsize=[2.95, 2.95]) # define a new figure with 1 plot
plot = ax.scatter(splitData["tracks1"][0]["x"], splitData["tracks1"][0]["y"], c=splitData["tracks1"][0]["t"], cmap=cmap, vmax=15, s=4)
cbar = fig.colorbar(plot, ax=ax, orientation='vertical')
cbar.ax.set_title('t (s)', fontsize=14)
cbar.ax.tick_params(labelsize=12)
plt.locator_params(axis='both', nbins=4)
plt.xlabel(r"x ($\rm \mu m$)", fontsize=14)
plt.ylabel(r"y ($\rm \mu m$)", fontsize=14)
plt.xticks(fontsize=12) 
plt.yticks(fontsize=12)
plt.tick_params(axis="x", which='major', direction="in", length=6, top=True)
plt.tick_params(axis="x", which='minor', direction="in", length=3, top=True)
plt.tick_params(axis="y", which='major', direction="in", length=6, right=True)
plt.tick_params(axis="y", which='minor', direction="in", length=3, right=True)
plt.tight_layout()

instr = input("Save figure? (yes/no): ")
if instr == "yes":
    figFolder = input("Folder path to save figure: ")
    plt.savefig(figFolder + "/figure1b_300.png", dpi=300)
    print("Figure saved!")
else:
    print("Figure was not saved!")

#%% MSD calculation for each data set

values = [] # list where all values from spliData dictionary are stored
dfs = [] # list where all Dataframes from "values" list are stored

for value in splitData.values(): # stores each splitData value to [values]
    values.append(value)   
    for df in value: # stores every Dataframe from [values] to [dfs]
        dfs.append(df)
        
nArrays = [] # list where all possible displacements of each Dataframe is stored
        
for segment in dfs: # stores each array (n) to [nArrays]
    n = np.array(range(1, len(segment), 1)) # array of all possible displacements for each Dataframe
    nArrays.append(n)

msdData = [] # list where msd values for all possible displacements (n) of each Dataframe are stored
row = np.array(range(len(nArrays))) # array with all index values of [nArrays]

for number in row: # reads every value in row array
    dtList = [] # list where all delays are stored
    nList = [] # list where all displacements are stored
    msdList = [] # list where all msd values are stored
    dfs[number].reset_index(drop=True, inplace=True) # resets index of each [dfs] Dataframe
    
    for i in nArrays[number]: # reads every possible displacement of each [nArrays] row
        dt = dfs[number]["t"][i] - dfs[number]["t"][0] # delay for every possible displacement
        nList.append(i) # stores every displacement in [nList]
        dtList.append(dt) # stores every delay in [dtList]
        jRange = np.array(range(len(dfs[number]) - i)) # sum range for each Dataframe and displacement
        total = 0 # initial value of sum
        for j in jRange: # reads every value in the sum range
            total += (dfs[number]["x"][j + i] - dfs[number]["x"][j])**2 + (dfs[number]["y"][j + i] - dfs[number]["y"][j])**2
        msd = total / (len(dfs[number]) - i)
        msdList.append(msd)
        
    ndtDataframe = pd.DataFrame({"n": nList, "dt": dtList, "msd": msdList})
    msdData.append(ndtDataframe)
    
#%% (optional) Plotting MSD over lag time for specific Nseg

plt.figure(figsize=(2.95, 2.95))

plt.plot(msdData[0]["dt"], msdData[0]["msd"], color='black', label='1000')
plt.plot(msdData[1]["dt"], msdData[1]["msd"], color=u'#1f77b4', label='500')
plt.plot(msdData[2]["dt"], msdData[2]["msd"], color=u'#1f77b4')
plt.plot(msdData[3]["dt"], msdData[3]["msd"], color=u'#ff7f0e', label='333')
plt.plot(msdData[4]["dt"], msdData[4]["msd"], color=u'#ff7f0e')
plt.plot(msdData[5]["dt"], msdData[5]["msd"], color=u'#ff7f0e')
plt.xlabel(r"$\rm \tau$ (s)", fontsize=14)
plt.ylabel(r"MSD ($\rm \mu m^2 /s$)", fontsize=14)
plt.xticks(np.arange(0, 1 + max(msdData[0]["dt"]), step=4), fontsize=12) 
plt.yticks(fontsize=12)
plt.tick_params(axis="x", which='major', direction="in", length=6, top=True)
plt.tick_params(axis="x", which='minor', direction="in", length=3, top=True)
plt.tick_params(axis="y", which='major', direction="in", length=6, right=True)
plt.tick_params(axis="y", which='minor', direction="in", length=3, right=True)
plt.legend(prop={'size': 8}, title=r"N$\rm _{seg}$", title_fontproperties={'size': 12, 'weight':'bold'}, loc='lower right') 
plt.tight_layout()
plt.show()

instr = input("Save figure? (yes/no): ")
if instr == "yes":
    figFolder = input("Folder path to save figure: ")
    plt.savefig(figFolder + "/figure1c_300.png", dpi=300)
    print("Figure saved!")
else:
    print("Figure was not saved!")

#%% Fitting: optimum number of data points

skipRows = 0
slopesData = {}

for i in Nseg:
    tracksRange = np.array(range(skipRows, skipRows + i))
    skipRows = skipRows + i
    
    slopeDataframe = pd.DataFrame()
    
    for selRange in tracksRange:
        limitList = []
        slopeList = []
        interList = []
        limitRange = np.array(range(1, len(msdData[selRange]), 1))
        for limit in limitRange:
            coef = np.polyfit(msdData[selRange]["dt"][:limit], msdData[selRange]["msd"][:limit], 1)
            slopeList.append(coef[0])
            interList.append(coef[1])
            limitList.append(limit)
        slopeDataframe["n"] = limitList
        slopeDataframe["D*{0}".format(selRange)] = slopeList
    
    slopeDataframeRange = range(len(slopeDataframe))
    slopeMeanList = []
    slopeStdList = []
    histMeanList = []
    histStdList = []
    
    for g in slopeDataframeRange:
        # based on mean and std
        slopeMean = np.mean(slopeDataframe.iloc[g, 1:])
        slopeStd = np.std(slopeDataframe.iloc[g, 1:])
        slopeMeanList.append(slopeMean)
        slopeStdList.append(slopeStd)
        # # based on histogram mean and std
        # hist, bins = np.histogram(slopeDataframe.iloc[g, 1:])
        # mids = 0.5*(bins[1:] + bins[:-1])
        # histMean = np.average(mids, weights=hist)
        # histStd = np.sqrt(np.average((mids - histMean)**2, weights=hist))
        # histMeanList.append(histMean)
        # histStdList.append(histStd)
    
    slopeDataframe["D*_mean"] = slopeMeanList
    slopeDataframe["D*_std"] = slopeStdList
    # slopeDataframe["D*_histMean"] = histMeanList
    # slopeDataframe["D*_histStd"] = histStdList
    slopeDataframe["s/D*_mean"] = [m /n for m, n in zip(slopeStdList, slopeMeanList)]
    # slopeDataframe["sHist/D*_histMean"] = [m /n for m, n in zip(histStdList, histMeanList)]
    slopesData["slope{0}".format(i)] = slopeDataframe

#%% Plotting & finding optimum number of fitting data points

plt.figure(figsize=(4.45, 2.95))
yMinPoints = []
indexPoints = []

for slp in slopesData.keys():
    if np.min(slopesData[slp]["s/D*_mean"]) == 0:
        print("N/A for", slp)
    else:
        plt.loglog(slopesData[slp]["n"], slopesData[slp]["s/D*_mean"], 'x', ms=3, label=Nseg[list(slopesData.keys()).index(slp)])
        mPList = []
        for mP in slopesData[slp]["s/D*_mean"][1:(len(slopesData[slp]["s/D*_mean"]) // 4)]:
            if mP > 0.005:
                mPList.append(mP)
        yMinPoint = np.min(mPList)
        indexPoint = slopesData[slp][slopesData[slp]["s/D*_mean"] == yMinPoint].index[0]
        yMinPoints.append(yMinPoint)
        indexPoints.append(indexPoint)
        
optPoints = pd.DataFrame({"index": indexPoints, "s/D*_mean": yMinPoints})
optPoint = optPoints["index"][optPoints["index"][optPoints["s/D*_mean"] == min(optPoints["s/D*_mean"])].index[0]]
print("Optimum number of fitting points =", optPoint + 1, ", with relative error =", round(min(optPoints["s/D*_mean"]), 3))

plt.xlabel("Fitting points n", fontsize=14)
plt.ylabel(r"$\rm s_{D^*}/\ \overline{D^*}$", fontsize=14)
plt.xticks(fontsize=12) 
plt.yticks(fontsize=12)
plt.tick_params(axis="x", which='major', direction="in", length=6, top=True)
plt.tick_params(axis="x", which='minor', direction="in", length=3, top=True)
plt.tick_params(axis="y", which='major', direction="in", length=6, right=True)
plt.tick_params(axis="y", which='minor', direction="in", length=3, right=True) 
plt.legend(prop={'size': 6}, title=r"$\rm N_{seg}/N_T$", ncol=2, title_fontproperties={'size': 10}, loc='lower right')
plt.tight_layout()
plt.show()

instr = input("Save figure? (yes/no): ")
if instr == "yes":
    figFolder = input("Folder path to save figure: ")
    plt.savefig(figFolder + "/figure2a_300.png", dpi=300)
    print("Figure saved!")
else:
    print("Figure was not saved!")
    
#%% (optional) Plotting MSD fitting with optimum n points

p = np.polyfit(msdData[3]['dt'][:6], msdData[3]['msd'][:6], 1)

plt.figure(figsize=(2.95, 2.95))
plt.scatter(msdData[3]['dt'][:6], msdData[3]['msd'][:6], color='black', label='n = 6')
plt.plot(msdData[3]['dt'][:6], p[0]*msdData[3]['dt'][:6]+p[1], '--', color='red', zorder=0, label='linear fit')
plt.xlabel(r"$\rm \tau $ (s)", fontsize=14)
plt.ylabel(r"MSD ($\rm \mu m^2 /s$)", fontsize=14)
plt.xticks(fontsize=12) 
plt.yticks(fontsize=12)
plt.locator_params(axis='both', nbins=7)
plt.tick_params(axis="x", which='major', direction="in", length=6, top=True)
plt.tick_params(axis="x", which='minor', direction="in", length=3, top=True)
plt.tick_params(axis="y", which='major', direction="in", length=6, right=True)
plt.tick_params(axis="y", which='minor', direction="in", length=3, right=True)
plt.legend(prop={'size': 8}, loc='lower right') 
plt.tight_layout()
plt.show()

instr = input("Save figure? (yes/no): ")
if instr == "yes":
    figFolder = input("Folder path to save figure: ")
    plt.savefig(figFolder + "/figure2b_300.png", dpi=300)
    print("Figure saved!")
else:
    print("Figure was not saved!")

#%% (optional) Plotting PDF of average slope for optimum points

plotSlopeList = [slopesData["slope3"]["D*3"][5], slopesData["slope3"]["D*4"][optPoint], slopesData["slope3"]["D*5"][optPoint]]

plotSlopeDF = pd.DataFrame({"3": plotSlopeList})
plotSlopeDF.plot.kde(figsize=(2.95, 2.95), legend=False)
plt.xlabel(r"$\rm D^*\ ({\mu m^2/s})$", fontsize=14)
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
    plt.savefig(figFolder + "/figure2c_300.png", dpi=300)
    print("Figure saved!")
else:
    print("Figure was not saved!")

#%% Plotting distributions of D as a function of Nseg

optDCoeffs = []

for i in Nseg:
    if len(slopesData["slope{0}".format(i)]) > optPoint:
        optSlope = slopesData["slope{0}".format(i)].iloc[optPoint, 1:(i + 1)]
        optDCoeffs.append(optSlope.values/4)

fig, ax = plt.subplots()

if len(optDCoeffs[0]) == 1:
    plt.axvline(x = optDCoeffs[0], color = 'gray', linestyle = '--', label='1')
    for entry in optDCoeffs[1:]:
        optDCoeffsDF = pd.DataFrame({'{0}'.format(len(entry)): entry})
        # optDCoeffs[entry].plot.hist(density=True, ax=ax)
        optDCoeffsDF.plot.kde(ax=ax, figsize=(2.95, 2.95))
else:
    for entry in optDCoeffs:
        optDCoeffsDF = pd.DataFrame({'{0}'.format(len(entry)): entry})
        # optDCoeffs[entry].plot.hist(density=True, ax=ax)
        optDCoeffsDF.plot.kde(ax=ax, figsize=(2.95, 2.95))

plt.xlabel(r"D ($\rm {\mu m^2/s})$", fontsize=14)
plt.ylabel("PDF", fontsize=14)
plt.xticks(fontsize=12) 
plt.yticks(fontsize=12)
plt.tick_params(axis="x", direction="in", length=6, top=True)
plt.tick_params(axis="y", direction="in", length=6, right=True)
plt.ylim(0)
plt.legend(prop={'size': 8}, title=r"$\rm N_{seg}/N_T$", title_fontproperties={'size': 10, 'weight':'bold'}, loc='upper right')
plt.tight_layout()
plt.show()

instr = input("Save figure? (yes/no): ")
if instr == "yes":
    figFolder = input("Folder path to save figure: ")
    plt.savefig(figFolder + "/figure3a_300.png", dpi=300)
    print("Figure saved!")
else:
    print("Figure was not saved!")
    
#%% Plotting D_mean and s_D as a function of Nseg

xDataList = []
yDataList = []
yErrorList = []
sigmaList = []
yErrorMaxList = []
yErrorMinList = []
sigmaMaxList = []
sigmaMinList = []


for i in Nseg:
    xData = len(data) // i
    xDataList.append(xData)
    
for dat in slopesData.keys():
    if len(slopesData[dat]) > optPoint:
        optDMean = slopesData[dat]["D*_mean"][optPoint] / 4
        optDStd = slopesData[dat]["D*_std"][optPoint] / 4
        yDataList.append(optDMean)
        yErrorList.append(optDStd)
    
for x in np.array(range(len(yDataList))):
    sigma = yDataList[x] * np.sqrt(2*(optPoint + 1) / (3*(xDataList[x] - (optPoint + 1))))
    sigmaList.append(sigma)
    sigmaMax = yDataList[x] + sigma
    sigmaMin = yDataList[x] - sigma
    yErrorMax = yDataList[x] + yErrorList[x]
    yErrorMin = yDataList[x] - yErrorList[x]
    sigmaMaxList.append(sigmaMax)
    sigmaMinList.append(sigmaMin)
    yErrorMaxList.append(yErrorMax)
    yErrorMinList.append(yErrorMin)
    
plt.figure(figsize=(2.95, 2.95))
if len(optDCoeffs[0]) == 1:
    plt.axhline(y = optDCoeffs[0], color = 'gray', linestyle = '--', label = r'$\rm {D_{T}}$')
    plt.scatter(xDataList[1:len(yDataList)], yDataList[1:], marker ='s', s=9, color = 'black', zorder=3)
    plt.errorbar(xDataList[1:len(yDataList)], yDataList[1:], yerr=yErrorList[1:], ls = 'none', ecolor='black', capsize=3, label=r'$\rm s_{D}$')
    plt.errorbar(xDataList[1:len(yDataList)], yDataList[1:], yerr=sigmaList[1:], ls = 'none', ecolor='red', capsize=3, label=r'$\rm \sigma_Q$')
    # plt.plot(xDataList[1:], yErrorMaxList[1:], linestyle='-.', color='black')
    # plt.plot(xDataList[1:], yErrorMinList[1:], linestyle='-.', color='black')
    # plt.plot(xDataList[1:], sigmaMaxList[1:], linestyle=':', color='red')
    # plt.plot(xDataList[1:], sigmaMinList[1:], linestyle=':', color='red')
    # plt.plot(1000, 0.16550430080922215, color='green', marker='*')
else:
    plt.scatter(xDataList, yDataList, marker ='s', color = 'black', zorder=3)
    plt.errorbar(xDataList, yDataList, yerr=yErrorList, ls = 'none', ecolor='black', capsize=3)
    plt.errorbar(xDataList, yDataList, yerr=sigmaList, ls = 'none', ecolor='red', capsize=3)

plt.xlabel(r"$\rm N_{seg}$", fontsize=14)
plt.ylabel(r"$\rm \overline{D} \ (\mu m^2/s)$", fontsize=14)
plt.xticks(np.arange(0, 600, step=100), fontsize=12) 
plt.yticks(fontsize=12)
plt.locator_params(axis='y', nbins=5)
plt.tick_params(axis="x", direction="in", length=6, top=True)
plt.tick_params(axis="y", direction="in", length=6, right=True)
plt.legend(prop={'size': 8}, loc='upper right')
plt.tight_layout()
plt.show()

instr = input("Save figure? (yes/no): ")
if instr == "yes":
    figFolder = input("Folder path to save figure: ")
    plt.savefig(figFolder + "/figure3b_300.png", dpi=300)
    print("Figure saved!")
else:
    print("Figure was not saved!")

print("D_T = ", optDCoeffs[0])