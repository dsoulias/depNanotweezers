import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
import os
import sys
sys.path.append("F:/PhD/Scripts/DEP/")
from mySG import mySG

np.set_printoptions(precision=2)  # For compact display
matplotlib.rcParams['pdf.fonttype'] = 42 # Nature's guidance for python graphs
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'

class ParticleTrace:
    
    def __init__(self, t, x, y, param):
        self.x = x - x[-1]
        self.y = y - y[-1]
        self.t = t - t[0]
        self.r = np.sqrt(self.x*self.x + self.y*self.y)
        self.rSmooth = self.r.copy()
        self.rDerSmooth = np.zeros(len(self.r))
        self.rDerSmoothSmooth = np.zeros(len(self.r))
        self.DEPForcerSmooth = np.zeros(len(self.r))

    def smoothData(self, window_size, order, mode='interp', cval=0.0, iStart=0, iStop=0):
        if iStop == 0:
            iStop = len(self.r)
        self.rtSmooth = mySG(np.array([self.t[iStart:iStop], self.r[iStart:iStop]]), window_size, order, deriv=0, delta=1.0, axis=0, mode=mode, cval=0.0)
        self.rSmooth[iStart:iStop] = self.rtSmooth[1]

    def smoothDerivative(self, window_size, order, mode='interp', cval=0.0, iStart=0, iStop=0):
        if iStop == 0:
            iStop = len(self.rSmooth)
        self.rtDerSmooth = mySG(np.array([self.t[iStart:iStop], self.rSmooth[iStart:iStop]]), window_size, order, deriv=1, delta=1.0, axis=0, mode=mode, cval=0.0)
        self.rDerSmooth[iStart:iStop] = self.rtDerSmooth[1]

    def smoothSmoothDerivative(self, window_size, order, mode='interp', cval=0.0, iStart=0, iStop=0):
        if iStop == 0:
            iStop = len(self.rSmooth)
        self.rtDerSmoothSmooth = mySG(np.array([self.t[iStart:iStop], self.rDerSmooth[iStart:iStop]]), window_size, order, deriv=0, delta=1.0, axis=0, mode=mode, cval=0.0)
        self.rDerSmoothSmooth[iStart:iStop] = self.rtDerSmoothSmooth[1]
    
    def calcDEPForce(self):
        self.DEPForce = -6*np.pi*param['particleRadius']*param['viscosity']*(param['xyConversion']/param['tConversion'])*self.rDerSmooth
        
    def smoothDEPForce(self):
        self.DEPForceSmooth = -6*np.pi*param['particleRadius']*param['viscosity']*(param['xyConversion']/param['tConversion'])*self.rDerSmoothSmooth

#%% Loading full trajectory

mainFolder = input("Main folder path: ") # input the address of the folder with the "...spots.csv"

for file in os.listdir(mainFolder):
    if file.endswith(".csv"):
        filename = os.path.join(mainFolder, file)
        print(filename)
        x, y, t, f = np.loadtxt(filename, delimiter=',', unpack=True, dtype='float', skiprows=1, usecols=(4,5,7,8)) # load data
        r = np.sqrt(x**2 + y**2)
    
#%% Plotting full trajectory with trapping region

indexStart = int(input("Frame DEP trapping starts: ")) - 1
indexStop = int(input("Frame DEP trapping stops: ")) - 1

plt.figure()
plt.plot(f, r)
plt.axvline(x = f[indexStart], color = 'red', linestyle = '--', label = r'$t_i$')
plt.axvline(x = f[indexStop], color = 'red', linestyle = '--', label = r'$t_f$')
plt.xlabel("Number of frames", fontsize=11, labelpad=16)
plt.ylabel("Displacement, r (\u03BCm)", fontsize=11, labelpad=16)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.tick_params(axis="x", which='major', direction="in", length=6, top=True)
plt.tick_params(axis="x", which='minor', direction="in", length=6, top=True)
plt.tick_params(axis="y", which='major', direction="in", length=6, right=True)
plt.tick_params(axis="y", which='minor', direction="in", length=6, right=True) 
plt.tight_layout()
plt.show()

#%% Exporting trajectory for trapping frames range

spotsDEP = {'x': x[indexStart:indexStop],
            'y': y[indexStart:indexStop],
            't': t[indexStart:indexStop]
               }

spotsDEPdf = pd.DataFrame(spotsDEP)

spotsDEPdf.to_csv(mainFolder + r'_spotsDEP.csv')

#%% Data smoothing for distance from tip

myFolder = input("Folder path: ") # input the address of the folder of interest

param = {'viscosity': (1.53*0.88e-3), 'particleRadius': 1e-6, 'xyConversion': 1e-6, 
         'tConversion': 1, 'particleDensity': 1055}

fig0 = plt.figure(0)
fig0.clf()

for file in os.listdir(myFolder):
    if file.endswith(".csv"):
        filename = os.path.join(myFolder, file)
        print(filename)
        x, y, t = np.loadtxt(filename, delimiter=',', unpack=True, dtype='float', skiprows=1, usecols=(1,2,3)) # load data from splitted trajectory
        # x, y, t = np.loadtxt(filename, delimiter=',', unpack=True, dtype='float', skiprows=1, usecols=(4,5,7)) # load data from full trajectory
        
        trace = ParticleTrace(t, x, y, param) # call class
        trace.smoothData(7, 2, mode='interp') # smooth data
        plt.figure(0)
        ax0 = fig0.subplots()
        p0 = ax0.plot(trace.t, trace.r, '.', color='black') # plot raw data
        p0, = ax0.plot(trace.t, trace.rSmooth, '-', color='green') # plot smoothed data
        plt.subplots_adjust(bottom=0.25, right=0.84) # create space in figure window for subplots
        plt.xlabel("Time (s)")
        plt.ylabel("Displacement (\u03BCm)")
        
        # Designing the slider & radio buttons

        ax0_window = plt.axes([0.25, 0.1, 0.6, 0.03]) # x-position, y-position, width & height
        ax0_order = plt.axes([0.25, 0.05, 0.3, 0.03]) # x-position, y-position, width & height
        ax0_mode = plt.axes([0.85, 0.4, 0.18, 0.35]) # x-position, y-position, width & height

        # Slider & radio buttons properties

        if (len(trace.r) % 2) == 0:
            windowMax = len(trace.r) - 1
        else:
            windowMax = len(trace.r)
        window_0 = Slider(ax0_window, 'Window size', valmin=7, valmax=windowMax, valinit=7, valstep=2)
        order_0 = Slider(ax0_order, 'Polynomial Order', valmin=2, valmax=4, valinit=2, valstep=1)
        mode_0 = RadioButtons(ax0_mode, ('interp', 'fit', 'wrap', 'skip', 'nearest', 'constant', 'mirror'), activecolor='red')

        # Updating the plot

        def update1(val):
            trace.smoothData(window_0.val, order_0.val, mode='interp')
            new_rSmooth = trace.rSmooth
            p0.set_ydata(new_rSmooth)
            fig0.canvas.draw() # redraw the figure

        window_0.on_changed(update1) # calling the function "update" when the window size is changed    
        order_0.on_changed(update1) # calling the function "update" when the polynomial order is changed
        
        def update2(label):
            modes = {'interp': 'interp', 'fit': 'fit', 'wrap': 'wrap', 'skip': 'skip', 
                     'nearest': 'nearest', 'constant': 'constant', 'mirror': 'mirror'}
            trace.smoothData(window_0.val, order_0.val, mode=modes[label])
            new_rSmooth = trace.rSmooth
            p0.set_ydata(new_rSmooth)
            fig0.canvas.draw()
            
        mode_0.on_clicked(update2) # calling the function "update" when the mode is changed        
        plt.show()
        
# Notes:
    # ALWAYS FOLLOW this sequence: select FIRST the window size, SECOND the polynomial order and THIRD the fitting mode
    # REMEMBER to input the desired values and mode in the following cell to save them
   
#%% Plotting distance from tip over time

trace.smoothData(input('Window size = '), input('Polynomial order = '), mode=str(input('Mode: '))) # smooth data

fig1 = plt.figure(1, figsize=(7, 7))
fig1.clf()

plt.plot(trace.t, trace.r, '.', color='black') # plot raw data
plt.plot(trace.t, trace.rSmooth, '-', color='green') # plot smoothed data
plt.xlabel("Time (s)", fontsize=16, labelpad=8)
plt.xticks(fontsize=16)
plt.ylabel("Distance from tip (\u03BCm)", fontsize=16, labelpad=8)
plt.yticks(fontsize=16)
plt.show()

#%% Data smoothing for DEP force

fig2 = plt.figure(2)
fig2.clf()

trace.smoothDerivative(7, 2, mode='interp') # smooth the first derivative of the data
trace.smoothSmoothDerivative(7, 2, mode='interp') # smooth the smoothed data
trace.calcDEPForce()
trace.smoothDEPForce()
ax2 = fig2.subplots()
p2 = ax2.plot(trace.rSmooth, trace.DEPForce,'.', color='black') # plot raw data
p2, = ax2.plot(trace.rSmooth, trace.DEPForceSmooth, '-', color='green') # plot smoothed data
plt.subplots_adjust(bottom=0.25, right=0.84) # create space in figure window for subplots
plt.xlabel("Distance from tip (\u03BCm)")
plt.ylabel("Force (N)")

# Designing the slider & radio buttons

ax2_window = plt.axes([0.25, 0.1, 0.6, 0.03]) # x-position, y-position, width & height
ax2_order = plt.axes([0.25, 0.05, 0.3, 0.03]) # x-position, y-position, width & height
ax2_mode = plt.axes([0.85, 0.4, 0.18, 0.35]) # x-position, y-position, width & height

# Slider & radio buttons properties

if (len(trace.r) % 2) == 0:
    windowMax = len(trace.r) - 1
else:
    windowMax = len(trace.r)
    
window_2 = Slider(ax2_window, 'Window size', valmin=7, valmax=windowMax, valinit=7, valstep=2)
order_2 = Slider(ax2_order, 'Polynomial Order', valmin=2, valmax=4, valinit=2, valstep=1)
mode_2 = RadioButtons(ax2_mode, ('interp', 'fit', 'wrap', 'skip', 'nearest', 'constant', 'mirror'), activecolor='red')

# Updating the plot

def update3(val):
    trace.smoothSmoothDerivative(window_2.val, order_2.val, mode='interp')
    trace.smoothDEPForce()
    new_DEPForce = trace.DEPForceSmooth
    p2.set_ydata(new_DEPForce)
    fig2.canvas.draw() # redraw the figure

window_2.on_changed(update3) # calling the function "update" when the window size is changed    
order_2.on_changed(update3) # calling the function "update" when the polynomial order is changed

def update4(label):
    modes = {'interp': 'interp', 'fit': 'fit', 'wrap': 'wrap', 'skip': 'skip', 
             'nearest': 'nearest', 'constant': 'constant', 'mirror': 'mirror'}
    trace.smoothSmoothDerivative(window_2.val, order_2.val, mode=modes[label])
    trace.smoothDEPForce()
    new_DEPForce = trace.DEPForceSmooth
    p2.set_ydata(new_DEPForce)
    fig2.canvas.draw()
    
mode_2.on_clicked(update4) # calling the function "update" when the mode is changed 
plt.show()

#%% Plotting DEP force over displacement

trace.smoothSmoothDerivative(input('Window size = '), input('Polynomial order = '), mode=str(input('Mode: ')))
trace.calcDEPForce()
trace.smoothDEPForce()

fig3 = plt.figure(3, figsize=(6.68, 5.65))
fig3.clf()

plt.plot(trace.rSmooth, trace.DEPForce*1e15, '.', color='black') # plot raw data
plt.plot(trace.rSmooth, trace.DEPForceSmooth*1e15, '-', color='green') # plot smoothed data
plt.xlabel("Displacement (\u03BCm)", fontsize=16, labelpad=8)
plt.xticks(fontsize=16)
plt.ylabel("Force (fN)", fontsize=16, labelpad=8)
plt.yticks(fontsize=16)
plt.tight_layout()
plt.show()

print("F_DEP_max (fN) =", round(max(trace.DEPForce*1e15), 2))
print("F_DEP_max_smooth (fN) =", round(max(trace.DEPForceSmooth*1e15), 2))

#%% Exporting DEP data

export_data = {'rSmooth(um)':trace.rSmooth,
               'DEPForceSmooth(fN)':trace.DEPForceSmooth*1e15}

df_export_data = pd.DataFrame(export_data)
df_export_data.to_csv(myFolder + r'/_depForce.csv')