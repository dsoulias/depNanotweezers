def mySG(arr, window_size, order, deriv=0, delta=1.0, axis=-1, mode='interp', cval=0.0):
    
    """
    arr:
    1 or 2-dimensional array. If arr is 1D, it contains the data to be smoothed, and the x-axis is considered to 
    be equally spaced with spaces of delta
    
    window_size:
    window_size has to be an odd number, so that half_window = (window_size-1)/2 is an integer. 
    For all data points other than the ones at the ends a polynomial is now fitted to half_window points to the left, 
    the point itself, and half_window points to the right. The value of this polynomial at the point is now taken as 
    the smoothed data point.
    
    order: 
    Order of the polynomial to be used to fit the data.
    
    deriv:
    Return the smoothed deriv'th derivative. The default is to return the smoothed version of the data.
    
    delta:
    Spacing if an equally-spaced x-axis is used (see arr).
    
    axis:
    If arr is 2D, axis defines which one is the x-axis. If axis=-1, then the x-axis is equally spaced with delta.
    All others are smoothed.  
    
    mode:
    The mode parameter defines how the data at the beginning and the end will be handled. 
    It is possible to use different modes at the leading and the tail end of the data. Use the notation:
        leadingMode:tailMode e.g. mode="fit:interp"
    Allowed modes: 
    - mirror: The data array is extended at the beginning and the end with half_window points as follows:
        p[half_window], ... p[2], p[1], p[0], p[1], p[2], ... p[half_window]
      and then all real data points can be fitted as normal
    - pointmirror: The data array is extended at the beginning and the end with half_window points as follows:
        -p[half_window], ... -p[2], -p[1], p[0], p[1], p[2], ... p[half_window]
      and then all real data points can be fitted as normal.
      Note this is similar to mirror other than that the vlaues are taken negative
    - constant: The data array is extended at the beginning and the end with half_window constant value points. 
      The value is cval.
        c, c, c, p[0], p[1], p[2], ...
    - nearest: Like constant, but instead of a cval the value of the first and last point, respectively, is used
        p[0], p[0], p[0], p[0], p[1], p[2], ..., p[n], p[n], p[n], p[n]
    - wrap: the last half_window data points are added to the beginning and the first half_window data points to 
      the end (the data are now circular)
        p[n-2], p[n-1], p[n], p[0], p[1], p[2], ..., p[n-2], p[n-1], p[n], p[0], p[1], p[2]
    - skip: the data array is not extended and the first and last half_window data points are not fitted but 
      are left as are
    - interp: the data array is not extended, and the first and last half_window data points are not fitted but 
      interpolated using the last full fit at the beginning and the end, respectively
    - fit: the data array is not extended, and the first and last half_window data points are fitted with reduced 
      number of points. For example, the first point is evaluated by fitting the the first half_window+1 data points,
      and then using the value of the poynom at position 0. 
      cval can be set to give more weight to points towards the beginning and the end, respectively (as there are now 
      points missing compared to the normal fits). If cval is positive, each point from the relevant position towards
      the beginning (or the end, repectivley) will get a weight of cval compared to 1 for all other half_window data 
      points. If cval is negative, the weight will increase with reducing number of points: 
        w = (-cval)*(half_window+1)/(no of points left to fit)
      The idea is that the combined weight of all points for the standard fitting is (2*half_window+1). Now, e.g. at 
      the beginning, if only the relevant point itself plus one more are left, these two get the all weight that 
      would normally be allocated to the left half_window+1 points. This almost forces the last point to be identical
      to its original. 
      
      An object of the same type as arr is returned, with all data arrays other than x-axis (defined by the axis
      parameter) smoothed.
    
    """

    import math
    import numpy as np

    
    try:
        window_size = np.abs(int(window_size))
        order = np.abs(int(order))
    except (ValueError, msg):
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2 
    
    
    if axis >= 0:
        t = arr[axis]
        yArr = np.copy(arr)
        yArr = np.delete(yArr, axis, 0)  # remove the x-axis from the yArr construct
        
    else:  # axis=-1 so use equidistant x-axis
        t = np.linspace(0, y.size*delta, y.size)
        if (arr.ndim) > 1:
            yArr = np.copy(arr)
        else:
            yArr = np.array([arr])

    modes = mode.split(":")
    if modes.__len__() == 1:
        bMode = modes[0]
        eMode = modes[0]
    elif modes.__len__() == 2:
        bMode = modes[0]
        eMode = modes[1]
    else:
        raise ValueError("mode is not one of the allowed types")

      
    # loop over all data arrays and smooth
    yOutArr = np.array([])

    for y in yArr:

        # work out how to treat the extreme values at the beginning
        if bMode=='mirror': # Repeats the values at the edges in reverse order. The value closest to the edge is not included.
            firstvals = y[1:half_window+1][::-1]
            firstTvals = (t[::-1][1:half_window+1]+t[0]-t[-1])[::-1]
        elif bMode=='pointmirror': # like mirror, but point-mirror at the first and last point
            firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
            firstTvals = (t[::-1][1:half_window+1]+t[0]-t[-1])[::-1]
        elif bMode=='constant': # The extension contains the value given by the cval argument.
            firstvals = np.full(half_window, cval)
            deltaT = (t[2*half_window+1]-t[0])/(2*half_window+1)
            firstTvals = np.linspace(-half_window*deltaT, -deltaT, half_window)
        elif bMode=='nearest': # The extension contains the nearest input value.
            firstvals = np.full(half_window, y[0])
            deltaT = (t[2*half_window+1]-t[0])/(2*half_window+1)
            firstTvals = np.linspace(-half_window*deltaT, -deltaT, half_window)
        elif bMode=='wrap': # The extension contains the values from the other end of the array
            firstvals = y[::-1][0:half_window][::-1]
            firstTvals = (t[::-1][1:half_window+1]+t[0]-t[-1])[::-1]
        elif bMode=='skip': # do not smooth the first and last half_window points
            firstvals = []
            firstTvals = []
        elif bMode=='interp':
            firstvals = []
            firstTvals = []
        elif bMode=='fit':
            firstvals = []
            firstTvals = []
        else:
            raise ValueError("mode for leading-end is not one of the allowed types")
        
        # work out how to treat the extreme values at the end
        if eMode=='mirror': # Repeats the values at the edges in reverse order. The value closest to the edge is not included.
            lastvals = y[::-1][1:half_window+1]
            lastTvals = t[1:half_window+1]-t[0]+t[-1]
        elif eMode=='pointmirror': # like mirror, but point-mirror at the first and last point
            lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
            lastTvals = t[1:half_window+1]-t[0]+t[-1]
        elif eMode=='constant': # The extension contains the value given by the cval argument.
            lastvals = np.full(half_window, cval)
            deltaT = (t[::-1][0]-t[::-1][2*half_window+1])/(2*half_window+1)
            lastTvals = np.linspace(t[::-1][0]+deltaT, t[::-1][0]+half_window*deltaT, half_window)
        elif eMode=='nearest': # The extension contains the nearest input value.
            lastvals = np.full(half_window, y[-1])
            deltaT = (t[::-1][0]-t[::-1][2*half_window+1])/(2*half_window+1)
            lastTvals = np.linspace(t[::-1][0]+deltaT, t[::-1][0]+half_window*deltaT, half_window)
        elif eMode=='wrap': # The extension contains the values from the other end of the array
            lastvals = y[0:half_window]
            lastTvals = t[1:half_window+1]-t[0]+t[-1]
        elif eMode=='skip': # do not smooth the first and last half_window points
            lastvals = []
            lastTvals = []
        elif eMode=='interp':
            lastvals = []
            lastTvals = []
        elif eMode=='fit':
            lastvals = []
            lastTvals = []
        else:
            raise ValueError("mode for tail-end is not one of the allowed types")

        y = np.concatenate((firstvals, y, lastvals))
        t = np.concatenate((firstTvals, t, lastTvals))

        yOut = np.array([])
        for i in range(half_window, y.size - half_window): # loop over all points to fit now
            tForFit = t[i-half_window:i+half_window+1]
            yForFit = y[i-half_window:i+half_window+1]
            param = np.polyfit(tForFit, yForFit, order, rcond=None, full=False, w=None, cov=False)
            yFit = 0
            for e in range(0, order+1-deriv):
                yFit = yFit + math.factorial(e+deriv)/math.factorial(e)*param[::-1][e+deriv]*t[i]**e
 
            yOut = np.append(yOut, yFit)
    
            if i==-1:
                #print(param)
                #print(tForFit, yForFit, np.polyval(param, tForFit))
                plt.plot(np.linspace(tForFit[0], tForFit[-1], 50), np.polyval(param, np.linspace(tForFit[0], tForFit[-1], 50)))

        
        if bMode=='skip': # do not smooth the first half_window points
            firstvals = y[0:half_window]
            yOut = np.concatenate((firstvals, yOut))
            
        if eMode=='skip': # do not smooth the last half_window points
            lastvals = y[::-1][0:half_window][::-1]
            yOut = np.concatenate((yOut, lastvals))
            

        if bMode=='interp': # use the last full fit to smooth the first and last half_window points
            tForFit = t[0:2*half_window+1]
            yForFit = y[0:2*half_window+1]
            param = np.polyfit(tForFit, yForFit, order, rcond=None, full=False, w=None, cov=False)
            # show fit at the start
            #plt.plot(np.linspace(tForFit[0], tForFit[-1], 50), np.polyval(param, np.linspace(tForFit[0], tForFit[-1], 50)),"g:")
            firstvals = np.array([])
            for i in range(0,half_window):
                yFit = 0
                for e in range(0, order+1-deriv):
                    #yFit = yFit + param[::-1][e]*t[i]**e
                    yFit = yFit + math.factorial(e+deriv)/math.factorial(e)*param[::-1][e+deriv]*t[i]**e

                firstvals = np.append(firstvals, yFit)
            yOut = np.concatenate((firstvals, yOut))


        if eMode=='interp': # use the last full fit to smooth the first and last half_window points
            tForFit = t[-(2*half_window+1):]
            yForFit = y[-(2*half_window+1):]
            param = np.polyfit(tForFit, yForFit, order, rcond=None, full=False, w=None, cov=False)
            #print(param)
            #print(tForFit, yForFit, np.polyval(param, tForFit))
            # show fit at the end
            #plt.plot(np.linspace(tForFit[0], tForFit[-1], 50), np.polyval(param, np.linspace(tForFit[0], tForFit[-1], 50)), "g:")
            lastvals = np.array([])
            for i in range(half_window,0,-1):
                yFit = 0
                for e in range(0, order+1-deriv):
                    yFit = yFit + math.factorial(e+deriv)/math.factorial(e)*param[::-1][e+deriv]*t[-i]**e

                lastvals = np.append(lastvals, yFit)
            yOut = np.concatenate((yOut, lastvals))
      

        if bMode=='fit': # smooth the first half_window points fitting with reduced number of points
            weightEnhancement = cval
            firstvals = np.array([])
            for i in range(0,half_window):
                tForFit = t[0:half_window+1+i]
                yForFit = y[0:half_window+1+i]
                if weightEnhancement<0:
                    weight = -(half_window+1)*weightEnhancement/(i+1)
                elif weightEnhancement>0:
                    weight = weightEnhancement
                else:
                    weight = 1

                wForFit = np.concatenate( (np.full((i+1), weight), np.full((half_window), 1)) )
                param = np.polyfit(tForFit, yForFit, order, rcond=None, full=False, w=wForFit, cov=False)
                # show fits at the start
                #plt.plot(np.linspace(tForFit[0], tForFit[-1], 50), np.polyval(param, np.linspace(tForFit[0], tForFit[-1], 50)), ":")

                yFit = 0
                for e in range(0, order+1-deriv):
                    yFit = yFit + math.factorial(e+deriv)/math.factorial(e)*param[::-1][e+deriv]*tForFit[i]**e
                firstvals = np.append(firstvals, yFit)
            yOut = np.concatenate((firstvals, yOut))

        
        if eMode=='fit': # smooth the first half_window points fitting with reduced number of points
            weightEnhancement = cval
            lastvals = np.array([])
            for i in range(half_window,0,-1):
                tForFit = t[-(half_window+i):]
                yForFit = y[-(half_window+i):]
                if weightEnhancement<0:
                    weight = -(half_window+1)*weightEnhancement/(i)
                elif weightEnhancement>0:
                    weight = weightEnhancement
                else:
                    weight = 1
                    
                wForFit = np.concatenate( (np.full((half_window), 1), np.full((i), weight)) )
                #print(weightEnhancement, tForFit, wForFit)
                param = np.polyfit(tForFit, yForFit, order, rcond=None, full=False, w=wForFit, cov=False)
                # show fits at the end
                #plt.plot(np.linspace(tForFit[0], tForFit[-1], 50), np.polyval(param, np.linspace(tForFit[0], tForFit[-1], 50)), ":")
            
                yFit = 0
                for e in range(0, order+1-deriv):
                    yFit = yFit + math.factorial(e+deriv)/math.factorial(e)*param[::-1][e+deriv]*t[-i]**e
                lastvals = np.append(lastvals, yFit)
            yOut = np.concatenate((yOut, lastvals))

            
        if yOutArr.size == 0:  # contains no elements yet
            yOutArr = np.array(yOut)
        else:
            yOutArr = np.append([yOutArr], [yOut], axis=0)

            

    if axis >= 0:
        if yOutArr.ndim > 1:
            return np.insert(yOutArr, axis, t, axis=0)
        else:
            if axis == 0:
                return np.array([t, yOutArr])
            else:
                return np.array([yOutArr, t])
    else:
        return yOutArr