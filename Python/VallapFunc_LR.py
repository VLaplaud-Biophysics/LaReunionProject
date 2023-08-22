# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 10:12:59 2023

@author: laplaud

Support functions for La Reunion data analysis
"""

import numpy as np
import pandas as pd
import mpmath as mpm

import matplotlib.pyplot as plt
import seaborn as sns

from scipy.spatial.distance import directed_hausdorff 
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d

import cv2 as cv

#  1. Compute normal to a vector

# Normal vector is normalized, and by default rotated counter clockwise from given vector
# the input vector is from (x1,y1) to (x2,y2)


def getNormal(x1,y1,x2,y2, **kwargs):
    
    rotation = 'CCW'
    
    for key, value in kwargs.items(): 
        if key == 'rotation':
            rotation = value
        else:            
            print('Unknown key : ' + key+ '. Kwarg ignored.')
            
    dx = x2 - x1
    dy = y2 - y1
    
    Norm = np.sqrt(np.square(dx) + np.square(dy))
    
    dxN = np.divide(dx,Norm)
    dyN = np.divide(dy,Norm)
    
    if rotation == 'CCW':
        
        x = -dyN
        y = dxN
        
    elif rotation == 'CW':
        
        x = dyN
        y = -dxN

    else:
        print('Wrong rotation parameter !! Should be ''CW'' or ''CCW''. Default is ''CCW'' with no input.')
            
    return(x,y) 
    
#  2. Plotting boxplots with data points on top 

# A function combining boxplot with seaborn's swarmplot for a better display of data

def boxswarmplot(Title,Ylabel,Data,facecolors,Labels,**kwargs):

    FS = (5,3)    

    for key, value in kwargs.items(): 
        if key == 'figsize':
            FS = value

    fig,ax = plt.subplots(dpi = 250,facecolor='white',figsize = FS)
    fig.suptitle(Title)
 
    
    cap= [None]*len(Data)
    med= [None]*len(Data)
    
    grouping = []
    
    for dat,col,lab,i in zip(Data,facecolors,Labels,range(len(Data))):
    
        # plots properties
        plotprops = {'color':'black'}
        boxprops = {'color':'black','facecolor':col}
        
        lab = lab + '\nn = ' + str(len(dat))
        
        Labels[i] = lab

        bp = ax.boxplot(dat, positions = [i], labels = [lab],patch_artist =True, boxprops=boxprops, capprops =plotprops,
                    showfliers=False,whiskerprops=plotprops,medianprops =plotprops)
        
        grouping = np.append(grouping,np.ones(len(dat))*i)
    
        cap[i] = bp['caps'][1].get_ydata(orig=True)[0]
        med[i] = bp['medians'][0].get_ydata(orig=True)[0]
    
    sns.swarmplot(x=grouping,y=pd.concat(Data),color = 'gray', size=2, ax = ax)
    
    ax.set_ylabel(Ylabel)
    
    ax.set_xticklabels(Labels)
    
    return(fig,ax,cap,med)
    

# 3. Coordinate conversion from cartesian to circular (in deg) an vice versa

def ToCirc(X,Y, **kwargs):
    
    Angle = 'rad'
    
    for key, value in kwargs.items(): 
        if key == 'angle':
            Angle = value
        else:            
            print('Unknown key : ' + key+ '. Kwarg ignored.')
    
    
    if Angle == 'deg':
        Alpha = np.rad2deg(np.arctan2(Y,X))
    elif Angle == 'rad':
        Alpha = np.arctan2(Y,X)
    else:
        print('Wrong angle unit : ' + Angle + '. Default to radians.')         
        Alpha = np.arctan2(Y,X)
        
    Radius = np.sqrt(np.square(X)+np.square(Y))
    
    return(Alpha,Radius)



def ToCart(Alpha,Radius, **kwargs):
    
    Angle = 'rad'
    
    for key, value in kwargs.items(): 
        if key == 'angle':
            Angle = value
        else:            
            print('Unknown key : ' + key + '. Kwarg ignored.')
    
    if Angle == 'deg':
        Alpharad = np.deg2rad(Alpha)
    elif Angle == 'rad':
        Alpharad = Alpha
    else:
        print('Wrong angle unit : ' + Angle + '. Default to radians.') 
        Alpharad = Alpha
    
    
    X = Radius*np.cos(Alpharad)
    Y = Radius*np.sin(Alpharad)
    
    return(X,Y)

# 4.1 Euclidian distance between two arrays of points in carthesian coordinates 
def dist(x1,y1,x2,y2):
    
    d = np.sqrt(np.square(x1-x2)+np.square(y1-y2))
    
    return(d)

# 5. simple ismember function, checks if A is within B
def ismember(A, B):
    bool_list = list(map(bool,[ np.sum(b == A) for b in B ]))
    return bool_list

# 6. R2 computation for a fit
def computeR2(Ydata,Yfit):
    # Ydata are the fitted data, Yfit the comuted value from the fit
    
    SumResidues = np.sum(np.square(np.subtract(Ydata,Yfit)))
    TotalVariance = np.sum(np.square(np.subtract(Ydata,np.mean(Ydata))))
    
    R2 = 1 - SumResidues/TotalVariance
    
    return R2
  
# 7. Creation of a list to define a mosaic subplot figure

def mosaicList(n):
    alphabet = 'abcdefghijklmn'
    list1 = [*alphabet[0]*n]
    list2 = [*alphabet[1:n+1]]
    mosaic = [list1[:]]
    for i in range(2):
        mosaic.append(list1)
    mosaic.append(list2)
    return(mosaic,list2) 

# 8. Function generating a summary of data and their variability for a specific dataframe column

def dataSummary(GDs,Ns,labels,Mult,col,name,unit):
    DataPooled = np.empty(0)
    nPooled = np.sum(Ns)
    
    print(name + ' : ')
    
    for GD,n,lab in zip(GDs,Ns,labels):
        DataMedian =  np.round(GD.loc[GD['Img'] == 0,col].median()*Mult*100)/100
        Var =  np.round(   np.mean(np.abs(GD.loc[GD['Img'] == 0,col]*Mult-DataMedian))*10000/DataMedian)/100
        DataPooled = np.append(DataPooled,GD.loc[GD['Img'] == 0,col].to_numpy())
        print(lab + ' -> ' + str(DataMedian) + ' ' + unit + ' ' + u"\u00B1" + ' ' + str(Var) + ' %' + ' % (n = ' + str(n) + ')')
    
    PooledMedian = np.round(np.median(DataPooled)*Mult*100)/100
    PooledVar = np.round(np.mean(np.abs(DataPooled*Mult - PooledMedian))*10000/PooledMedian)/100
    print('Pooled -> ' + str(PooledMedian) + ' ' + unit + ' ' + u"\u00B1" + ' ' + str(PooledVar) + ' %' + ' % (n = ' + str(nPooled) + ')' )
    