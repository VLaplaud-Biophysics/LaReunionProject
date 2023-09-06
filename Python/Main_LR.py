# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 10:12:59 2023

@author: laplaud

Main codes for La Reunion data analysis
"""

# Imports 

import numpy as np
import pandas as pd

import os

import VallapFunc_LR as vf


from scipy.interpolate import LinearNDInterpolator as Linterp, NearestNDInterpolator as Ninterp


from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'inline')



def MonthAvg(dfs,AllVars,monthlist):
    
    dffinal = pd.DataFrame(columns=['STATION','ALTITUDE','LATITUDE','LONGITUDE','M_List'])

    for df,Vars in zip(dfs,AllVars):

        # Month selection
        dfmonth = df.loc[vf.ismember(monthlist,df['MONTH'].values)].copy(deep=True)

        StationList = np.unique(df['STATION'])

        for Var in Vars:
            for S in StationList:

                dfmonth.loc[(dfmonth['STATION'] == S)&(dfmonth['MONTH'] == monthlist[0]),Var+'_AVG'] = np.nanmean(dfmonth.loc[df['STATION'] == S][Var])
                dfmonth.loc[(dfmonth['STATION'] == S)&(dfmonth['MONTH'] == monthlist[0]),'M_List'] = str(monthlist)
        dffinal = dfmonth.loc[dfmonth['MONTH'] == monthlist[0],
                                                  ['STATION','ALTITUDE','LATITUDE','LONGITUDE','M_List']+
                                                  [ V + '_AVG' for V in Vars]].merge(dffinal,how='outer',
                                                                                     on=['STATION','ALTITUDE','LATITUDE','LONGITUDE','M_List'])

    return dffinal



def WeatherInterp(DF_plants,DF_weather,columns,**kwargs):
    
    interpmethod = 'linear'
    
    for key, value in kwargs.items(): 
        if key == 'interpmethod':
            interpmethod = value 
    
    if interpmethod == 'linear':
        Interp = Linterp
    elif interpmethod == 'nearest':
        Interp = Ninterp
    else:
        print('Wrong method : ' + interpmethod + '. Option are ''nearest'' or ''linear''. Using ''nearest''.')
        Interp = Ninterp

    xref = DF_weather['Latitude (km)'].values
    yref = DF_weather['Longitude (km)'].values
    zref = DF_weather['Altitude (km)'].values
        
    xq = DF_plants['Latitude (km)'].values
    yq = DF_plants['Longitude (km)'].values
    zq = DF_plants['Altitude (km)'].values
        
    for column in columns:
        
        values = DF_weather[column].values

        interp_weather = Interp(list(zip(xref,yref,zref)),values)  

        DF_plants[column] = interp_weather(xq,yq,zq)
    
    return DF_plants 