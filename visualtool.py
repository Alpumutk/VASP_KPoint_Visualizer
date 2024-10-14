#!/usr/bin/env python
# coding: utf-8

#import libraries

from __future__ import print_function
import pandas as pd
import numpy as np 
from numpy import *
import sys, os
import matplotlib.pyplot as plt
import time
import ipywidgets as widgets
import plotly
import plotly.graph_objects as go
import warnings
#warnings.filterwarnings("ignore", category=DeprecationWarning) #ignore Depreciation Warning


def Transmatrix(transmatrix_file):
    """Reads in transmatrix file and creates dataframe with values.
    
    :param transmatrix_file: Path of transmatrix file
    :return trans_df: Dataframe with transmatrix values
    """
    
    trans_df = pd.read_csv(transmatrix_file,sep='\s+',names = ['kweight','cband','vband','Ec','Ev','px','ipx','py','ipy','pz','ipz','occ'])
    return trans_df


def Initialize(outcarfile, trans_df, Enmin, Enmax, cbandmin, cbandmax,vbandmin,vbandmax, klist):
    """This function needs to be run once to initialize data, which modifies the input transmatrix dataframe to include coordinates and kpoints
    
    :param outcarfile: Path of outcar file
    :param trans_df: Dataframe with transmatrix values
    :param cbandmin: Minimum conduction band value to cut the transmatrix by
    :param cbandmax: Maximum conduction band value to cut the transmatrix by
    :param vbandmin: Minimum valence band value to cut the transmatrix by
    :param vbandmax: Maximum valence band value to cut the transmatrix by
    :param klist: A list of k-points to cut the transmatrix by
    :return trans_df_cut: Modified dataframe with transmatrix values
    :return k_df: Dataframe with positions of k-points 
    :return VBM: Band number of the valence band maximum
    :return kpoint_number: Number of k-points
    :return band_number: Number of bands
    """
    trans_df['Ediff'] = trans_df['Ec']-trans_df['Ev'] #creates a column that contains the energy difference between the conduction and valence bands
    with open(outcarfile) as f: #open OUTCAR file
        lines = list(f) #Read all lines in as a list
        kpoint_number = [int(line.split()[3]) for n,line in enumerate(lines,1) if "NKPTS" in line][0] #Obtain number of k-points
        band_number = [int(line.split()[14]) for n,line in enumerate(lines,1) if "NKPTS" in line][0] #Obtain number of bands
        kpt_ln = [n for n,line in enumerate(lines,1) if "Following reciprocal coordinates:" in line][0] #Obtain the line that preceeds the coordinates of all k-points
        E_header_line_number = next((n for n, line in enumerate(lines, 1) if " band No." in line)) #Obtain the first line that contains "Band No.", will be used to get VBM
  
    
    VBM=[float(np.loadtxt(outcarfile,usecols=2,skiprows=E_header_line_number+a,max_rows=1)) for a in range(band_number)].index(0) #read through all band occupations to find the index of the first unoccupied band, this will equal VBM as band no. 1 = 0th index
    tm_linenumber = VBM*(band_number-VBM) #calculate the number of lines in the transmatrix per kpoint
    k_df=pd.read_csv(outcarfile,skiprows=kpt_ln+1,nrows=kpoint_number,sep='\s+',usecols=[0,1,2],names=["xcoord","ycoord","zcoord"]) #read in the locations of k-points as a dataframe
    k_df = k_df.add(0.5) #shift all coordinates by 0.5 to move Γ-point to (0.5,0.5,0.5)
    kpt_no = [i for i in range(len(k_df))]
    k_df['k index'] = kpt_no

    trans_df['k index'] = np.repeat(np.array([i for i in range(kpoint_number)]) , tm_linenumber) #add a column to trans_df that contains k-point no.
    trans_df['xcoord'], trans_df['ycoord'], trans_df['zcoord']  = np.repeat(k_df['xcoord'].to_numpy(), tm_linenumber), np.repeat(k_df['ycoord'].to_numpy(), tm_linenumber),  np.repeat(k_df['zcoord'].to_numpy(), tm_linenumber) #add columns to trans_df that contain coordinates
    
    
    trans_df_cut =  trans_df[(trans_df['Ediff']<=Enmax) & (trans_df['Ediff']>=Enmin) & (trans_df['cband']<=cbandmax) & (trans_df['cband']>=cbandmin) & (trans_df['vband']<=vbandmax)& (trans_df['vband']>=vbandmin)] #Cuts trans_df according to user input ranges 
    trans_df_cut = trans_df_cut[trans_df_cut['k index'].isin(klist)]
    prob = np.sqrt(trans_df_cut['px']**2+trans_df_cut['ipx']**2)+np.sqrt(trans_df_cut['py']**2+trans_df_cut['ipy']**2)+np.sqrt(trans_df_cut['pz']**2+trans_df_cut['ipz']**2) #Calculate transition probability by summing x,y,z components
    trans_df_cut.insert(len(trans_df_cut.columns),"prob",prob) #append transition probability as a column for trans_df
    trans_df_cut = trans_df_cut.reset_index(drop=True) #resets indices after modification
    
    
    return trans_df_cut, k_df,VBM, kpoint_number,band_number #return outputs


def Dielectric(trans_df_cut, Enmin, Enmax, kpoint_number, Volume = 200.0, sigma = 0.10, GRID = 10000, scissor = 0.0, SOC = 'FALSE'): 
    
    """This function modifies the trans_df according to user inputted range on bands, kpoints and energies. It calculates imaginary dielectric spectrum and energy in this range and in a grid determined by user input values.
    
    :param trans_df_cut: Dataframe with transmatrix values
    :param Enmin: Minimum energy value of the calculated dielectric spectrum
    :param Enmax: Maximum energy value of the calculated dielectric spectrum
    :param kpoint_number: Number of k-points
    :param Volume: Volume in reciprocal space
    :param sigma: Width of smearing
    :param GRID: Number of grid points in calculated dielectric spectrum
    :param scissor: Scissor correction for the bandgap
    :param SOC: Option for enabling spin-orbit coupling
    :return Energy: List with energy values
    :return tot_curve_ave: List with dielectric spectrum
    """
    
    #Define constants
    hr = 27.2116
    bohr = 0.529177249
    pi = 3.14159
    KB = 8.61733E-5

    start_time = time.time() #keeps track of time

    ######## Constants and Coefficients for Dielectric Function Calc ########
    fac = 8.*pi**2*hr**3*bohr**3/Volume
    hbareta = 4.*0.69314718/(sigma**2)
    lfac = 2*fac/pi
    ##########################################################################
    ######## Prep Grids for dielectric function ########
    Energy = linspace(Enmin,Enmax,GRID)
    tot_curve_ave = [0]*len(Energy)
    ###################################################
    ######## Get VBM and CBM from Transmatrix File ########
    VBM_en, CBM_en   = trans_df_cut['Ev'].max(), trans_df_cut['Ec'].min()
    print('Valence band maximum and conduction band minimum from Transmatrix are:   ',VBM_en, CBM_en)
    print('(CBM-VBM)/2 = ',(VBM_en+CBM_en)/2.)
    ########################################################

    #Get matrix elements and dielectric function
    for i in range(len(trans_df_cut)): #loop over modified transmatrix
        PX, PY, PZ = np.sqrt(trans_df_cut['px'][i]**2+trans_df_cut['ipx'][i]**2), np.sqrt(trans_df_cut['py'][i]**2+trans_df_cut['ipy'][i]**2), np.sqrt(trans_df_cut['pz'][i]**2+trans_df_cut['ipz'][i]**2)
        dE = trans_df_cut['Ec'][i]-trans_df_cut['Ev'][i] #calculate energy difference between conduction and valence bands
        ddE = dE+scissor
        MX2, MY2, MZ2 = trans_df_cut['kweight'][i]*((ddE/dE)*(ddE/dE))*(PX/ddE)**2, trans_df_cut['kweight'][i]*((ddE/dE)*(ddE/dE))*(PY/ddE)**2, trans_df_cut['kweight'][i]*((ddE/dE)*(ddE/dE))*(PZ/ddE)**2 
        MA2 = (MX2+MY2+MZ2)/3 
        av = np.imag(MA2*lfac*ddE/(ddE**2-(Energy+1j*sigma)**2))
        tot_curve_ave += av
       


    if SOC=='TRUE': #check for spin-orbit coupling condition
        tot_curve_ave = tot_curve_ave/2.


    ### Write spectrum to file ##########################
    file = open('outputspectrum', "w")
    for index in range(len(Energy)):
        file.write(str(Energy[index])+" "+str(tot_curve_ave[index])+"\n")
    file.close()

    print("--- %s seconds to completion ---" % (time.time() - start_time))
    ###############################################
    
   
    
    return Energy, tot_curve_ave #return outputs

def PlotDielectric(Energy, tot_curve_ave):
    """This function plots the dielectric spectrum vs. energy
    
    :param Energy: List with energy values
    :param tot_curve_ave: List with dielectric spectrum
    """
    
    plt.plot(Energy, tot_curve_ave)
    plt.xlabel('Energy (eV)')
    plt.ylabel('ε2')
    plt.show()
    

def RunTransition(trans_df, ev,deltaE, option, kpt_number, k_df):
    """This function transforms the trans_df such that only transitions within an energy range of (E-deltaE, E+deltaE) are considered. Transition probabilities are calculated for each transition and added to trans_df.
    
    :param trans_df: Dataframe with transmatrix values
    :param ev: Transition energy value to cut the transmatrix by 
    :param deltaE: Determines the range of transition energies(E-deltaE, E+deltaE)  to be considered
    :param option: When set to 'yes', sums all of the transition probabilities for a given k-point
    :param kpt_number: Number of k-points
    :param k_df: Dataframe with positions of k-points 
    :return trans_df_plot: Modified transmatrix dataframe for plotting
    :return max_prob: Maximum oscillator strength for plotting
    """

  
        
  
    max_prob = max(trans_df['prob']) #calculate the maximum oscillator strength in the transmatrix for plotting purposes
    trans_df_cut =  trans_df[(trans_df['Ediff']<=(ev+(deltaE))) & (trans_df['Ediff']>=(ev-(deltaE)))]
    trans_df_cut = trans_df_cut.reset_index(drop=True) #Reset index number 
    columns_to_remove= ['kweight','Ec', 'Ev', 'px','ipx','py','ipy','pz','ipz','occ', 'Ediff']  # List of excess columns to remove for plotting
    trans_df_plot = trans_df_cut.drop(columns=columns_to_remove) #Remove columns that will not be used in plotting
    
        
    if option == 'yes': #check for condition where all transition probabilities at a k-point are summed
        summedproblist = [probsum for i in range(kpt_number) if (probsum := sum(trans_df_plot['prob'][j] for j in range(len(trans_df_plot['k index'])) if i == trans_df_plot['k index'][j])) != 0] #loop over all kpoints and check if they match the k indices in trans_df_plot. Probabilities of repeating k-points are then summed and appended to a list if it does not equal zero.
        max_prob = max(summedproblist) #calculate the maximum oscillator strength in the transmatrix for plotting purposes
        coord = [[trans_df_plot['xcoord'][j], trans_df_plot['ycoord'][j], trans_df_plot['zcoord'][j]] for i in range(kpt_number) for j in range(len(trans_df_plot['k index'])) if i == trans_df_plot['k index'][j]] #loop over all kpoints and check if they match the k indices in trans_df_plot. Append x,y,z coordinates of those that match to a list.
        coord_unique = np.array([coord[i] for i in range(len(coord)) if i == coord.index(coord[i])]) #Filter repeating coordinates the same kpoint
        coord_unique = np.reshape(coord_unique, (len(coord_unique),3)) #Reshape array 
        kptindex = list(set(i for i in range(kpt_number) for j in range(len(trans_df_plot['k index'])) if i == trans_df_plot['k index'][j])) #loop over all kpoints and append ones that match the k indices in trans_df_plot.
        lowerbnd = ['.' for i in range(len(kptindex))] #Empty data for band information in the plotting
        upperbnd = ['.' for i in range(len(kptindex))] #Empty data for band information in the plotting
        trans_df_plot = pd.DataFrame(data = {'xcoord': coord_unique[:,0], 'ycoord': coord_unique[:,1], 'zcoord': coord_unique[:,2], 'vband': lowerbnd, 'cband': upperbnd, 'prob': summedproblist, 'k index': kptindex}) #create new trans_df plot that satisfies the sum condition

    #Add information about every k-point, so that both dataframes can be combined to show all k-points alongside the transition-likely ones              
    k_df['vband'] = ['.' for i in range(kpt_number)]
    k_df['cband'] = ['.' for i in range(kpt_number)]
    k_df['prob'] = [0.02 for i in range(kpt_number)] #small number such that zero probability points are shown in the plot
    trans_df_plot = pd.concat([trans_df_plot,k_df], ignore_index=True) #combine both dataframes
    

 

    return trans_df_plot, max_prob #return modified trans_df

def Plot(trans_df_plot,max_prob): 
    """This function plots the modified dataframe from RunTransition.
    
    :param trans_df_plot: modified transmatrix dataframe for plotting
    :param max_prob: Maximum oscillator strength 
    :return: figure output for plotting
    """
    label = [] #list that will contain the displayed hovertext information
    for i in range(len(trans_df_plot['k index'])): #loop over length of the dataframe
        label.append(['Upper band: ' + str(trans_df_plot['cband'][i]), 'Lower band: ' + str(trans_df_plot['vband'][i]), 'K point: ' + str(trans_df_plot['k index'][i])]) #append information for the hovertext
    fig= go.Figure(go.Scatter3d(x=trans_df_plot['xcoord'], y=trans_df_plot['ycoord'], z=trans_df_plot['zcoord'], mode='markers', hovertext = label,hoverinfo="text",marker=dict(color="green"))) #plot in 3D
    fig.update_traces(marker=dict(color=trans_df_plot['prob'], size=trans_df_plot['prob']*100,colorscale='viridis', colorbar=dict(thickness=10, title = 'Oscillator Strength'), cmin = 0, cmax = max_prob)) #modify color/size of kpoints and add colorscale/colorbar.
    
    # define initial camera position
    x_0 = 2
    y_0 = 2
    z_0 = 1

    fig.update_layout(title='Probability',width=600,height=600,scene_camera_eye=dict(x=x_0, y=y_0, z=z_0))
    return fig
    # fig.show()

def Interactiveplot(trans_df, tr_emin, tr_emax, tr_deltaE, deltaE, option, kpt_number, k_df ):
    """This function creates an interactive plot for transition energies selected.
    
    :param trans_df: Dataframe with transmatrix values
    :param tr_emin: Minimum transition energy value to cut the transmatrix by
    :param tr_emax: Maximum transition energy value to cut the transmatrix by
    :param tr_deltaE: Energy step size (in units eV) of the interactive slider.
    :param deltaE: Determines the range of transition energies(E-deltaE, E+deltaE)  to be considered
    :param option: When set to 'yes', sums all of the transition probabilities for a given k-point
    :param cbandmin: Minimum conduction band value to cut the transmatrix by
    :param cbandmax: Maximum conduction band value to cut the transmatrix by
    :param vbandmin: Minimum valence band value to cut the transmatrix by
    :param vbandmax: Maximum valence band value to cut the transmatrix by
    :param kmin: Minimum k-point value to cut the transmatrix by
    :param kmax: Maximum k-point value to cut the transmatrix by
    :param kpt_number: Number of k-points
    :param k_df: Dataframe with positions of k-points
    
    """
    
    sfig = go.Figure()
    for e in np.arange(tr_emin,tr_emax,tr_deltaE):
        trans_df_plot, max_prob = RunTransition(trans_df, e,deltaE, option, kpt_number, k_df)
        #print(max(trans_df_plot['prob']), e)
        f=Plot(trans_df_plot,max_prob)
        sfig.add_trace(f.data[0])

    # Make 10th trace visible
    sfig.data[1].visible = True

    # Create and add slider
    steps = []
    for i in range(len(sfig.data)):
        step = dict(
            method="update",
            args=[{"visible": [False] * len(sfig.data)},
                  {"title": "Slider switched to energy: " + str(tr_emin+i*tr_deltaE) + "eV"}],  # layout attribute
        )
        step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
        steps.append(step)

    sliders = [dict(
        len = 0.9,
        active=10,
        currentvalue={"prefix": "Frequency: "},
        pad={"t": 50},
        steps=steps
    )]
    sfig.update_layout(
        sliders=sliders,
        height = 700,
        width = 700
    )
    sfig.show()


