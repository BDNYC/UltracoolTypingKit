import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table, Column
import os
import time

plt.rcParams['keymap.save'] = ''

'''
Ellianna Schwab, Kelle Cruz

TypeFinder.py contains typing_kit, which provides plots to qualitatively spectral type L-type brown dwarfs in NIR regime.
Typing methods are taken from Cruz et al. 2017 and offer comparison to NIR spectral standards.

To use input path to file_name in a string. A band-by-band grid of Cruz et al 2017 templates
will be overplotted with the target spectrum. To compare to the NIR spectral standards as well, 
key in the spectral type number, '0' - '8'. 
Cruz et al. 2017 templates are shown band-by-band followed by a comparison to the overall NIR spectrum of NIR 
spectral standards as defined in Table 9 of Cruz et al 2017.
'''

class Data:

    def __init__(self, file_name):

        '''
        This function splits the data into J-H-K bands and normalizes each band by its average.
        It also normalizes a separate copy of the overall spectra to Kirkpatrick '10 normalization specs.
        '''

        hdulist = fits.open(file_name)
        spectrum = hdulist[0]
        wavelength = spectrum.data[0]
        flux = spectrum.data[1]
        
        #Create J, H and K bands for spectrum

        wavelength_J = []
        flux_J = []
        wavelength_H = []
        flux_H = []
        wavelength_K = []
        flux_K = []
        for jj in range(len(wavelength)):
            if wavelength[jj] >= 0.87 and wavelength[jj] <= 1.39:
                wavelength_J.append(wavelength[jj])
                flux_J.append(flux[jj])
            elif wavelength[jj] >= 1.41 and wavelength[jj] <= 1.89:
                wavelength_H.append(wavelength[jj])
                flux_H.append(flux[jj])
            elif wavelength[jj] >= 1.91 and wavelength[jj] <= 2.39:
                wavelength_K.append(wavelength[jj])
                flux_K.append(flux[jj])


        #Normalize Each Band        

        wavelength_J = np.array(wavelength_J)
        flux_J = np.array(flux_J)
        wavelength_H = np.array(wavelength_H)
        flux_H = np.array(flux_H)
        wavelength_K = np.array(wavelength_K)
        flux_K = np.array(flux_K)
    
        flux_J = flux_J/np.mean(flux_J)
        flux_H = flux_H/np.mean(flux_H)
        flux_K = flux_K/np.mean(flux_K)

        #Normalize the Overall Spectra to Kirkpatrick Normalization
        norm_flux_array=[]
        for jj in range(len(wavelength)):
            if wavelength[jj] >= 1.28 and wavelength[jj] <= 1.39:
                norm_flux_array.append(flux[jj])
        norm_flux_array=np.array(norm_flux_array)

        wavelength=np.array(wavelength)
        flux=np.array(flux)
        norm_flux=flux/np.mean(norm_flux_array)

        self.wavelength_J = wavelength_J
        self.flux_J = flux_J
        self.wavelength_H = wavelength_H
        self.flux_H = flux_H
        self.wavelength_K = wavelength_K
        self.flux_K = flux_K
        self.wavelength = wavelength
        self.norm_flux = norm_flux


class Templates:

    def __init__(self):
        '''
        This function transfers all the template data from their ascii files into astropy tables.
        Templates are those compiled for the analysis in Cruz et al. 2017.
        '''

        #Beginning with the field templates. J-H-K bands are all diff lengths and so need their own table
        field_templates = defaultdict(dict)

        for band in ['J', 'H', 'K']:
            for spectype in range(9):
                field_templates[band][spectype] = Table.read('templates/Cruz2017_Templates.hdf5', \
                                                  path='field/L{}{}_template'.format(str(spectype), band))


        #Next with beta templates
        beta_templates = defaultdict(dict)

        for band in ['J', 'H', 'K']:
            for spectype in range(2):
                beta_templates[band][spectype]= Table.read('templates/Cruz2017_Templates.hdf5', \
                                                path='beta/L{}{}_template'.format(str(spectype), band))


        #Next with gamma templates
        gamma_templates = defaultdict(dict)

        for band in ['J', 'H', 'K']:
            for spectype in range(5):
                gamma_templates[band][spectype]= Table.read('templates/Cruz2017_Templates.hdf5', \
                                                path='gamma/L{}{}_template'.format(str(spectype), band))


        self.gamma_templates = gamma_templates
        self.beta_templates = beta_templates
        self.field_templates = field_templates


class MakePlot:

    def __init__(self, file_name, type_number, NIR_standards): #This defines often reused variables

        '''
        This initializes all the variables necessary to make the on_key_press plots.
        '''

        #Can define templates in each one! and can define all the pieces in each one and make multiple plot types :D

        self.wavelength_J = Data(file_name).wavelength_J
        self.flux_J = Data(file_name).flux_J
        self.wavelength_H = Data(file_name).wavelength_H
        self.flux_H = Data(file_name).flux_H
        self. wavelength_K = Data(file_name).wavelength_K
        self.flux_K = Data(file_name).flux_K
        self.wavelength = Data(file_name).wavelength
        self.norm_flux = Data(file_name).norm_flux

        #Define the other variables
        self.type_number = type_number
        self.NIR_standards = NIR_standards
        self.file_name = file_name


    def bracketed_plot(self, gravity_type): #This is for bracketed plots that are not first or last

        '''
        This creates the bracketed plots that pop up after key-press, as long as they are bracketed by L dwarfs.
        '''

        if gravity_type == 'f':
            templates = Templates().field_templates
        elif gravity_type == 'b':
            templates = Templates().beta_templates
        elif gravity_type == 'g':
            templates = Templates().gamma_templates

        #Create the Plots
        fig2, axes2 = plt.subplots(figsize=(15,45),
            nrows=3, ncols=4, sharex=False, sharey=False, 
            gridspec_kw={'width_ratios':[1,1,1,3]}, facecolor='white'
            )

        for jj, ii in zip([0,1,2], [self.type_number-1, self.type_number, self.type_number+1]):

                ##J Band
                axes2[jj, 0].plot((templates['J'][ii])['wavelength'], (templates['J'][ii])['flux'], c='red')
                axes2[jj, 0].fill_between((templates['J'][ii])['wavelength'], (templates['J'][ii])['upper'], \
                                          (templates['J'][ii])['lower'], color='#c6c6c6')
                axes2[jj, 0].plot(self.wavelength_J, self.flux_J, c='k')
                if gravity_type == 'f':
                    axes2[jj, 0].annotate('L{}'.format(ii), xy=(0.1, 0.8), xycoords='axes fraction', color='k')
                elif gravity_type == 'b':
                    axes2[jj, 0].annotate(r'L{}$\beta$'.format(ii), xy=(0.1, 0.8), xycoords='axes fraction', color='k')
                else:
                    axes2[jj, 0].annotate(r'L{}$\gamma$'.format(ii), xy=(0.1, 0.8), xycoords='axes fraction', color='k')
                axes2[jj, 0].set_ylim([-0.5, 2])
                axes2[jj, 0].axis('off')

                #H Band
                axes2[jj, 1].plot((templates['H'][ii])['wavelength'], (templates['H'][ii])['flux'], c='red')
                axes2[jj, 1].fill_between((templates['H'][ii])['wavelength'], (templates['H'][ii])['upper'], \
                                          (templates['H'][ii])['lower'], color='#c6c6c6')
                axes2[jj, 1].plot(self.wavelength_H, self.flux_H, c='k')
                axes2[jj, 1].axis('off')


                #K Band
                axes2[jj, 2].plot((templates['K'][ii])['wavelength'], (templates['K'][ii])['flux'], c='red')
                axes2[jj, 2].fill_between((templates['K'][ii])['wavelength'], (templates['K'][ii])['upper'], \
                                          (templates['K'][ii])['lower'], color='#c6c6c6')
                axes2[jj, 2].plot(self.wavelength_K, self.flux_K, c='k')
                axes2[jj, 2].axis('off')


                #All Together Now! This is where the spectral standard comes in.
                temp_hdulist = fits.open(self.NIR_standards[ii])
                temp_spectrum = temp_hdulist[0]
                temp_wavelength = temp_spectrum.data[0]
                temp_flux = temp_spectrum.data[1]

                temp_norm_flux_array=[]
                for kk in range(len(temp_wavelength)):
                    if temp_wavelength[kk] >= 1.28 and temp_wavelength[kk] <= 1.39:
                        temp_norm_flux_array.append(temp_flux[kk])
                temp_norm_flux_array=np.array(temp_norm_flux_array)

                temp_wavelength=np.array(temp_wavelength)
                temp_flux=np.array(temp_flux)
                temp_norm_flux=temp_flux/np.mean(temp_norm_flux_array)


                axes2[jj, 3].plot(temp_wavelength, temp_norm_flux, c='red')
                axes2[jj, 3].plot(self.wavelength, self.norm_flux, c='k')
                axes2[jj, 3].axis('off')
                
        save_name = 'new_types/' + self.file_name[12:-5] + r'_L{}{}'.format(self.type_number, gravity_type)

        def save_figure(event):
            if event.key == 'w':
                plt.savefig(save_name)
                save_types(self.file_name, self.type_number, gravity_type)

        fig2.canvas.mpl_connect('key_press_event', manipulate_figure)
        fig2.canvas.mpl_connect('key_press_event', save_figure)


    def first_bracketed_plot(self, gravity_type): 

        '''
        This creates the bracketed plots for L0 at all gravities that pop up after key-press.
        '''

        if gravity_type == 'f':
            templates = Templates().field_templates
        elif gravity_type == 'b':
            templates = Templates().beta_templates
        elif gravity_type == 'g':
            templates = Templates().gamma_templates

        #Create the Plots
        fig2, axes2 = plt.subplots(figsize=(15,45),
            nrows=3, ncols=4, sharex=False, sharey=False, 
            gridspec_kw={'width_ratios':[1,1,1,3]}, facecolor='white'
            )

        #Begin with the M9 template and then iterate over the rest
        axes2[0, 0].annotate('M9', xy=(0.1, 0.8), xycoords='axes fraction', color='k')
        axes2[0, 0].axis('off')
        axes2[0, 1].axis('off')
        axes2[0, 2].axis('off')
        
        temp_hdulist = fits.open("spectra/M9V_LHS2924.fits") #!!! For now M9 field template used for all gravities
        temp_spectrum = temp_hdulist[0]
        temp_wavelength = temp_spectrum.data[0]
        temp_flux = temp_spectrum.data[1]
        temp_norm_flux_array=[]
        for kk in range(len(temp_wavelength)):
            if temp_wavelength[kk] >= 1.28 and temp_wavelength[kk] <= 1.39:
                temp_norm_flux_array.append(temp_flux[kk])
        temp_norm_flux_array = np.array(temp_norm_flux_array)

        temp_wavelength=np.array(temp_wavelength)
        temp_flux=np.array(temp_flux)
        temp_norm_flux=temp_flux/np.mean(temp_norm_flux_array)
        axes2[0, 3].plot(temp_wavelength, temp_norm_flux, c='red')
        axes2[0, 3].plot(self.wavelength, self.norm_flux, c='k')
        axes2[0, 3].axis('off')


        for jj, ii in zip([1,2], [self.type_number, self.type_number+1]):

                ##J Band
                axes2[jj, 0].plot((templates['J'][ii])['wavelength'], (templates['J'][ii])['flux'], c='red')
                axes2[jj, 0].fill_between((templates['J'][ii])['wavelength'], (templates['J'][ii])['upper'], \
                                          (templates['J'][ii])['lower'], color='#c6c6c6')
                axes2[jj, 0].plot(self.wavelength_J, self.flux_J, c='k')
                if gravity_type == 'f':
                    axes2[jj, 0].annotate('L{}'.format(ii), xy=(0.1, 0.8), xycoords='axes fraction', color='k')
                elif gravity_type == 'b':
                    axes2[jj, 0].annotate(r'L{}$\beta$'.format(ii), xy=(0.1, 0.8), xycoords='axes fraction', color='k')
                else:
                    axes2[jj, 0].annotate(r'L{}$\gamma$'.format(ii), xy=(0.1, 0.8), xycoords='axes fraction', color='k')                
                axes2[jj, 0].set_ylim([-0.5, 2])
                axes2[jj, 0].axis('off')

                #H Band
                axes2[jj, 1].plot((templates['H'][ii])['wavelength'], (templates['H'][ii])['flux'], c='red')
                axes2[jj, 1].fill_between((templates['H'][ii])['wavelength'], (templates['H'][ii])['upper'], \
                                          (templates['H'][ii])['lower'], color='#c6c6c6') 
                axes2[jj, 1].plot(self.wavelength_H, self.flux_H, c='k')
                axes2[jj, 1].axis('off')


                #K Band
                axes2[jj, 2].plot((templates['K'][ii])['wavelength'], (templates['K'][ii])['flux'], c='red')
                axes2[jj, 2].fill_between((templates['K'][ii])['wavelength'], (templates['K'][ii])['upper'], \
                                          (templates['K'][ii])['lower'], color='#c6c6c6')
                axes2[jj, 2].plot(self.wavelength_K, self.flux_K, c='k')
                axes2[jj, 2].axis('off')


                #All Together Now! This is where the spectral standard comes in.
                temp_hdulist = fits.open(self.NIR_standards[ii])
                temp_spectrum = temp_hdulist[0]
                temp_wavelength = temp_spectrum.data[0]
                temp_flux = temp_spectrum.data[1]

                temp_norm_flux_array=[]
                for kk in range(len(temp_wavelength)):
                    if temp_wavelength[kk] >= 1.28 and temp_wavelength[kk] <= 1.39:
                        temp_norm_flux_array.append(temp_flux[kk])
                temp_norm_flux_array=np.array(temp_norm_flux_array)

                temp_wavelength=np.array(temp_wavelength)
                temp_flux=np.array(temp_flux)
                temp_norm_flux=temp_flux/np.mean(temp_norm_flux_array)


                axes2[jj, 3].plot(temp_wavelength, temp_norm_flux, c='red')
                axes2[jj, 3].plot(self.wavelength, self.norm_flux, c='k')
                axes2[jj, 3].axis('off')
                
        save_name = 'new_types/' + self.file_name[12:-5] + r'_L{}{}'.format(self.type_number, gravity_type)

        def save_figure(event):
            if event.key == 'w':
                plt.savefig(save_name)
                save_types(self.file_name, self.type_number, gravity_type)

        fig2.canvas.mpl_connect('key_press_event', manipulate_figure)
        fig2.canvas.mpl_connect('key_press_event', save_figure)


    def last_bracketed_plot(self, gravity_type):

        '''
        This creates the bracketed plots at all gravities that are final in their gravity and that pop up after key-press.
        '''

        if gravity_type == 'f':
            templates = Templates().field_templates
        elif gravity_type == 'b':
            templates = Templates().beta_templates
        elif gravity_type == 'g':
            templates = Templates().gamma_templates

    	#Create the Plots
        fig2, axes2 = plt.subplots(figsize=(15,45),
            nrows=3, ncols=4, sharex=False, sharey=False, 
            gridspec_kw={'width_ratios':[1,1,1,3]}, facecolor='white'
            )

        for jj, ii in zip([0,1], [self.type_number-1, self.type_number]):


                ##J Band
                axes2[jj, 0].plot((templates['J'][ii])['wavelength'], (templates['J'][ii])['flux'], c='red')
                axes2[jj, 0].fill_between((templates['J'][ii])['wavelength'], (templates['J'][ii])['upper'], \
                                          (templates['J'][ii])['lower'], color='#c6c6c6')
                axes2[jj, 0].plot(self.wavelength_J, self.flux_J, c='k')
                if gravity_type == 'f':
                    axes2[jj, 0].annotate('L{}'.format(ii), xy=(0.1, 0.8), xycoords='axes fraction', color='k')
                elif gravity_type == 'b':
                    axes2[jj, 0].annotate(r'L{}$\beta$'.format(ii), xy=(0.1, 0.8), xycoords='axes fraction', color='k')
                else:
                    axes2[jj, 0].annotate(r'L{}$\gamma$'.format(ii), xy=(0.1, 0.8), xycoords='axes fraction', color='k')                
                axes2[jj, 0].set_ylim([-0.5, 2])
                axes2[jj, 0].axis('off')

                #H Band
                axes2[jj, 1].plot((templates['H'][ii])['wavelength'], (templates['H'][ii])['flux'], c='red')
                axes2[jj, 1].fill_between((templates['H'][ii])['wavelength'], (templates['H'][ii])['upper'], \
                                          (templates['H'][ii])['lower'], color='#c6c6c6') 
                axes2[jj, 1].plot(self.wavelength_H, self.flux_H, c='k')
                axes2[jj, 1].axis('off')


                #K Band
                axes2[jj, 2].plot((templates['K'][ii])['wavelength'], (templates['K'][ii])['flux'], c='red')
                axes2[jj, 2].fill_between((templates['K'][ii])['wavelength'], (templates['K'][ii])['upper'], \
                                          (templates['K'][ii])['lower'], color='#c6c6c6')
                axes2[jj, 2].plot(self.wavelength_K, self.flux_K, c='k')
                axes2[jj, 2].axis('off')


                #All Together Now! This is where the spectral standard comes in.
                temp_hdulist = fits.open(self.NIR_standards[jj])
                temp_spectrum = temp_hdulist[0]
                temp_wavelength = temp_spectrum.data[0]
                temp_flux = temp_spectrum.data[1]

                temp_norm_flux_array=[]
                for kk in range(len(temp_wavelength)):
                    if temp_wavelength[kk] >= 1.28 and temp_wavelength[kk] <= 1.39:
                        temp_norm_flux_array.append(temp_flux[kk])
                temp_norm_flux_array=np.array(temp_norm_flux_array)

                temp_wavelength=np.array(temp_wavelength)
                temp_flux=np.array(temp_flux)
                temp_norm_flux=temp_flux/np.mean(temp_norm_flux_array)


                axes2[jj, 3].plot(temp_wavelength, temp_norm_flux, c='red')
                axes2[jj, 3].plot(self.wavelength, self.norm_flux, c='k')
                axes2[jj, 3].axis('off')
                
        #Finish with T0 template
        axes2[2, 0].annotate('T0', xy=(0.1, 0.8), xycoords='axes fraction', color='k')
        axes2[2, 0].axis('off')
        axes2[2, 1].axis('off')
        axes2[2, 2].axis('off')
        
        temp_hdulist = fits.open("spectra/T0_2M0423_T0.fits")
        temp_spectrum = temp_hdulist[0]
        temp_wavelength = temp_spectrum.data[0]
        temp_flux = temp_spectrum.data[1]
        temp_norm_flux_array=[]
        for kk in range(len(temp_wavelength)):
            if temp_wavelength[kk] >= 1.28 and temp_wavelength[kk] <= 1.39:
                temp_norm_flux_array.append(temp_flux[kk])
        temp_norm_flux_array=np.array(temp_norm_flux_array) 
        temp_wavelength=np.array(temp_wavelength)
        temp_flux=np.array(temp_flux)
        temp_norm_flux=temp_flux/np.mean(temp_norm_flux_array)
        
        axes2[2, 3].plot(temp_wavelength, temp_norm_flux, c='red')
        axes2[2, 3].plot(self.wavelength, self.norm_flux, c='k')
        axes2[2, 3].axis('off')

        save_name = 'new_types/' + self.file_name[12:-5] + r'_L{}{}'.format(self.type_number, gravity_type)

        def save_figure(event):
            if event.key == 'w':
                plt.savefig(save_name)
                save_types(self.file_name, self.type_number, gravity_type)

        fig2.canvas.mpl_connect('key_press_event', manipulate_figure)
        fig2.canvas.mpl_connect('key_press_event', save_figure)


def save_types(file_name, type_number, gravity_type):

    '''
    saves output csv file with filename in column 1, and spectral type in column 2
    '''

    now_date = time.strftime("%d-%m-%Y")

    if os.path.isfile(r'new_types/{}_SpecTypes.csv'.format(now_date)) == True :
        f = open(r'new_types/{}_SpecTypes.csv'.format(now_date), 'a')
        f.write(r'{}, L{}{}'.format(file_name[12:], type_number, gravity_type))
        f.write('\n')
        f.close()

    else :
        f = open(r'new_types/{}_SpecTypes.csv'.format(now_date), 'w')
        f.write(r'{}, L{}{}'.format(file_name[12:], type_number, gravity_type))
        f.write('\n')        
        f.close()


def manipulate_figure(event):

    '''
    closes selected plot window when q is pressed on the keyboard
    '''

    if event.key == 'q':
        plt.close(event.canvas.figure)

### Function definition
def typing_kit(file_name, make_templates=False) :

    '''
    Ellianna Schwab, Kelle Cruz

    typing_kit provides plots to qualitatively spectral type L-type brown dwarfs in NIR regime.
    Typing methods are taken from Cruz et al. 2017 and offer comparison to NIR spectral standards.

    To use input path to file_name in a string. A band-by-band grid of Cruz et al 2017 templates
    will be overplotted with the target spectrum. To compare to the NIR spectral standards as well, 
    key in the spectral type number, '0' - '8'. 
    Cruz et al. 2017 templates are shown band-by-band followed by a comparison to the overall NIR spectrum of NIR 
    spectral standards as defined in Table 9 of Cruz et al 2017.
    '''


    ##===============
    #If Necessary, Convert Ascii templates to HDF5
    ##===============

    if make_templates == True:

        #First for field templates
        for band in ['J', 'H', 'K']:
            for spectype in range(9):
                band_data = ascii.read("templates/L{}{}_f.txt".format(spectype, band))
                T = Table([band_data['col1'], band_data['col2'], band_data['col4'], band_data['col5']], \
                          names = ('wavelength', 'flux', 'upper', 'lower'))
                T.write('templates/Cruz2017_Templates.hdf5', path='field/L{}{}_template'.format(str(spectype), band), append=True)

        #Next with beta templates
        for band in ['J', 'H', 'K']:
            for spectype in range(2):
                band_data = ascii.read("templates/L{}{}_b.txt".format(spectype, band))
                T = Table([band_data['col1'], band_data['col2'], band_data['col4'], band_data['col5']], \
                          names = ('wavelength', 'flux', 'upper', 'lower'))
                T.write('templates/Cruz2017_Templates.hdf5', path='beta/L{}{}_template'.format(str(spectype), band), append=True)

        #Next with gamma templates
        for band in ['J', 'H', 'K']:
            for spectype in range(5):
                band_data = ascii.read("templates/L{}{}_g.txt".format(spectype, band))
                T = Table([band_data['col1'], band_data['col2'], band_data['col4'], band_data['col5']], \
                          names = ('wavelength', 'flux', 'upper', 'lower'))
                T.write('templates/Cruz2017_Templates.hdf5', path='gamma/L{}{}_template'.format(str(spectype), band), append=True)



    ##===============
    #This defines the template arrays
    ##===============


    gamma_templates = Templates().gamma_templates
    beta_templates = Templates().beta_templates
    field_templates = Templates().field_templates


    ##===============
    #This defines the input spectrum data arrays
    ##===============


    wavelength_J = Data(file_name).wavelength_J
    flux_J = Data(file_name).flux_J
    wavelength_H = Data(file_name).wavelength_H
    flux_H = Data(file_name).flux_H
    wavelength_K = Data(file_name).wavelength_K
    flux_K = Data(file_name).flux_K


    ##===============
    #This defines the initial grid plot
    ##===============
    

    fig1, axes1 = plt.subplots(9, 9, figsize=(16,11), facecolor='white')


    #Create the plots for alpha gravity

    for ii in range(9):

		##J Band
    	axes1[ii, 0].plot((field_templates['J'][ii])['wavelength'], (field_templates['J'][ii])['flux'], c='red')
    	axes1[ii, 0].fill_between((field_templates['J'][ii])['wavelength'], (field_templates['J'][ii])['upper'], \
                                  (field_templates['J'][ii])['lower'], color='#c6c6c6')
    	axes1[ii, 0].plot(wavelength_J, flux_J, c='k') #target object
    	axes1[ii, 0].annotate('L{}'.format(ii), xy=(0.1, 0.9), xycoords='axes fraction', color='k')
    	axes1[ii, 0].axis('off')

		#H Band
    	axes1[ii, 1].plot((field_templates['H'][ii])['wavelength'], (field_templates['H'][ii])['flux'], c='red') 
    	axes1[ii, 1].fill_between((field_templates['H'][ii])['wavelength'], (field_templates['H'][ii])['upper'], \
                                  (field_templates['H'][ii])['lower'], color='#c6c6c6') 
    	axes1[ii, 1].plot(wavelength_H, flux_H, c='k') #target object
    	axes1[ii, 1].axis('off')

		#K Band
    	axes1[ii, 2].plot((field_templates['K'][ii])['wavelength'], (field_templates['K'][ii])['flux'], c='red')
    	axes1[ii, 2].fill_between((field_templates['K'][ii])['wavelength'], (field_templates['K'][ii])['upper'], \
                                  (field_templates['K'][ii])['lower'], color='#c6c6c6')
    	axes1[ii, 2].plot(wavelength_K, flux_K, c='k') #target object
    	axes1[ii, 2].axis('off')
            
            
    #Create the plots for beta gravity
            
    for ii in range(2):

		##J Band
    	axes1[ii, 3].plot((beta_templates['J'][ii])['wavelength'], (beta_templates['J'][ii])['flux'], c='red')
    	axes1[ii, 3].fill_between((beta_templates['J'][ii])['wavelength'], (beta_templates['J'][ii])['upper'], \
		                          (beta_templates['J'][ii])['lower'], color='#c6c6c6')
    	axes1[ii, 3].plot(wavelength_J, flux_J, c='k')
    	axes1[ii, 3].annotate(r'L{}$\beta$'.format(ii), xy=(0.1, 0.9), xycoords='axes fraction', color='k')
    	axes1[ii, 3].axis('off')

		#H Band
    	axes1[ii, 4].plot((beta_templates['H'][ii])['wavelength'], (beta_templates['H'][ii])['flux'], c='red') 
    	axes1[ii, 4].fill_between((beta_templates['H'][ii])['wavelength'], (beta_templates['H'][ii])['upper'], \
		                          (beta_templates['H'][ii])['lower'], color='#c6c6c6') 
    	axes1[ii, 4].plot(wavelength_H, flux_H, c='k')
    	axes1[ii, 4].axis('off')

		#K Band
    	axes1[ii, 5].plot((beta_templates['K'][ii])['wavelength'], (beta_templates['K'][ii])['flux'], c='red')
    	axes1[ii, 5].fill_between((beta_templates['K'][ii])['wavelength'], (beta_templates['K'][ii])['upper'], \
		                          (beta_templates['K'][ii])['lower'], color='#c6c6c6')
    	axes1[ii, 5].plot(wavelength_K, flux_K, c='k')
    	axes1[ii, 5].axis('off')
            
    for ii in range(2,9):
        for jj in [3,4,5]:
            axes1[ii, jj].axis('off')

    #Create the plots for gamma gravity
            
    for ii in range(5):

		##J Band
    	axes1[ii, 6].plot((gamma_templates['J'][ii])['wavelength'], (gamma_templates['J'][ii])['flux'], c='red')
    	axes1[ii, 6].fill_between((gamma_templates['J'][ii])['wavelength'], (gamma_templates['J'][ii])['upper'], \
		                          (gamma_templates['J'][ii])['lower'], color='#c6c6c6')
    	axes1[ii, 6].plot(wavelength_J, flux_J, c='k')
    	axes1[ii, 6].annotate(r'L{}$\gamma$'.format(ii), xy=(0.1, 0.9), xycoords='axes fraction', color='k')
    	axes1[ii, 6].axis('off')

		#H Band
    	axes1[ii, 7].plot((gamma_templates['H'][ii])['wavelength'], (gamma_templates['H'][ii])['flux'], c='red')
    	axes1[ii, 7].fill_between((gamma_templates['H'][ii])['wavelength'], (gamma_templates['H'][ii])['upper'], \
		                          (gamma_templates['H'][ii])['lower'], color='#c6c6c6') 
    	axes1[ii, 7].plot(wavelength_H, flux_H, c='k')
    	axes1[ii, 7].axis('off')

		#K Band
    	axes1[ii, 8].plot((gamma_templates['K'][ii])['wavelength'], (gamma_templates['K'][ii])['flux'], c='red')
    	axes1[ii, 8].fill_between((gamma_templates['K'][ii])['wavelength'], (gamma_templates['K'][ii])['upper'], \
		                          (gamma_templates['K'][ii])['lower'], color='#c6c6c6')
    	axes1[ii, 8].plot(wavelength_K, flux_K, c='k')
    	axes1[ii, 8].axis('off')
            
    for ii in range(5,9):
        for jj in [6,7,8]:
            axes1[ii, jj].axis('off')            
            
    fig1.tight_layout()
    
    
    ##===============
    #This now defines what happens after key_press
    ##===============
    
    
    def on_key_press(event):

        #List the NIR_standards, so they can be called in the plots
        NIR_standards = ["spectra/nir/U20165_0345+2540_new_w.fits", "spectra/nir/spex_prism_2130-0845_080713.fits", \
                         "spectra/nir/U10244_0408-1450.fits", "spectra/nir/u11291_1506+1321_050323.fits", \
                         "spectra/nir/U12101_2158-1550_davy.fits", "spectra/nir/spex_prism_2137+0808_U20909.fits", \
                         "spectra/nir/U10880.fits", "spectra/nir/u10721_050323.fits", "spectra/nir/2M1632.fits"]
        
        beta_NIR_standards = ["spectra/nir/spex_prism_1552+2948_080730.fits", "spectra/nir/U20114_0227-1624.fits"]

        gamma_NIR_standards = ["spectra/nir/spex_prism_0141-4633_040905.fits", "spectra/nir/U10381_0518-2756.fits", \
                               "spectra/nir/U10397.fits", "spectra/nir/spex_prism_2208+2921_U40004.fits", \
                               "spectra/nir/U20198_2M0501-0010_dagny.fits"]
        
        ##==============
        #First this begins with field gravity options, uses integer keys
        ##==============
        
        # If the key is in an easily bracketed part, or 1,2,3,4,5,6,7
        if event.key in '1234567':

            type_number = int(event.key)

            field_vars = MakePlot(file_name, 
            	                  type_number, 
            	                  NIR_standards)

            field_vars.bracketed_plot('f')


        # For objects bracketed by M dwarf above
        elif event.key == '0':
            type_number = int(event.key)
           
            field_vars = MakePlot(file_name, 
            	                  type_number, 
            	                  NIR_standards)

            field_vars.first_bracketed_plot('f')


        # For objects bracketed by T dwarf below            
        elif event.key == '8':
            type_number = int(event.key)

            #List the NIR_standards, so they can be called
            NIR_standards = ["spectra/nir/u10721_050323.fits", "spectra/nir/2M1632.fits"] 

            field_vars = MakePlot(file_name, 
            	                  type_number, 
            	                  NIR_standards)

            field_vars.last_bracketed_plot('f')

        
        ##==============
        #For beta gravity options, uses ctrl+int
        ##==============
        
        elif event.key == 'ctrl+0':
            type_number = int(event.key[-1])
           
            beta_vars = MakePlot(file_name, 
            	                 type_number, 
            	                 beta_NIR_standards)

            beta_vars.first_bracketed_plot('b')


        elif event.key == 'ctrl+1':
            type_number = int(event.key[-1])

            beta_vars = MakePlot(file_name, 
            	                 type_number, 
            	                 beta_NIR_standards)

            beta_vars.last_bracketed_plot('b')

            
        ##==============
        #For gamma gravity options, uses alt+int
        ##==============
        
        elif event.key == 'alt+0':
            type_number = int(event.key[-1])
           
            gamma_vars = MakePlot(file_name, 
            	                  type_number, 
            	                  gamma_NIR_standards)

            gamma_vars.first_bracketed_plot('g')

                    
        elif event.key == 'alt+1' or event.key == 'alt+2' or event.key =='alt+3':
            type_number = int(event.key[-1])
            
            gamma_vars = MakePlot(file_name, 
            	                  type_number, 
            	                  gamma_NIR_standards)

            gamma_vars.bracketed_plot('g')

                    
        elif event.key == 'alt+4':
            type_number = int(event.key[-1])

            gamma_vars = MakePlot(file_name, 
            	                  type_number, 
            	                  gamma_NIR_standards)

            gamma_vars.last_bracketed_plot('g')


        #Raise an error message if any other key is pressed
        elif event.key in 'abcdefghijklmnoprstuvxyz`-=[];,.':
            raise Exception('Invalid Key. Try any int in range(10) for field, ctrl + int in range(1) for beta, or alt + int in range(4) for gamma.')    
            

    ##==============
    #Close the key_press loop and draw the new figures
    ##==============

    fig1.canvas.draw()
    
    fig1.canvas.mpl_disconnect(fig1.canvas.manager.key_press_handler_id)
    fig1.canvas.mpl_connect('key_press_event', on_key_press) #This leaves it open to a second try

    ##==============
    #Close the plot via the keyboard
    ##==============

    fig1.canvas.mpl_connect('key_press_event', manipulate_figure)