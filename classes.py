# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none' # when exporting svg, keeps the text as text, not path
import os
import sys
import numpy as np
from hapi import *
db_begin('data')  # sets the folder for the hapi tables


class Observer():
    """class representing an observer as on spectracalc.com
    
    Parameters
        ----------
        unit : String
            'wav' : wavenumber [1/cm] or 'lam' : wavelength [nm]
        lower : Float
            
        upper: Float

        Returns
        -------
        None.
    
    """

    def __init__(self,
                 unit = 'wav', 
                 lower=3065.0,
                 upper=3069.0,
                 ):
        self.unit = unit
        self.line_list = []
        if self.unit == 'wav':
            self.lower_wav = lower
            self.upper_wav = upper
            self.lower_lam = Helpers.wav2lam(upper)
            self.upper_lam = Helpers.wav2lam(lower)
        elif self.unit == 'lam':
            self.lower_wav = Helpers.lam2wav(upper)
            self.upper_wav = Helpers.lam2wav(lower)
            self.lower_lam = lower
            self.upper_lam = upper
        else:
            print(f'{self.unit} type not defined. Please use "wav" for wavenumbers [1/cm] or "lam" for wavelength [nm].')
            return



class Gas():
    def __init__(self,
                 gas_name="H2O",
                 VMR=0.1):
        self.gas_name = gas_name
        self.VMR = VMR
        self.M = Helpers.hitran_molecule_number(
            self.gas_name)  # hitran molecule ID
        self.I = 1  # hitran Isotope number. Not yet taken into account, thus just 1


class Gas_Cell():
    """Class representing a gas cell as on spectracalc.com 
       After initialization, use the method .add_gas to add gasses to the cell."""
    def __init__(self,
                 name,
                 temperature,
                 pressure,
                 length,
                 no_gasses,
                 ):
        self.name = name
        self.temperature = temperature  # Kelvin
        self.pressure = pressure  # atm
        self.length = length  # cm
        #self.diluent = diluent # dictionary, like: {'air':0.8, 'O2' :0.2}
        self.diluent = {'air': 1,
                        'self':0}
        self.no_gasses = no_gasses  # number of gasses in the cell
        self.gasses = []
        # wavenumber
        self.nu = []
        self.coef = []
        self.absorp = []
        self.absorbance = []
        self.prop_const = []
        # wavelength
        self.lam = []
        self.coef_lam = []
        self.absorp_lam = []
        self.absorbance_lam = []
        self.prop_const_lam = []

    def add_gas(self, gas_name, VMR, *args):
        """
        Adds a new gas to the gas_cell

        Parameters
        ----------
        gas_name : String
            Has to fit the hitran gas names. see here: https://hitran.org/docs/molec-meta/
        VMR : Float
            Volume mixing ratio [0, 1]

        Returns
        -------
        None.

        """
        gas = Gas(gas_name=gas_name,
                  VMR=VMR)

        if len(self.gasses) < self.no_gasses:
            self.gasses.append(gas)
            self.diluent['air'] -= VMR
            self.diluent['self'] += VMR
            
            print(
                '####################################### \n' +
                f'Added {gas.gas_name} to the gas cell. \n*' +
                f'{len(self.gasses)}/{self.no_gasses} gasses are now in the cell: \n' +
                f'{[gas.gas_name for gas in self.gasses]}')
        else:
            print(
                f"Can't add more gasses to the cell. \n" +
                f'Gas cell already consists {self.no_gasses}: {[gas.gas_name for gas in self.gasses]}')


class Spectra():
    """Parent class representing containing all the information needed to plot one gas spectrum.
        containins following objects:
            - An observer (wavenumber range).
            - Gas cells (with specific environment conditions and gasses).
        
        contains following methods:
            - plot: quick way to visualize the results.
            - print: (__repr__) -> summarizes all the relevant parameter.
    """

    def __init__(self, name='my spectrum'):
        self.name = name
        self.observer = None
        self.gas_cells = []
        self.gas_cell_number = 0  # to count and adress the cells

    def add_gas_cell(self,
                     temperature=296,
                     pressure=1,
                     length=10,
                     no_gasses=1,
                     ):
        """
        Adds a gas cell to the spectra.
        Do not confuse with the gasses in one cell. Their absorption will be calculated combined.

        Parameters
        ----------
        temperature : float, optional
            [Kelvin]. The default is 296 K.
        pressure : float, optional
            [mbar]. The default is 1 atm.
        length : float, optional
            [cm]. The default is 10 cm.
        no_gasses : int, optional
            Number of gasses in the cell. The default is 1.


        Returns
        -------
        None.

        ! hapi yet only supports parallel gas cells. I.e. they will both be plotted,
        but the absorption will not be calculated compbined (i.e. in a row).
        """
        self.gas_cells.append(
            Gas_Cell(
                name = self.gas_cell_number,
                temperature = temperature,
                pressure = pressure,
                length = length,
                #diluent = diluent, # implemented such, that air - VMR. This way, selfbroadenin is correct.
                no_gasses = no_gasses))
        self.gas_cell_number += 1
        return None

    def __repr__(self, cell = 'all'):
        gas_cell_str = ""
        if cell == 'all': 
            cell = slice(0,len(self.gas_cells))
        else: cell = slice(cell, cell+1)
        for gas_cell in self.gas_cells[cell]:
            gas_cell_str += f'Gas cell {gas_cell.name}: \n'
            gas_cell_str += f'\t length: {gas_cell.length} cm |'
            gas_cell_str += f' temp: {gas_cell.temperature} K|'
            gas_cell_str += f'pressure: {gas_cell.pressure} atm|'
            gas_cell_str += f'gas matrix: {gas_cell.diluent} \n'
            gas_cell_str += '\t Gasses (VMR): \n'
            for gas in gas_cell.gasses:
                gas_cell_str += f'\t \t {gas.gas_name}: {gas.VMR} \n'
        
        if self.observer.unit == 'wav':
            observer_str = (f'\t lower: {self.observer.lower_wav} [1/cm] \n' +
                      f'\t upper: {self.observer.upper_wav} [1/cm] \n' )
        elif self.observer.unit == 'lam': 
            observer_str = (f'\t lower: {self.observer.lower_lam} [nm] \n' +
                            f'\t upper: {self.observer.upper_lam} [nm] \n' )

        return_str = (f'Summary of the spectum {self.name}: \n' +
                      observer_str + 
                      gas_cell_str
                      )
        return return_str

    def download(self, line_list=True, min_intensity = 5E-23, step = 0.002, HITRAN_units = False):
        """fetches data from HItran and calculates spectra based on the gas cells.
        Equivalent to the "calculate" button on spectracalc.

         Parameters
        ----------
        line_list : BOOL, optional
            if set True, a line list for the gasses in the cells will be downloaded.

        min_intensity : FLOAT, optional
            Minimum intensity for a gas-line to be shown. Default: 5E-23. This is only for the sticky-line plot, the absorption sepctra is calculated on all lines.
            
        step: float, optional
            Step of the x-Axis (wav, lam). Default: 0.002
        
        HITRAN_units: BOOL
            Whether to use Hitran Units (cgs) or SI. Default is False (i.e. SI units)

        Returns
        -------
        None.

        """

        for gas_cell in self.gas_cells:

            # fetch data into data folder
            # getHelp(fetch)
            for gas in gas_cell.gasses:
                try:
                    fetch(TableName=gas.gas_name,
                          M=gas.M,   # Hitran molecule number
                          I=gas.I,       # Isotopes number
                          numin=self.observer.lower_wav,
                          numax=self.observer.upper_wav
                          )
                except Exception as e:
                    print('no lines within the set interval.')

                if line_list:
                    _x, _y, = getStickXY(gas.gas_name)
                    _y[_y < min_intensity] = 0
                    self.observer.line_list.append({'nu': _x,
                                                    'y': _y,
                                                    'lam': np.flip(np.asarray([Helpers.wav2lam(wav) for wav in _x])),
                                                    'y_lam': np.flip(_y),
                                                    'label': gas.gas_name,
                                                    })

            # absorption coefficient per gas_cell (can contain multiple gasses)
            # getHelp(absorptionCoefficient_Voigt)
            try:
                print('to check: ----------------------- \n',
                      [(gas.M, gas.I, gas.VMR) for gas in gas_cell.gasses])
                
                gas_cell.nu, gas_cell.coef = absorptionCoefficient_Voigt(
                    Components=[(gas.M, gas.I, gas.VMR) for gas in gas_cell.gasses], 
                    SourceTables=[gas.gas_name for gas in gas_cell.gasses], 
                    HITRAN_units=False,
                    Environment={'T': gas_cell.temperature, 'p': gas_cell.pressure},
                    Diluent= gas_cell.diluent,
                    OmegaStep= step,
                    OmegaRange = [self.observer.lower_wav, self.observer.upper_wav])
                # if HITRAN_units:  # to do. see https://github.com/hapijs/hapi/issues/4206
                
                # wavelength conversion
                gas_cell.lam = np.flip(np.asarray([Helpers.wav2lam(wav) for wav in gas_cell.nu]))
                gas_cell.coef_lam = np.flip(gas_cell.coef)
            except IndexError:
                print(
                    f'Error: \n' +
                    f'Your gas cell {gas_cell.name} has not all gasses added. Currently {len(gas_cell.gasses)}/{gas_cell.no_gasses} gasses.')
                return

            # absorption spectrum
            #getHelp(absorptionSpectrum)
            _, gas_cell.absorp = absorptionSpectrum(
                gas_cell.nu, gas_cell.coef, Environment={
                    'l': gas_cell.length})
            gas_cell.absorp_lam = np.flip(gas_cell.absorp)
            
            # absorbance spectrum
            gas_cell.absorbance = [-np.log(1 - val) for val in gas_cell.absorp]
            gas_cell.absorbance_lam = np.flip(gas_cell.absorbance)
            
            # prop constants
            gas_cell.prop_const = gas_cell.gasses[0].VMR * gas_cell.length *1E4 / np.array(gas_cell.absorbance)
            gas_cell.prop_const_lam = gas_cell.gasses[0].VMR * gas_cell.length *1E4 / gas_cell.absorbance_lam
            
            
        return None
    def export(self, directory = "exports"):
        """ Method to export the spectral data. It will export the absorption data of your gas cells in the "directory" directorty.
        
         Parameters
        ----------
        dir : string, optional
            directory where to save the txt file. default = 'exports'


        Returns
        -------
        None.
        
        """
        # check for directory
        try: os.mkdir(os.path.join(os.getcwd(), directory)) 
        except FileExistsError as error: pass
        
        for gas_cell in self.gas_cells:
            filename = str(self.name) + '_gas_cell_' + str(gas_cell.name)
            
            with open(directory + os.sep + filename+'.txt', 'w') as f:
                # write header
                f.write(self.__repr__(cell = gas_cell.name))
                
                # write data (dependent on chosen unit)
                if self.observer.unit == 'wav':    
                    f.writelines('\n wavenumber [1/cm] \t absorption\n')
                    [f.writelines(str(nu) +'\t' +str(absorp) + '\n') for nu, absorp in zip(gas_cell.nu, gas_cell.absorp)]   
                elif self.observer.unit == 'lam':
                    f.writelines('\n wavelenth [nm] \t absorption\n')
                    [f.writelines(str(lam) +'\t' +str(absorp_lam) + '\n') for lam, absorp_lam in zip(gas_cell.lam, gas_cell.absorp_lam)] 
                else: print(f'{unit} unit not known. Please use "wav" for wavenumber [1/cm] or "lam" for wavelength [nm] as argument for the observer')
        print(f'Exported the absorption data to: {directory}')
    
    def plot(self, 
             ylim=None, 
             ylog=True, 
             export=False, 
             figsize = (12,6), 
             filetype = 'svg', 
             fontsize = 10, 
             color = None, 
             language = "English",
             absorbance = False,
             prop_const = False):
        """
        Simple plot function. You may adjust it to your needs.

        Parameters
        ----------
        ylim : List, optional
            You can set the ylimits, e.g. [0,1]. The default is automatic.

        ylog : BOOL, optional
            if True, the y axis of the line list plot will be logarithmic.
            
        export : BOOL, optional
            if True, the plot will be saved in exports/
        filetype : STRING, optional
            file extension of the save plot. For example 'svg', or 'pdf'
            
        color: STRING, optional
            String of the colorpalette to use or a single color. for example 'seaborn-v0_8-dark-palette', or 'orange'. See        
            https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html

        language : STRING, optional
            Either 'English' or 'German'. Sets the axis labels accordingly
        
        absorbance : BOOL, optional
            If set to True, instead of the absorption the aborbance will be plotted.
        
        prop_const :BOOL, optional
            If set to True, a 2nd y-axis with the proportionality factor FAC is drawn:  conc = FAC * Absorance
            Only yet works for one gas in a gas cell.
        -------
        None.

        """
        # set color palette
        if color == None: color = 'seaborn-v0_8-whitegrid'
        plt.style.use(color) 
        
        # parameter
        fontsize_subplot_title = fontsize +2
        fontsize_ax_label = fontsize
        fontsize_ticks = fontsize
        plt.rcParams['font.size']= fontsize
        plt.rcParams['axes.labelsize'] = fontsize
        plt.rcParams['xtick.labelsize'] = fontsize
        plt.rcParams['ytick.labelsize'] = fontsize

        


        
        '''
        # to change default colormap
        plt.rcParams["image.cmap"] = "Set1"
        # to change default color cycle
        plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.Set1.colors)
        '''
        #langauge stuff
        if language == 'German': 
            name_wavlen = 'Wellenlänge' 
            name_wavnum = 'Wellenzahl'
            if absorbance: name_absorp = 'Absorbanz'
            else: name_absorp = 'Absorption'
            name_cell = 'GasZelle'
            for_str = 'für'
        else:
            name_wavlen = 'wavelength' 
            name_wavnum = 'wavenumber'
            if absorbance: name_absorp = 'Absorbance'
            else: name_absorp = 'Absorption'
            name_cell = 'cell'
            for_str = 'for'
        
        
        # check if line list is available
        if self.observer.line_list:
            nrows = 2
        else:
            nrows = 1
        
        # create figure
        fig, axs = plt.subplots(nrows=nrows, ncols=1, figsize=figsize)
        #fig.suptitle(self.name)

        # set x range for all subplots
        if self.observer.unit == 'wav': custom_xlim = (self.observer.lower_wav, self.observer.upper_wav)
        else: custom_xlim = (self.observer.lower_lam, self.observer.upper_lam)
        plt.setp(axs, xlim=custom_xlim)
        
        # line list plot (if data is downdloaded)
        if self.observer.line_list:
            if self.observer.unit == 'wav': 
                [axs[0].plot(gas['nu'], gas['y'], label=gas['label']) 
                 for gas in self.observer.line_list]
                #axs[0].set_xlim([self.observer.lower_wav, self.observer.upper_wav])
            elif self.observer.unit == 'lam': 
                [axs[0].plot(gas['lam'], gas['y_lam'], label=gas['label'])
                 for gas in self.observer.line_list]
                #axs[0].set_xlim([self.observer.lower_lam, self.observer.upper_lam])
            
            # legend line list plot
            axs[0].legend(loc='upper right', shadow=True, fontsize=fontsize)
            axs[0].set_ylabel(
                r"intensity [$cm^{-1} * mol^{-1}*cm^2 $]",
                fontsize=fontsize_ax_label)
            axs[0].set_xticklabels(())
            #axs[0].set_title('Linelist', fontsize=fontsize_subplot_title)
            axs[0].tick_params(labelsize=fontsize_ticks)
            #axs[0].grid()
            if ylog:
                axs[0].set_yscale('log')

        
        # gas cell plot
        if self.observer.line_list:
            ax1 = axs[1]
        else:
            ax1 = axs
        if prop_const: 
            ax3 = ax1.twinx()
            ax3.grid(None)
            ax3.set_ylim(-10,430)
        for gas_cell in self.gas_cells:
            label_str = name_cell +' '+ str(gas_cell.name)+': '
            for gas in gas_cell.gasses:
                # convert gasVMR in nice units
                if gas.VMR < 1 and gas.VMR > 0.001: gas_VMR_str = f'{str(gas.VMR *100)} %'
                elif gas.VMR < 0.001 and gas.VMR > 10E-8: gas_VMR_str = f'{str(gas.VMR *1E6)} ppm'
                else: gas_VMR_str = str(gas.VMR)
                label_str += str(gas.gas_name) + ': ' + \
                    gas_VMR_str +' ' + for_str +' ' + str(gas_cell.length) + ' cm @' + str(gas_cell.pressure) +' atm'
            
            # handle what to plot
            if absorbance and self.observer.unit== 'lam': x, y = gas_cell.lam, gas_cell.absorbance_lam
            elif absorbance and self.observer.unit== 'wav': x, y = gas_cell.nu, gas_cell.absorbance
            elif not absorbance and self.observer.unit== 'wav': x, y = gas_cell.nu, gas_cell.absorp
            elif not absorbance and self.observer.unit== 'lam': x, y = gas_cell.lam, gas_cell.absorp_lam
            ax1.plot(x, y, label=label_str)
            if prop_const: 
                if self.observer.unit== 'lam': y2 = gas_cell.prop_const_lam
                elif self.observer.unit== 'wav': y2 = gas_cell.prop_const
                #y2[y2 > 400] = 0 # for low absorbance, we get otherwise really high values which destroy the plot
                ax3.plot(x,y2, '--', color='grey')    
                ax3.set_ylabel('absorption coefficient $k_{\\nu}$ [ppm*m]', fontsize = fontsize_ax_label)
                ax3.tick_params(axis='y', labelsize=fontsize_ticks)
       
        
        # legend gas cell plot
        #ax1.legend(loc='upper right', shadow=True, fontsize=fontsize, facecolor = 'white')
        ax1.legend(bbox_to_anchor = (0.5, 1.2), loc='lower center', shadow=True, fontsize=fontsize) # for specific positioning
        #https://stackoverflow.com/questions/44413020/how-to-specify-legend-position-in-matplotlib-in-graph-coordinates
        if self.observer.unit == 'wav':    ax1.set_xlabel(name_wavnum + ' [1/cm]', fontsize=fontsize_ax_label,
                                                         labelpad = int(fontsize*0.5))
        elif self.observer.unit == 'lam':  ax1.set_xlabel(name_wavlen + ' [nm]', fontsize=fontsize_ax_label,
                                                         labelpad = int(fontsize*0.5))
        ax1.set_ylabel(name_absorp, fontsize=fontsize_ax_label)
        #ax1.set_title('Gas cells', fontsize=fontsize_subplot_title)
        ax1.tick_params(labelsize=fontsize_ticks)
        #ax1.grid()

        # other unit on top
        if self.observer.unit == 'wav': 
            ax2 = ax1.secondary_xaxis(
                'top', functions=( Helpers.wav2lam, Helpers.lam2wav))
            ax2.set_xlabel(name_wavlen + ' [nm]', fontsize=fontsize_ax_label, 
                           labelpad= int(2 * fontsize*0.5))
        elif self.observer.unit == 'lam':
            ax2 = ax1.secondary_xaxis(
                'top', functions=( Helpers.lam2wav, Helpers.wav2lam))
            ax2.set_xlabel(name_wavnum + ' [1/cm]', fontsize=fontsize_ax_label, 
                           labelpad= int(2 * fontsize*0.5))
        fig.tight_layout()
        ax2.tick_params(labelsize=fontsize_ticks)
        if ylim:
            ax1.ylim(ylim)
            
        if export: 
            plt.rcParams['svg.fonttype'] = 'none'
            fig.patch.set_alpha(0.)
            plt.savefig(os.path.join('exports', self.name +'_plot.' + filetype))
            plt.savefig(os.path.join('exports', self.name +'_plot.pdf'), dpi=300)

        return ax1

class Helpers():
    @staticmethod
    def wav2lam(wn):
        # cm^{-1} to nm
        with np.errstate(divide='ignore'):
            lam = 1.0e7 / wn
        return lam

    @staticmethod
    def lam2wav(lam):
        # nm to cm^{-1}
        with np.errstate(divide='ignore'):
            wav = 1.0e7 / lam
        return wav

    @staticmethod
    def hitran_molecule_number(name):
        """given molecule name, it returns the hitran id"""
        hitran_molecule_dic = {
                            'H2O': 1,  # Water
                            'CO2': 2,  # Carbon Dioxide
                            'O3': 3,  # Ozone
                            'N2O': 4,  # Nitrous Oxide
                            'CO': 5,  # Carbon Monoxide
                            'CH4': 6,  # Methane
                            'O2': 7,  # Oxygen
                            'NO': 8,  # Nitric Oxide
                            'SO2': 9, # Sulfur Dioxide
                            'NO2': 10, # Nitrogen Dioxide
                            'NH3': 11, # Ammonia
                            'HNO3': 12, # Nitric Acid
                            'OH': 13, # Hydroxyl
                            'HF': 14, # Hydrogen Fluoride
                            'HCl': 15, # Hydrogen Chloride
                            'HBr': 16, # Hydrogen Bromide
                            'HI': 17, # Hydrogen Iodide
                            'ClO': 18, # Chlorine Monoxide
                            'OCS': 19, # Carbonyl Sulfide
                            'H2CO': 20, # Formaldehyde
                            'HOCl': 21, # Hypochlorous Acid
                            'N2': 22, # Nitrogen
                            'HCN': 23, # Hydrogen Cyanide
                            'CH3Cl': 24, # Methyl Chloride
                            'H2O2': 25, # Hydrogen Peroxide
                            'C2H2': 26, # Acetylene
                            'C2H6': 27, # Ethane
                            'PH3': 28, # Phosphine
                            'COF2': 29, # Carbonyl Fluoride
                            'SF6': 30, # Sulfur Hexafluoride
                            'H2S': 31, # Hydrogen Sulfide
                            'HCOOH': 32, # Formic Acid
                            'HO2': 33, # Hydroperoxyl
                            'O': 34, # Oxygen Atom
                            'ClONO2': 35, # Chlorine Nitrate
                            'NO+': 36, # Nitric Oxide Cation
                            'HOBr': 37, # Hypobromous Acid
                            'C2H4': 38, # Ethylene
                            'CH3OH': 39, # Methanol
                            'CH3Br': 40, # Methyl Bromide
                            'CH3CN': 41, # Acetonitrile
                            'CF4': 42, # PFC-14
                            'C4H2': 43, # Diacetylene
                            'HC3N': 44, # Cyanoacetylene
                            'H2': 45, # Hydrogen
                            'CS': 46, # Carbon Monosulfide
                            'SO3': 47, # Sulfur trioxide
                            'C2N2': 48, # Cyanogen
                            'COCl2': 49, # Phosgene
                            'SO': 50, # Sulfur Monoxide
                            'CH3F': 51, # Methyl fluoride
                            'GeH4': 52, # Germane
                            'CS2': 53, # Carbon disulfide
                            'CH3I': 54, # Methyl iodide
                            'NF3': 55, # Nitrogen trifluoride
                        }

        return hitran_molecule_dic[name]


class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
