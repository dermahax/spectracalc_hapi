# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
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
        self.no_gasses = no_gasses  # number of gasses in the cell
        self.gasses = []
        # wavenumber
        self.nu = []
        self.coef = []
        self.absorp = []
        # wavelength
        self.lam = []
        self.coef_lam = []
        self.absorp_lam = []

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
            gas_cell_str += f'pressure: {gas_cell.pressure} atm \n'
            gas_cell_str += '\t Gasses: \n'
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

    def download(self, line_list=True, min_intensity=5E-23):
        """fetches data from HItran and calculates spectra based on the gas cells.
        Equivalent to the "calculate" button on spectracalc.

         Parameters
        ----------
        line_list : BOOL, optional
            if set True, a line list for the gasses in the cells will be downloaded.

        min_intensity : FLOAT, optional
            Minimum intensity for a gas-line to be shown. Default: 5E-23.

        Returns
        -------
        None.

        """

        for gas_cell in self.gas_cells:

            # fetch data into data folder
            # getHelp(fetch)
            for gas in gas_cell.gasses:
                fetch(TableName=gas.gas_name,
                      M=gas.M,   # Hitran molecule number
                      I=gas.I,       # Isotopes number
                      numin=self.observer.lower_wav,
                      numax=self.observer.upper_wav
                      )

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
                    Components=[
                        (gas.M, gas.I, gas.VMR) for gas in gas_cell.gasses], SourceTables=[
                        gas.gas_name for gas in gas_cell.gasses], HITRAN_units=False, Environment={
                        'T': gas_cell.temperature, 'p': gas_cell.pressure})
                
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
    
    def plot(self, ylim=None, ylog=True, export=False):
        """
        Simple plot functio. You may adjust it to your needs.

        Parameters
        ----------
        ylim : List, optional
            You can set the ylimits, e.g. [0,1]. The default is automatic.

        ylog : BOOL, optional
            if True, the y axis of the line list plot will be logarithmic.

        Returns
        -------
        None.

        """
        # parameter
        fontsize_subplot_title = 10
        fontsize_ax_label = 8
        fontsize_ticks = 8

        # check if line list is available
        if self.observer.line_list:
            nrows = 2
        else:
            nrows = 1

        # create figure
        fig, axs = plt.subplots(nrows=nrows, ncols=1)
        fig.suptitle(self.name)

        # line list plot (if data is downdloaded)
        if self.observer.line_list:
            if self.observer.unit == 'wav': 
                [axs[0].plot(gas['nu'], gas['y'], label=gas['label'])
                 for gas in self.observer.line_list]
            elif self.observer.unit == 'lam': 
                [axs[0].plot(gas['lam'], gas['y_lam'], label=gas['label'])
                 for gas in self.observer.line_list]
            axs[0].legend(loc='upper right', shadow=True, fontsize='small')
            axs[0].set_ylabel(
                r"intensity [$cm^{-1} * mol^{-1}*cm^2 $]",
                fontsize=fontsize_ax_label)
            axs[0].set_xticklabels(())
            axs[0].set_title('Linelist', fontsize=fontsize_subplot_title)
            axs[0].tick_params(labelsize=fontsize_ticks)
            axs[0].grid()
            if ylog:
                axs[0].set_yscale('log')

        # gas cell plot
        if self.observer.line_list:
            ax1 = axs[1]
        else:
            ax1 = axs
        for gas_cell in self.gas_cells:
            label_str = ''
            for gas in gas_cell.gasses:
                label_str += str(gas.gas_name) + ': ' + \
                    str(gas.VMR / gas_cell.length) + ' *m'
            # wav | lam
            if self.observer.unit == 'wav':
                ax1.plot(gas_cell.nu, gas_cell.absorp,
                         label=label_str,  # names of all the gasses in the cell
                         )
            elif self.observer.unit == 'lam':
                ax1.plot(gas_cell.lam, gas_cell.absorp_lam,
                         label=label_str,  # names of all the gasses in the cell
                         )
        ax1.legend(loc='upper right', shadow=True, fontsize='small')
        if self.observer.unit == 'wav':    ax1.set_xlabel('wavenumber [1/cm]', fontsize=fontsize_ax_label)
        elif self.observer.unit == 'lam':  ax1.set_xlabel('wavelength [nm]', fontsize=fontsize_ax_label)
        ax1.set_ylabel('absorption', fontsize=fontsize_ax_label)
        ax1.set_title('Gas cells', fontsize=fontsize_subplot_title)
        ax1.tick_params(labelsize=fontsize_ticks)
        ax1.grid()

        # other unit on top
        if self.observer.unit == 'wav': 
            ax2 = ax1.secondary_xaxis(
                'top', functions=( Helpers.wav2lam, Helpers.lam2wav))
            ax2.set_xlabel('wavelength [nm]', fontsize=fontsize_ax_label)
        elif self.observer.unit == 'lam':
            ax2 = ax1.secondary_xaxis(
                'top', functions=( Helpers.lam2wav, Helpers.wav2lam))
            ax2.set_xlabel('wavenumber [1/cm]', fontsize=fontsize_ax_label)
        fig.tight_layout()
        ax2.tick_params(labelsize=fontsize_ticks)

        if ylim:
            ax1.ylim(ylim)
            
        if export: plt.savefig(self.name +'_plot.png')


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
            'NO': 8,  #	Nitric Oxide
         	'SO2': 9, #	Sulfur Dioxide
         	'NO2': 10, #	Nitrogen Dioxide
         	'NH3': 11, #	Ammonia
        }
        """
         these can be added if needed...
         12	HNO3	Nitric Acid
         13	OH	Hydroxyl
         14	HF	Hydrogen Fluoride
         15	HCl	Hydrogen Chloride
         16	HBr	Hydrogen Bromide
         17	HI	Hydrogen Iodide
         18	ClO	Chlorine Monoxide
         19	OCS	Carbonyl Sulfide
         20	H2CO	Formaldehyde
         21	HOCl	Hypochlorous Acid
         22	N2	Nitrogen
         23	HCN	Hydrogen Cyanide
         24	CH3Cl	Methyl Chloride
         25	H2O2	Hydrogen Peroxide
         26	C2H2	Acetylene
         27	C2H6	Ethane
         28	PH3	Phosphine
         29	COF2	Carbonyl Fluoride
         30	SF6	Sulfur Hexafluoride
         31	H2S	Hydrogen Sulfide
         32	HCOOH	Formic Acid
         33	HO2	Hydroperoxyl
         34	O	Oxygen Atom
         35	ClONO2	Chlorine Nitrate
         36	NO+	Nitric Oxide Cation
         37	HOBr	Hypobromous Acid
         38	C2H4	Ethylene
         39	CH3OH	Methanol
         40	CH3Br	Methyl Bromide
         41	CH3CN	Acetonitrile
         42	CF4	PFC-14
         43	C4H2	Diacetylene
         44	HC3N	Cyanoacetylene
         45	H2	Hydrogen
         46	CS	Carbon Monosulfide
         47	SO3	Sulfur trioxide
         48	C2N2	Cyanogen
         49	COCl2	Phosgene
         50	SO	Sulfur Monoxide
         51	CH3F	Methyl fluoride
         52	GeH4	Germane
         53	CS2	Carbon disulfide
         54	CH3I	Methyl iodide
         55	NF3	Nitrogen trifluoride
         """

        return hitran_molecule_dic[name]


class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
