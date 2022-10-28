# -*- coding: utf-8 -*-

import os,sys
from hapi import *
db_begin ('data') # sets the folder for the hapi tables
import matplotlib.pyplot as plt


class Observer():
    """class representing an observer as on spectracalc.com"""
    def __init__(self, 
                 nu_min = 3065.0,
                 nu_max = 3069.0,
                 ):
        self.nu_min = nu_min
        self.nu_max = nu_max


class Gas():
    def __init__(self,
                 gas_name = "H2O",
                 VMR = 0.1):
        self.gas_name = gas_name
        self.VMR = VMR
        self.M = Helpers.hitran_molecule_number(self.gas_name) # hitran molecule ID
        self.I = 1 # hitran Isotope number. Not yet taken into account, thus just 1
    

class Gas_Cell():
    """Class representing a gas cell as on spectracalc.com """
    def __init__(self, 
                 name = '',
                 temperature = 296, 
                 pressure = 1,
                 length = 10,
                 no_gasses = 1
                 ):
        self.name = name
        self.temperature = temperature # Kelvin
        self.pressure = pressure # atm
        self.length = length # cm
        self.no_gasses = no_gasses # number of gasses in the cell
        self.gasses = []
        self.nu = []
        self.coef = []
        self.absorp = []
    
    def add_gas(self, gas_name, VMR):
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
        gas = Gas(gas_name = gas_name, 
                  VMR = VMR)
        
        if len(self.gasses)< self.no_gasses:
            self.gasses.append(gas)
            print('####################################### \n' +
                  f'Added {gas.gas_name} to the gas cell. \n*' +
                  f'{len(self.gasses)}/{self.no_gasses} gasses are now in the cell: \n'+ 
                  f'{[gas.gas_name for gas in self.gasses]}')
        else: 
            print(f"Can't add more gasses to the cell. \n" +
                  f'Gas cell already consists {self.no_gasses}: {[gas.gas_name for gas in self.gasses]}')

                   

class Spectra():
    """Parent class, containing an observer and gas_cells"""
    def __init__(self, name = 'my spectrum'):
        self.name = name
        self.observer = None
        self.gas_cells = []
        self.gas_cell_number = 0 # to count and adress the cells
        
    def add_gas_cell(self, 
                     temperature = 296, 
                     pressure = 1,
                     length = 10,
                     no_gasses = 1,
                     ):
        """
        Adds a gas cell to the spectra.
        ! hapi yet only does parallel gas cells. I.e. they will both be plotted,
        but the absorption will not be calculated compbined. 
        (Do not confuse with the gasses in one cell. Their absorption will be calculated combined.)

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

        """
        self.gas_cells.append(Gas_Cell(self.gas_cell_number, temperature, pressure, 
                                       length, no_gasses))
        self.gas_cell_number += 1
        return None
     
    def __repr__(self):
        gas_cell_str =''
        for gas_cell in self.gas_cells:
            gas_cell_str += f'Gas cell {gas_cell.name}: \n'
            gas_cell_str += f'\t length: {gas_cell.length} cm |'
            gas_cell_str += f' temp: {gas_cell.temperature} K|'
            gas_cell_str += f'pressure: {gas_cell.pressure} atm \n'
            gas_cell_str += '\t Gasses: \n'
            for gas in gas_cell.gasses:
                gas_cell_str += f'\t \t {gas.gas_name}: {gas.VMR} \n'
            
            
        return_str = ('##################################### \n' +
                      f'Summary of the spectum {self.name}: \n' +
                      f'\t nuMin: {self.observer.nu_min} \n' +
                      f'\t nuMax: {self.observer.nu_max} \n' +
                      gas_cell_str +
                      '##################################### \n'
                      )
        return return_str
          
    def download(self):
        """fetches data, calculates spectra and plots them.
        Equivalent to the "calculate" button on spectracalc."""
        

        for gas_cell in self.gas_cells:
            
            # fetch data into data folder
            # getHelp(fetch)
            for gas in gas_cell.gasses:
                fetch( TableName = gas.gas_name,
                       M = gas.M,   # Hitran molecule number
                       I = gas.I,       # Isotopes number
                       numin = self.observer.nu_min,
                       numax = self.observer.nu_max
                       )
        
            # absorption coefficient per gas_cell (can contain multiple gasses)
            #getHelp(absorptionCoefficient_Voigt)
            try:
                print('to check: ----------------------- \n', [(gas.M , gas.I , gas.VMR) for gas in gas_cell.gasses])
                gas_cell.nu, gas_cell.coef =  absorptionCoefficient_Voigt( Components = [(gas.M , gas.I , gas.VMR) for gas in gas_cell.gasses],
                                                    SourceTables = [gas.gas_name for gas in gas_cell.gasses] ,
                                                    HITRAN_units = False, 
                                                    Environment ={'T':gas_cell.temperature, 'p':gas_cell.pressure})
            except IndexError: 
                print(f'Error: \n' +
                      f'Your gas cell {gas_cell.name} has not all gasses added. Currently {len(gas_cell.gasses)}/{gas_cell.no_gasses} gasses.')
                return
            
            # absorption spectrum
            getHelp(absorptionSpectrum)
            _, gas_cell.absorp = absorptionSpectrum(gas_cell.nu, gas_cell.coef, Environment = {'l' : gas_cell.length})
        return None
            

    def plot(self, ylim = None):
        """
        Simple plot functio. You may adjust it to your needs.

        Parameters
        ----------
        ylim : List, optional
            You can set the ylimits, e.g. [0,1]. The default is automatic.

        Returns
        -------
        None.

        """
        f = plt.figure()
        fig, ax1 = plt.subplots()
        for gas_cell in self.gas_cells:
            label_str = ''
            for gas in gas_cell.gasses:
                label_str += str(gas.gas_name) + ': ' + str(gas.VMR / gas_cell.length) + ' *m'
            ax1.plot(gas_cell.nu, gas_cell.absorp, 
                     label= label_str,  # names of all the gasses in the cell
                     )
        ax1.legend(loc='upper right', shadow=True, fontsize='small')
        ax1.set_xlabel('wavenumber [1/cm]')
        ax1.set_ylabel('absorption')
        ax1.grid()
        ax1.set_title(self.name)
        
        # wavelength on top
        ax2 = ax1.secondary_xaxis('top', functions=(Helpers.wav2lam, Helpers.lam2wave))
        ax2.set_xlabel('wavelength [nm]')
        fig.tight_layout()
            
        if ylim:
            ax1.ylim(ylim)
        


        
class Helpers():
    @staticmethod
    def wav2lam(wn):
        # cm^{-1} to nm
        return 1.0e7 / wn
    @staticmethod
    def lam2wave(lam):
        # nm to cm^{-1} 
        return 1.0e7 / lam
     
    @staticmethod
    def hitran_molecule_number(name):
        """given molecule name, it returns the hitran id"""
        hitran_molecule_dic = {
         'H2O' : 1,  # Water
         'CO2' : 2,	# Carbon Dioxide	
         'O3' : 3,	# Ozone	
         'N2O' : 4,	# Nitrous Oxide	
         'CO' : 5,	# Carbon Monoxide	
         'CH4' : 6,	# Methane	
         'O2' : 7,	# Oxygen
         }
        """         
         8	NO	Nitric Oxide	
         9	SO2	Sulfur Dioxide	
         10	NO2	Nitrogen Dioxide	
         11	NH3	Ammonia	
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