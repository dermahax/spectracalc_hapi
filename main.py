# -*- coding: utf-8 -*-
"""
#mber 27.10.2022
OOP spectracalc representation using hapi. see spectracalc.com

# worklflow:
    # 1. set up spectra class
    # 2. add observer
    # 3. add gass_cells & gasses
    # 4. download data
    # 5. plot
"""


#db_begin ('data') # sets the folder for the hapi tables
from classes import Spectra, Observer, HiddenPrints


#%% 1. create a spectrum
my_spectrum = Spectra(name = 'Methane_Line_1_')

#%% 2. add an observer
my_spectrum.observer = Observer(nu_min =  3055.0,    # [1/cm]lower wavenumber
                                nu_max =  3069.0,    # [1/cm]upper wavenumber
                                )

#%% 3. add gas_cells (may consist of multiple gasses) with gasses
# cell 0
my_spectrum.add_gas_cell(temperature    = 296,          # K
                         pressure       = 1,            # atm
                         length         = 100,          # cm
                         no_gasses      = 1)            # number of gasses
# add gass(es) to the cell (to the last cell added)
my_spectrum.gas_cells[-1].add_gas(gas_name  = "H2O",    # as in hitran
                                  VMR       = 1E-2)
# more gasses could be added in one cell.

# cell 1
my_spectrum.add_gas_cell()
my_spectrum.gas_cells[-1].add_gas(gas_name  = "CH4",
                                  VMR       = 100E-6)

# cell 2...

#%% 4. download (running hapi functions)
with HiddenPrints():            # disables hapi stdout
    my_spectrum.download(line_list = True) # when disabling line list, it will not be plotted

#%% 5. plot
my_spectrum.plot()

#%% 6. summary
print(my_spectrum)
