# SPSPy
SPSPy is a Python based package of tools for use with the Super-Enge Split-Pole Spectrograph at FSU. Much of the code here is based on Java programs originally written at Yale University by D.W. Visser, C.M. Deibel, and others. Currently the package contains SPSPlot, a tool aimed at informing users which states should appear at the focal plane of the SESPS, and SPANC, a tool for calibrating the position spectra from the focal plane.

Updated 10/2025 to fix incompatibilities of outdated libraries with Python v3.12.

## Installation
### Clone Repository 
	git clone https://github.com/jce10/spspy.git
	cd spspy
### Create and Load Python Virtual Environment and Install Required Python Packages
	python3 -m venv .venv
	source .venv.bin/activate
	pip install -r requirements.txt
Note: `pip` is highly recommended as the package installer. 
Note x2: You can name your virtual enviornment something other than ".venv" -- edit this to whatever you'd like. To deactivate the environment, use the command `deactivate` in the active terminal. 

## Nuclear Data
Based on the inputs from the user, `NuclearData.py` utilizes the Livechart Data Download [API](https://www-nds.iaea.org/relnsd/vcharthtml/api_v0_guide.html) services from IAEA. Currently, the API allows for CSV formatting, with development for other format types under construction. `NuclearData.py` fetches levels for the inputted reaction, the CSV is generated, and levels saved to a `pandas` `DataFrame` for future use. 

## SPSPlot
This tool is intended to be used for guiding the settings of the SPS to show specific states on the focal plane detector. The user inputs reaction information, and the program runs through the kinematics to calculate the energies of ejecta into the the SESPS via the [PyCatima](https://github.com/hrosiak/pycatima) module. Updated 10/2025 with pyCatima v1.98. 

SPSPlot displays excitation energies of the state, the kinetic energy of the ejectile, or the focal plane z-offset for a state. SPSPlot can also export the calculated reaction information to a csv file. 
Note: that since levels are obtained from NNDC, SPSPlot requires an internet connection.

## SPANC
SPANC is the program used to calibrate SESPS focal plane spectra. It works by the user specifying a target, reaction, calibration peaks, and output peaks. The target is a description of the physical target foil used in the SPS, which is used to calculate energy loss effects. 

The target must contain the isotope used as the target in the reaction description. The reaction indicates to the program what type of ejecta are expected, as well as the settings of the spectrograph. Calibration data is given as centroids from a spectrum with correspoding excitation energies, as well as associated uncertainties. The calibration peaks are then fit using the scipy ODR package (see scipy ODR for more documentation). The fit is plotted, and the results are shown in a table. Additionally, residuals are plotted and shown in a table. The user can then feed the program an output peak, or a peak for which the user would like to calculate the excitation energy of a state using the calibration fit. The peak excitation energy will then be reported, with uncertainty. The user can also give a FWHM to be converted from focal plane position to energy. 


## Running the tools
Activate the environment containing the requirements and then simply `python main.py` from the repository. This will bring up a launcher from which you can select a tool.

### Known issues
	1. Saving SPSPlot and SPANC projects for future use is unavailable right now. This still being fixed. 
