# Hyperspectral Data Processing

## Matlab codes
The group of matlab codes does nothing functionally differently than the current hyperspectral processing codes in the Link / Landes Group. I have recreated the matlab codes in python. I verified that the code exactly reproduces matlab. To check yourself, follow along with the Jupyter Notebook, `Match Matlab code and Python`.


## How to run Matlab
* First run `Hyper_Analysis.m`. This file does that data processing that turns the raw data into the correctly normalized data.
	* To do this, update line 8 `fold` to the folder that contains your raw data.
	* Create a file in that folder called `mydata.txt` that contains the white counts, dark counts, and data file.
	* A few files now should be greated in this folder.

* Next, run `plot_particle_spectra_basic.m`.
	* Update lines 8-9 to be the folder that contains the raw data, and another folder that you'll make of the same name as the data file.
	* Make the above folder inside the data folder.
	* Make a file called positions.txt, you'll populate it shortly.
	* Run the first two blocks in matlab. One by one, identify the coordinates of each NP and place them in the positions text file.
	* Once this is done, run the rest of the script.


## How to run Python
