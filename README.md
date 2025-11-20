Creates synthetic data from Exponentially modified Gaussian functions found in "Improved peak shape fitting in alpha spectra" by Pommé et al. 
Values for data taken from "Novel device to study double-alpha decay at the FRS Ion Catcher" by Varga et al. 
Synthetic data is placed into histogram which is then fitted using ROOT in C++. https://root.cern.ch/

Improvements need to be made to fit. Particularly in the first and last peaks. Suspect you need to make improvements to the initial guess parameters for the fit.
