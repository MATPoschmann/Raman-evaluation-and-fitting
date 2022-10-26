# Raman-evaluation-and-fitting

The program opens every .csv (wavelength, intensity, no header) file in its own folder and creates meanspectra of these.
In between it shows a graphic with each spectrum, so you can compare and exclude single ones from the evaluation. 
The program results in meanspectra with standard deviation as graphic- and text-file. 

The program is specialized to evaluate data from carbon materials.
It finds every maxima in the spectrum using the numpy function "argrelexxtrema".
It takes intensity values of D and G peak and calculates I_G over I_D ratio.
It fits a Breit-Wigner-Fano and a Lorenz peak to the G and D peak and returns the fitting values.
Results are saved as text and graphic file.
