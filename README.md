# impedanceAdjustment
Code that corresponds to "An impedance matching algorithm for common-mode interference removal in vagus nerve recordings"

impedance_adjustment.m  - impedance adjustment algorithm

reproduce_figures.m     - code that reproduces the figures in the paper given the complete recordings
                          this code will not run completely because the full recordings are not provided
                          but a user can see how impedance_adjustment is called with the sample data

ecgArtifacts.mat         - sample ECG artifact data
stimulationArtifacts.mat - sample stimulation artifact data

MessageUpdater.m        - a class that displays the status used in impedance_adjustment
assoiate_points.m       - a dependency that associates ECG artifacts peaks between channels

This third party function is also required to use this code
https://www.mathworks.com/matlabcentral/fileexchange/45577-inverse-short-time-fourier-transform-istft-with-matlab
The Matlab R2019a istft function has not been tested with this impedance adjustment algorithm