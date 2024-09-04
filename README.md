# WaveRayModel
The WaveRayModel is a coastal wave model that uses forward and backward tracking of wave rays (orthogonals to wave crest) over a bathymetric depth grid. In backtracking mode offshore wave spectra can be transformed to inshore points.

## Licence
The code is provided as Open Source code (issued under a BSD 3-clause License).

## Requirements
WaveRayModel is written in Matlab(TM) and requires v2016b, or later. In addition, WaveRayModel requires the _dstoolbox_ and the _muitoolbox_. Meshes make use of mesh2d by Darren Engwirda and this is included as a submodule in the repository, linked to https://github.com/dengwirda/mesh2d.

## Background
WaveRayModel is a MatlabTM App to examine wave refraction and shoaling over coastal bathymetries. Using a ray tracing method forward and backward tracking rays can be generated for a range of wave periods (frequencies), directions and water levels. Forward tracking provides a spatial visualisation of how waves approach the shore. Backward tracking is used to construct transfer tables that allow inshore wave conditions to be estimated from an offshore timeseries of waves (height, period and direction). Several different definitions of the offshore wave spectrum are provided and depth saturation can also be included (TMA spectrum). The App includes options to plot the offshore and inshore spectrums and the various transfer coefficients as a function of wave period, mean directions and water level. Note, however, that this implementation does not include processes of wave reflection, wave diffraction and wave breaking, although the latter is implicitly captured by including depth saturation using the TMA spectrum.

<img width="454" alt="waveraymodel_images" src="https://github.com/user-attachments/assets/13f41271-44d4-489e-b8ba-4ffbc6780bb9">

## Manual
The WaveRayModel manual in the app/doc folder provides further details of setup and configuration of the model. The files for the example use case can be found in the app/example folder. 

## See Also
The repositories for _dstoolbox_, _muitoolbox_, _muiAppLIb_ and https://github.com/dengwirda/mesh2d.