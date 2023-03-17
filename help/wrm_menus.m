%% Menu Options
% Summary of the options available for each drop down menu.

%% File
% * *New*: clears any existing model (prompting to save if not already saved) and a popup dialog box prompts for Project name and Date (default is current date). 
% * *Open*: existing Asmita models are saved as *.mat files. User selects a model from dialog box.
% * *Save*: save a file that has already been saved.
% * *Save as*: save a file with a new or different name.
% * *Exit*: exit the program. The close window button has the same effect.

%% Tools
% * *Refresh*: updates Cases tab.
% * *Clear all > Project*: deletes the current project, including all Setup data and all Cases.
% * *Clear all > Figures*: deletes all results plot figures (useful if a large number of plots have been produced).
% * *Clear all > Cases*: deletes all Cases listed on the Cases tab but does not affect the model setup.

%% Project
% * *Project Info*: edit the Project name and Date
% * *Cases > Edit Description*: user selects a Case to edit the Case description.
% * *Cases > Edit Data Set*: initialises the Edit Data UI for editing data sets.
% * *Cases > Save*: user selects a data set to be saved from a list box of Cases and the is then prompted to name the file. The data are written to an Excel spreadsheet. 
% * *Cases > Delete*: user selects Case(s) to be deleted from a list box of Cases and results are then deleted (model setup is not changed).
% * *Cases > Reload*: user selects a Case to reload as the current parameter settings.
% * *Cases > View settings*: user selects a Case to display a table listing the parameters used for the selected Case. 
% * *Export/Import > Export*: user selects a Case class instance to export as a mat file.
% * *Export/Import > Import*: user selects an exported Case class instance (mat file) to be loaded.
%%
% *NB*: to export the data from a Case for use in another application 
% (eg text file, Excel, etc), use the *Project>Cases>Edit Data Set* option 
% to make a selection and then use the ‘Copy to Clipboard’ button to paste 
% the selection to the clipboard.


%% Setup > Import Data > Bathymetry
% Options to Load, Add etc are similar to timeseries data. This option
% loads a grid of x, y, z data.

%% Setup > Import Data for timeseries
% Load data for Waves, Water Levels, and Winds.To create a new instance (e.g. for a
% different location or data source) use Load. To add data to an
% existing data set, use Add.
%%
% * *Load data*: prompts user for file format to be loaded. The options available vary with Data type. The user then selects one or more files to be loaded. Once files have been read, user is prompted for a description (working title) for the data set. 
% * *Add data*: prompts user for file to be loaded (only one file at a time can be added). Only files with the format used to create the data set can be used to Add data to a data record.
% * *Delete data*: prompts user for Case to be deleted.
% * *Quality control*: runs a series of checks on the data. This is only available for waves and water levels in this version. See <matlab:ct_open_manual manual> for details of the checks applied. Once run, the last column of the Data tab is updated to show that qc has been applied.lly, the default values are appropriate but these can be adjusted and saved with the project if required.
%% Setup > Grid Tools
% * *Grid Tools > Translate Grid*: interactively translate grid x-y
% coordinates;
% * *Grid Tools > Rotate Grid*: interactively flip or rotate grid;   
% * *Grid Tools > Re-Grid*: regrid a gridded dataset to match another grid or to user
% specified dimensions;
% * *Grid Tools > Sub-Grid*: interactively define a subgrid and save grid as a new Case;               
% * *Grid Tools > Combine Grids*: superimpose one grid on another based on maximum
% or minimum set of values;
% * *Grid Tools > Add Surface*: add horizontal surface to an extisting
% grid;
% * *Grid Tools > To curvilinear*: map grid from cartesian to curvilinear coordinates; 
% * *Grid Tools > From curvilinear*: map grid from curvilinear to cartesian
% coordinates;
% * *Grid Tools > Display Dimensions*: display a table with the dimensions
% of a selected grid;
% * *Grid Tools > Difference Plot*: generate a plot of the difference
% between two grids;
% * *Grid Tools > Plot Sections*: interactively define sections and plot
% them on a figure;
% * *Grid Tools > Digitise Line*: interactively digitise a line (with
% option to add z values) using selected grid as base map;
% * *Grid Tools > Export xyz Grid*: select a Case and export grid as xyz
% tuples;

%% Setup > Run Parameters
% * *Run Parameters > Run Conditions*: used to define the range of wave periods and water levels, together with the number of intervals to use. The ray cutoff water depth and the shoreline angle are also defined in this menu.
% * *Run Parameters > Forward Tracking*: used to define the coordinates of a transect along which the rays start, together with the wave direction and the number of rays.
% * *Run Parameters > Backward Tracking*: used to define the range of directions and number of rays, along with the coordinates of the start point.

%% Setup > Data Cleanup
% * *Concatenate two timeseries*: allows two timeseries data sets to be joined. Two records, or two variables, of a similar type can be joined to form a single timeseries or timeseries collection. 
% * *Resample timeseries*: allows a selected timeseries to be resampled at user specified interval, using a user specified method (e.g. mean/max/min over the interval). 
% * *Patch timeseries*: allows gaps in a selected timeseries to be patched using the data from another timeseries that overlaps the primary timeseries (at least for some or all of the gaps). 
% * *Trim timeseries*: clip the ends of a time series between specified dates.

%% Setup (other)
% * *Grid Parameters*: dialogue to set dimensions of default grid
% * *Model Constants*: a number of constants are used in the model. Generally, the default values are appropriate but these can be adjusted and saved with the project if required.


%% Run
% * *Check Start Points*: utility generates a plot showing the bathymetry, start points and initial ray directions.
% * *Forward Rays*: uses the defined run parameters (Setup>Run Parameters>Run Conditions) and start points (Setup>Run Parameters>Forward Tracking) to generate a set of rays. These can then be viewed on the Q-Plot tab, or using the Analysis>Ray Plots option.
% * *Backward Rays*: uses the defined run parameters (Setup>Run Parameters>Run Conditions) and start point (Setup>Run Parameters>Backward Tracking) to generate a set of rays. These can then be viewed on the Q-Plot tab, or using the Analysis>Ray Plots option.
% * *Transfer Table*: generates the Transfer Table from a set of Backward Tracking rays.  Summary output from the tables can then viewed on the Q-Plot tab, or using the Analysis>Spectral Plots>Transfer Table option.
% * *Run Property Timeseries*: uses a selected Transfer Table and imported wind or wave data to compute the inshore conditions which are saved as a timeseries 
% * *Run Spectra Timeseries*: uses a selected Transfer Table and imported
% directional wave spectra data to compute the inshore spectra, which are
% saved as spectra and property data (with the same definitions as the
% input timeseries).
% * *Derive Output*: initialises the Derive Output UI to select and define manipulations of the data or call external functions and load the result as new data set.

%% Analysis
% * *Plots*: initialises the Plot UI to select variables and produce various types of plot. The user selects the Case, Dataset and Variable to used, along with the Plot type and any Scaling to be applied from a series of drop down lists, 
% * *Statistics*: initialiss the Statistics UI to select data and run a range of standard statistical methods.
% * *Ray Plots*: used to select a case and plot the rays for a selected wave period and water level, or all rays for all periods and levels. There is also the option to plot the celerity or group celerity variation along a single set of rays.
% * *Spectral Plots > Transfer Table*: uses a Transfer Table created using _Run>Spectral Transfer>Create Table_ to plot data from the inshore or offshore Transfer Tables. For a selected water level the options are for celerity or group celerity inshore and direction, celerity, group celerity, offshore depth (seaward end of ray), average depth along rays and minimum depth along rays.
% * *Spectral Plots > Transfer Coefficients*: uses a selected Transfer Table and a unit wave height (Hs=1m) to derive the inshore wave parameters and return a set of transfer coefficients. These can then be plotted as a function of wave period and mean offshore direction for a selected water level.
% * *Spectral Plots > O/I Spectrum*: uses a selected Transfer Table and user defined input conditions for wind or wave forcing to derive a 2D offshore spectrum and the associated inshore spectrum, which are then plotted along with details of the conditions.
% * *Spectral Plots>O/I Animation*: uses a selected Transfer Table and an imported time series of wind or wave conditions to generate an animation of offshore and inshore spectrum (frames are similar to the above image). The user is prompted to select the type of spectrum to use and the range of data to use. If the record selected has more than 5000 data points the user is asked to confirm. This is because the plotting is data intensive and can quickly create large arrays. However, it does provide a useful way of checking what is happening in the transfer process.

%% Help
% * *Help*: access the online documentation for CoastalTools.

%% See Also
% The <matlab:wrm_open_manual manual> provides further details of setup and 
% configuration of the model.