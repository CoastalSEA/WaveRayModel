%% WaveRayModel Classes and Functions
% The WaveRayModel App is built using the <matlab:doc('muitoolbox') muitoolbox>, 
% the <matlab:doc('dstoolbox') dstoolbox> and a number of model specific classes.

%% WaveRayModel classes
% * *WaveRayModel*: defines the behaviour of the main UI.
% * *MS_RunParams*: UI to set pareameters used for the skill score
% * *Ray*: creates individual ray objects. Called by RayTracks and uses arc_ray function.
% * *RayTracks*: constructing array of wave ray tracks as a function of wave direction, wave period and water level, working in forward or backward ray tracking mode. 
% * *SpectralTransfer*: builds the offshore and inshore Transfer Tables from a backward ray tracking data set (class RayTracks) for use in WRM_WaveModel. Also has a method to create plots of the transfer coefficients for a unit wave height.
% * *WRM_Bathy*: generate idealised bathymetries using a linear slope or Dean's profile on a linear or crenulate shoreline. Also includes an option for a mound on a linear slope.
% * *WRM_BT_Params*: parameters define initial positions and direction of backward tracking rays.
% * *WRM_FT_Params*: parameters define initial positions and direction of forward tracking rays.
% * *WRM_Mesh*: generate triangular mesh of nearshore area from a Cartesian grid.
% * *WRM_RunParams*: run parameters for the WaveRayModel.
% * *WRM_WaveModel*: Class for wave refraction using backward ray transfer function. Constructs inshore time series from an offshore timeseries. Also includes methods to plot the offshore and inshore spectra.

%% WaveRayModel functions
% Functions to construct rays, wave sprectra and generate plots
%%
% * _arc_ray_ – compute the exit position and direction of a ray entering a triangular element at the position and direction defined by the incoming ray.
% * _celerity_grid_ – calculate the celerity and group celerity over a bathymetry grid for a range of wave periods and a range of water levels.
% * _celerity_mesh_ - calculate the celerity and group celerity over a bathymetry mesh for a range of wave periods and a range of water levels
% * _compass2trig_ – convert compass directions to trigonometric angles, or vice versa.
% * _curvspace_ -  evenly spaced points along an existing curve in 2D or 3D., Author: Yo Fukushima, 2005, https://www.mathworks.com/matlabcentral/fileexchange/7233-curvspace 
% * _datawell_directional_spectra_ - estimates the directional distribution of a wave spectrum for directions, dirs, given the mean, spread, skewness and kurtosis parameters as output by datawell buoys SPT file format.
% * _directional_spreading_ – sample a directional spreading function at selected direction intervals.
% * _get_element_ – define triangular polyshape based on quadrant being entered by ray.
% * _get_inshore_spectrum_ - construct the offshore and inshore spectra for given wave conditions or wave buoy spectral data.
% * _get_inshore_wave_ - integrate the 2-D spectra to obtain wave parameters and transfer coefficients. 
% * _get_intersection_ – find intersection of a triangle element and a line segment that can be a straight line or an arc segment.
% * _get_quadrant_ – find the quadrant that the start point lies in or on. Once first intersection has been found subsequent quadrants are defined in next_element, which calls get_quadrant if ray direction is aligned to axis.
% * _interpolate_cbrewer_ – interpolate a colorbrewer map to ncolors levels, Charles Robert, 2011, part of  cbrewer from Matlab(TM) Exchange Forum.
% * _is_axis_point_ – find whether point lies on an axis and is travelling in the direction of that axis.
% * _isangletol_ – boolean check of whether an angle lies between upper and lower bounds defined as specific angles, or a tolerance.
% * _isclosetol_ – boolean check of whether an x,y point lies close to another x,y point
% * _mesh_arc_ray_ – compute the exit position and direction of a ray entering a trian.gular mesh element at the position and direction defined by the incoming.
% * _mesh2d_ – generate triangular mesh, Author: Darren Engwirda, 2017, https://github.com/dengwirda/mesh2d. NB: the mesh2d folder in the App is empty and this set of functions needs to be downloaded from the source and installed in the designated sub-folder of the WaveRayModel App.
% * _PoolWaitbar_ - initiallises and updates a waitbar when running a loop using parfor, Matlab(TM) Forum: Edric Ellis, 2019
% * _trigradient_ – compute the gradients in the x and y directions using triangular mesh as input, Author: Mick Warehime, 2013, https://www.mathworks.com/matlabcentral/fileexchange/36837-trigradient-m .
% * _wave_spectrum_ – calculate the spectral energy at a number of frequencies using a selection of spectrum definitions (Bretschneider open ocean, Pierson-Moskowitz fully developed, JONSWAP fetch limited, and TMA shallow water).
% * _wave_spectrum_params_ - integrate a 2-D spectra to obtain wave parameters
% * _wrm_animation_ - animation of model spectra timeseries.
% * _wrm_runmovie_ - callback function for animation figure buttons and slider modified from muiPlots to handle two subplots.


%% Additional Coastal Classes
% In the folder ../muiAppCoastalClasses.
%%
% * *ctWaterLevelData* – import water level data
% * *ctWaveData* – import wave data
% * *ctWindData* – import wind data

%% Additional Coastal Functions
% In the folder ../muiAppCoastalFcns.
%%
% * _celerity_ - calculate the wave celerity using Hunt's equation.
% * _crenulate_bay_ - generate the shoreline for an equilibrium crenulate bay using the method of Hsu and Evans, 1989 (Hsu and Evans, 1989).
% * _data_cleanup_ – set of functions to help clean up timeseries
% * _deanbeachprofile_ - find the bed slope across the surf zone.
%   the profile is based on a user defined slope between Hw and SWL (0mOD)
%   and a Dean profile below this level. 

%% Grid Classes
% Classes used to manipulate cartesian grids can be found in the
% _muiAppGridFcns_ folder and include the following:
%%
% * *GD_GridProps*: class inherits <matlab:doc('muipropertyui') muiPropertyUI> 
% abstract class, providing an interface to define the extent and intervals
% of a cartesian grid. 
% * *FGDinterface*: an abstract class to support classes that need additional
% functionality to handle grids. The class inherits *GDinterface*, which
% in turn inherits <matlab:doc('muidataset') muiDataSet> 
% to provide an extensive set of methods to handle datasets
% of various types (eg from models or imported files). 
% * *GD_ImportData*: class inherits <matlab:doc('fgdinterface') FGDinterface> abstract class (see above)
% to load xyz data from a file.

%% Grid Functions
% Functions used to manipulate cartesian grids can be found in the
% _muiAppGridFcns_ folder and include the following:
%%
% * _gd_ax_dir_
% - check direction of grid axes and reverse if descending, OR
% find grid orientation using ishead and direction of x-axis, OR
% check a grid axis direction by prompting user.
% * _gd_basin_hypsometry_
% - compute area and volume hypsometry from gridded elevation data.
% * _gd_basin_indices_
% - get the indices of the grid x-axis that fall within the basin or channel,
% when the mouth is offset from the grid origin. (NB: assumes basin/channel
% is aligned iwth the x-axis). Also returns the index of mouth position on 
% the x-axis.
% * _gd_basin_properties_
% - use the basin hypsometry from gd_basin_hypsometry to compute several 
% along-channel/x-axis morphological properties.
% * _gd_colormap_
% - check if Mapping toolbox is installed to use land/sea colormap, or call
% _cmap_selection_ (see <matlab:doc('psfunctions') Plotting and statistical functions> 
% in the <matlab:doc('muitoolbox') muitoolbox>) if not available.
% * _gd_digitisepoints_
% - creates figure to interactively digitise points on a grid and add
% elevations if required.
% * _gd_dimensions_
% - get the grid dimsnions for a grid struct (as used in GDinterface).
% * _gd_grid_line_
% - create a grid from scattered data points input as xyz tuples.
% * _gd_gross_properties_
% - compute the gross properties of a gridded bathymetry.
% * _gd_plan_form_ 
% - compute planform variation along the x-axis at specified planar levels.
% * _gd_plotgrid_
% - create pcolor plot of gridded surface.
% * _gd_plotsections_
% - display grid and allow user to interactively define start and
% end points of a section line to be plotted in a figure.
% * _gd_property_plots_
% - plots displayed on Proprety tab or stand-alone figure in Apps that use 
% GDinterface, such as ChannelForm and ModelSkill.
% * _gd_section_properties_
% - compute the width, cross-sectional area and prism along channel.
% * _gd_selectpoints_
% - accept figure to interactively select one or more points on a grid.
% * _gd_setpoint_
% - interactively select a point on a plot and return the point
% coordinates. Includes an option to enter an additional value at the
% selected point (e.g. for elevation).
% * _gd_startendpoints_
% - accept figure to interactively select start and end points on a grid.
% * _gd_subdomain_
% - accept figure to interactively select a subdomain of a grid.
% * _gd_property_plots_ 
% - plots displayed on Proprety tab in ChannelForm model and on a figure 
% in ModelSkill.
% * _gd_xy2sn_
% - map grid from cartesian to curvilinear coordinates with option to return 
% the elevations on the source cartesian grid, or as a curvilinear grid.
% * _gd_sn2xy_
% - map grid from curvilinear to cartesian coordinates.
% * _getconvergencelength_ 
% - least squares fit using fminsearch to % find the convergence length of 
% a channel from a distance-width xy data set.
% * _getsubgrid_
% - extract a subdomain from a grid and return the extracted grid and the 
% source grid indices of the bounding rectangle.
% * _a_star_
% - implements the A* search algorithm to find the shortest path given
% constraints (inaccessible cells) and a cost function (e.g. water depths).
% Author: Alex Ranaldi, 2022, https://github.com/alexranaldi/A_STAR
% * _InterX_ 
% - intersection of two curves. MATLAB Central File Exchange, 
% Author: NS, 2010, https://www.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections.
% * _xy2sn_ 
% - Bart Vermeulen,2022, Cartesian to Curvilinear 
%   coordinate forward and backward transformation. 
%   https://www.mathworks.com/matlabcentral/fileexchange/55039-cartesian-to-curvilinear-coordinate-forward-and-backward-transformation 
% * _sn2xy_ 
% - as above.

%%
% Further details can be found in <matlab:doc('grid_class_fcns') Grid classes and functions>
% 

%% See Also 
% The ModelSkill <matlab:wrm_open_manual manual>, which provides further details 
% of setup and configuration of the model.