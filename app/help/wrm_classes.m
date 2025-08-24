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
% * *WRM_SedimentTransport*: compute alongshore and cross-shore transport for an array of points on the edge of the surf zone. Uses the output of a Batch Wave Run as the nearshore wave input.
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
% * _wrm_transport_plots_ - uses the sediment transport results for a set of points along the coast to 
% examine drift rates, the divergence of drift and the Peclet number (indicates balance of advection and diffusion).

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
% * _extract_wave_data_ - extract Hs, Tp and Dir from a dataset that does not use default naming
%   convention (e.g. Copernicus re-analysis data).

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
%% Core grid functions
% * *gd_ax_dir*
% - check direction of grid axes and reverse if descending, OR
% find grid orientation using ishead and direction of x-axis, OR
% check a grid axis direction by prompting user.
% * *gd_centreline.m*
% - create a centreline of a channel using function _a_star_ to trace the
% shortest path between start and end points whilst finding the deepest
% points (i.e. a thalweg).
% * *gd_colormap*
% - check if Mapping toolbox is installed to use land/sea colormap, or call
% _cmap_selection_ if not available (see <matlab:doc('psfunctions') Plotting and statistical functions> 
% in the <matlab:doc('muitoolbox') muitoolbox>).
% * *gd_convergencelength* 
% - least squares fit using fminsearch to % find the convergence length of 
% a channel from a distance-width xy data set.
% * *gd_digitisepoints*
% - creates figure to interactively digitise points on a grid and add
% elevations if required.
% * *gd_dimensions*
% - get the grid dimsnions for a grid struct (as used in GDinterface).
% * *gd_getpoint.m*
% - interactively select a point on a plot and return the point
% coordinates.
% * *gd_grid_line*
% - create a grid from scattered data points input as xyz tuples.
% * *gd_plotgrid*
% - create pcolor plot of gridded surface.
% * *gd_plotsections*
% - display grid and allow user to interactively define start and
% end points of a section line to be plotted in a figure.
% * *gd_pnt2vec.m*
% - convert an array of structs with x,y (and z) fields to a [Nx2] or [Nx3] 
% array of points, or a single stuct with vectors for the x, y (and z)
% fields.
% * *gd_read_image.m*
% - read cdata for an image from an ASCII text file, with the positional
% information.
% * *gd_readshapefile.m*
% - read the x and y coordinates from a shape file. Lines are concatenated
% and separated by NaNs in single x and y vectors. Suitable for reading
% boundaries or sections into a single array.
% * *gd_selectpoints*
% - accept figure to interactively create a specified number of x,y points
% on a grid.
% * *gd_setpoint*
% - interactively select a single point on a plot and return the point
% coordinates. Includes an option to enter an additional value at the
% selected point (e.g. for elevation).
% * *gd_setpoints.m*
% - interactively create a set of points on a plot and return the point
% coordinates. Includes an option to enter an additional value at the
% selected points (e.g. elevation).
% * *gd_startendpoints*
% - accept figure to interactively select start and end points on a grid.
% * *gd_subdomain*
% - accept figure to interactively select a subdomain of a grid.
% * *gd_subgrid*
% - extract a subdomain from a grid and return the extracted
% grid and the source grid indices of the bounding rectangle.
% * *gd_write_image.m*
% - write image to a jpg or tif file.
% * *gd_write_image_ascii.m*
% - write cdata from an image to an ASCII text file, with the positional 
% information.
% * *gd_write_tiff.m*
% - write the data from an image to a tif file including the position
% information.
% * *gd_xy2sn*
% - map grid from cartesian to curvilinear coordinates with option to return 
% the elevations on the source cartesian grid, or as a curvilinear grid.
% * *gd_sn2xy*
% - map grid from curvilinear to cartesian coordinates.

%% Additional utility functions
% * *gd_lineongrid_plot*
% - plots a defined line onto a countour or surface plot of a grid (e.g a
%   channel centre-line).
% * *gd_user_function*
% - function for user to define bespoke use of grids and grid tools.


%% Functions from Matlab(TM) Exchange Forum
%%
% * *a_star*
% - implements the A* search algorithm to find the shortest path given
% constraints (inaccessible cells) and a cost function (e.g. water depths).
% Author: Alex Ranaldi, 2022, https://github.com/alexranaldi/A_STAR
% * *InterX* 
% - intersection of two curves. MATLAB Central File Exchange, 
% Author: NS, 2010, https://www.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections.
% * *xy2sn* 
% - Bart Vermeulen,2022, Cartesian to Curvilinear 
%   coordinate forward and backward transformation. 
%   https://www.mathworks.com/matlabcentral/fileexchange/55039-cartesian-to-curvilinear-coordinate-forward-and-backward-transformation 
% * *sn2xy* 
% - as above.

%%
% Further details can be found in <matlab:doc('grid_class_fcns') Grid classes and functions>
% 

%% See Also 
% The ModelSkill <matlab:wrm_open_manual manual>, which provides further details 
% of setup and configuration of the model.