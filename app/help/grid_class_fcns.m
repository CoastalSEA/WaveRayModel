%% Grid classes and functions
% The _muiAppGridFcns_ folder contains classes and functions to support the 
% handling of gridded data, for Apps built using the 
% <matlab:doc('muitoolbox') muitoolbox>.

%% Grid classes
% The Apps that use grids can make use of the following classes:
%%
% *GD_GridProps*: class inherits <matlab:doc('muipropertyui') muiPropertyUI> 
% abstract class, providing an interface to define the extent and intervals
% of a cartesian grid. Methods within the class include:
%%
% * setInput - calls the default interface to input property values for grid
% extend and grid intervals in the x and y directions.
% * setGridDimensions - pass the x and y dimensions of a grid to updatge the
% property settings.
% * getGridDimensions - retrieve the x and y dimensions as vectors and the
% grid intervals (delx, dely).

%%
% *GDinterface*: an abstract class to support classes that need additional
% functionality to handle grids. The class inherits <matlab:doc('muidataset') muiDataSet> 
% and together they provide an extensive set of methods to handle datasets
% of various types (eg from models or imported files). The methods in the
% <matlab:doc('gdinterface') GDinterface> include the following:
%%
% Methods to add, delete and update grids and propetry tables:
%%
% * setGrid - load the grid results into a dstable. 
% * getGrid - retrieve a grid.  
% * addGrid - add a grid to an existing <matlab:doc('muitoolbox') muitoolbox> Case table.
% * deleteGrid - delete a grid to an existing <matlab:doc('muitoolbox') muitoolbox> Case table.

%%
% Methods to manipulate grids:
%%
% * translateGrid - interactively translate grid x-y coordinates.
% * rotateGrid - interactively flip or rotate grid.
% * reGridData - regrid a gridded dataset to match another grid or to user
% specified dimensions.
% * subGridData - interactively define a subgrid and save grid as a new Case.
% * addSurface - add horizontal surface to an extisting grid. 
% * curvilinear_xy2sn - map grid from cartesian to curvilinear coordinates.
% * curvilinear_sn2xy - map grid from curvilinear to cartesian coordinates.
% * exportGrid - select a Case and export grid as xyz tuples.

%%
% Static methods to manipulate grids
%%
% * setCombinedGrids - superimpose one grid on another based on maximum
% or minimum set of values.
% * diffGridsPlot - generate a plot of the difference between two grids.
% * getGridLine - interactively digitise a line and save to a file.
% * gridMenuOptions - handles the menu calls ffor the default menu of Grid
% Tools. 

%%
% *GD_ImportData* class inherits GDinterface abstract class (see above)
% to load xyz data from a file. Methods within the class include:
%%
% * loadData - set up a new dataset by loading data from a file. 
% * addData - add additional data to an existing dataset.
% * tabPlots - plots the grid as a filled contour plot (can be overloaded if
% required for a specific application).

%%
% These methods are supplemented by the methods available in
% <matlab:doc('muidataset') muiDataSet> as well as the
% <matlab:doc('gdinterface') GDinterface> summarised above.

%%
% *FGDinterface*: an abstract class to support classes that need additional
% functionality to handle grids. The class inherits *GDinterface*, to extend 
% the grid handling methods, to extract inlet/channel properties
% from the grid. The setGrid, getGrid, addGrid and deleteGrid methods in
% *GDinterface* are overloaded to handle inlet/channel specific properties
% such as distance to mouth of inlet, orientation of channel, etc.
% Other methods in the <matlab:doc('fgdinterface') FGDinterface> include the following:
%%
% * setProperties - initialises the tables for derived properties of gridded data set
% * addProperties - add properties to an existing set of grid property tables
% * delProperties - delete a row of properties from existing grid property
% tables
% * addFormProps - add form properties to a gridded data set
% * delFormProperties - delete ALL property tables associated with a 
% selected gridded data set
% * setModelFormProps -add a set of Hypsometry, Section and Gross form properties 
% * editGridInletData - UI to edit, definition of channel head, x-distance 
% to mouth and definition of centre-line (if used).

%%
% Static methods to manipulate grids
%%
% * addValleyBase - use an xyz definition of the channel thalweg to create 
% a valley base this is used to modify an existing valley defined above MHW 
% to include the valley form down to the pre-Holocene surface (i.e. below the
% existing channel bed).
% * addShore - add a strip to shore side of a grid by extrapolating from the
% neighbouring strip.
% * formMenuOptions - handles the menu calls ffor the default menu of Form
% Tools. 

%%
% *FGD_ImportData* class inherits FGDinterface abstract class (see above)
% to load xyz data from a file. Methods within the class include:
%%
% * loadData - set up a new dataset by loading data from a file. 
% * addData - add additional data to an existing dataset.
% * tabPlots - plots the grid as a filled contour plot (can be overloaded if
% required for a specific application).

%%
% These methods are supplemented by the methods available in
% <matlab:doc('muidataset') muiDataSet> as well as the <matlab:doc('fgdinterface') FGDinterface> 
% and <matlab:doc('gdinterface') GDinterface> summarised above.


%% Grid property functions
% Summary of functions that derive and manipulte properties of the grid 
% such as basin/channel dimensions and available in the _muiAppGridFcns_ 
% folder. Use the Matlab(TM) help function in the command window to get 
% further details of each function.
%%
% *gd_basin_hypsometry*
% - compute area and volume hypsometry from gridded elevation data and
% return as a cell array of along-channel values, containing _z-intervals_, 
% _zsurf(z)_ and _zvol(z)_, and a dstable of the surface area hyspometry, _SArea(x,z)_.
%%
%   [hyps,hdst] = gd_basin_hypsometry(grid,wl,histint,limits,isHW); %see below for explantion of input variables

%%
% *gd_basin_properties*
% - uses the hypsometry, _hdst_, obtained from gd_basin_hypsometry to compute a number
% of along-channel/x-axis morphological properties. Returns _cprops_, a table of
% cumulative (from head) along-channel/x-axis properties and _xprops_, a table of 
% incremental properties for each x-interval.
%%
%   [cprops,xprops] = gd_basin_properties(grid,wl,hdst); %see below for explantion of input variables

%%
% *gd_gross_properties*
% - compute the gross properties of a gridded bathymetry
%%
%   grossprops = gd_gross_properties(grid,wl,props); %see below for explantion of input variables

%%
% *gd_section_properties*
% - compute the width, cross-sectional area and prism along channel.
%%
%   props = gd_section_properties(grid,wl,hdst); %see below for explantion of input variables
%%
% where _grid_ is a struct containing the gridded data, _wl_ is a struct 
% for water levels at 3 elevations (high water, mean tide and low water), _histint_ is the 
% vertical interval to be used for the hypsometry, _limits_ defines the 
% upper and lower limits for vertical range of the hypsometry, _isHW_ is a flag to apply mask to data to 
% apply an upper bound to a given data set at high water within limits that 
% might be set for multiple grids. _props_ is the table returned by
% gd_section_properties and _hdst_ is the dstable returned by gd_basin_hypsometry.

%%
% *gd_property_plots*
% - plots displayed on Proprety tab or stand-alone figure in Apps that use 
% GDinterface, such as ChannelForm and ModelSkill.
%%
%   gd_property_plots(obj,irec,src)
%%
% where _obj_ is an instance of any form model that uses the GDinterface
% abstract class, _irec_ the row/timestep to use to access grid data and
% _src_ is a handle to the calling menu, create Figure button, or a figure
% handle.

%% Grid utility functions
% Summary of utility functions available in the _muiAppGridFcns_ folder. Use the Matlab(TM)
% help function in the command window to get further details of each
% function.

%%
% *gd_ax_dir*
% - check direction of grid axes and reverse if descending, OR
% find grid orientation using ishead and direction of x-axis, OR
% check a grid axis direction by prompting user.
%%
%   output = gd_ax_dir(ax,x,y,ischeck);
%%
% where (i) if _ax_ is a graphic axes handle the functions checks the direction
% of _x_ and _y_ and reverses _ax.XDir_ or _ax.YDir_ if descending, or (ii) 
% if _ax_ is a grid struct then uses grid.ishead and grid.x to determine
% orientation, or (iii) if _ax_ is empty the _x_ and _y_ values are assigned
% to as vecotrs, then if _ischeck_ is empty or omitted, user is 
% prompted to confirm each axis direction, whearas if _ischeck_ is a [1x2]
% logical array, for x and y axes respectively, reversing the axis 
% direction if true.
% To use the function for a single axis, assign the unwanted axis input to 
% be empty. See function Help for further details.

%%
% *gd_basin_indices*
% - get the indices of the grid x-axis that fall within the basin or channel,
% when the mouth is offset from the grid origin. (NB: assumes basin/channel
% is aligned with the x-axis). Also returns the index of mouth position, ixM, and 
% head, or tidal limit, position, ixH, on the x-axis.
%%
%   [idx,ixM,ixH] = gd_basin_indices(grid,Lt);    %Lt is optional
%% 
% where _grid_ is a struct of x, y, z (e.g. as used in getGrid in the
% GDinterface) and _Lt_ is the distance from the mouth to the head or tidal
% limit.

%%
% *gd_colormap*
% - check if Mapping toolbox is installed to use land/sea colormap, or call
% _cmap_selection_ (see <matlab:doc('psfunctions') Plotting and statistical functions> 
% in the <matlab:doc('muitoolbox') muitoolbox>) if not available.
%%
%   gd_colormap(zlimits);  %zlimits - minimum and maximum elevations [minz,maxz]

%%
% *gd_digitisepoints*
% - accept figure to interactively digitise points on a grid and add
% elevations if required.
%% 
%   points = gd_digitisepoints(grid,promptxt,isxyz,isdel);
%% 
% where _grid_ is a struct of x, y, z (e.g. as used in getGrid in the
% GDinterface); _promptxt_ is a cell array of prompts to be used for each
% point being defined; _isxyz_ is a logical flag true to input z values -
% optional, default is false; and _isdel_ is a logical flag true to delete
% figure on completion - optional, default is false.

%%
% *gd_dimensions*
% - get the grid dimensions for a grid struct (as used in classes that 
% inherit GDinterface). Outputs a table of grid dimensions including: 
% delx, dely - grid interval in the x and y dimensions;
% xsgn, ysgn - direction of x, y axes (-ve is descending);
% xint, yint - number of grid intervals in the x and y dimensions;
% xmin, xmax, ymin, ymax, zmin, zmax - minimum and maximum x, y and z
% dimensions.
%%
%   gdims = gd_dimensions(grid)

%%
% *gd_grid_line*
% - create a grid from scattered data points input as xyz tuples.
%%
%   [X,Y,Z] = gd_grid_line(xin,yin,zin,missing); %where xin, yin, zin are vectors
%   [X,Y,Z] = gd_grid_line(xyz,missing);         %where xyz is a an nx3 array     

%%
% *gd_plan_form* 
% - compute planform variation along the x-axis at specified planar levels.
%%
%   yz = gd_plan_form(grid,wl);   %see below for explantion of input variables

%%
% *gd_plotgrid*
% - create pcolor plot of gridded surface.
%%
%   ax = gd_plotgrid(hfig,grid);  %hfig is figure handle and ax is axes handle

%%
% *gd_plotsections*
% - display grid and allow user to interactively define start and
% end points of a section line to be plotted in a figure.
%%
%   gd_plotsections(grid);  %where grid is a struct of x, y, z

%%
% *gd_selectpoints*
% - accept figure to interactively select one or more points on a grid.
%%
%   points = gd_selectpoints(grid,npts,promptxt,isdel) %uses graphical selection

%%
% *gd_setpoint*
% - interactively select a point on a plot and return the point
% coordinates. Includes an option to enter an additional value at the
% selected point (e.g. elevation).
%%
%   point = gd_setpoint(ax,promptxt,isxyz);
%%
% where _ax_ is the figure axes to be uses,  _promptxt_ is the user prompt
% for the point being defined and _isxyz_ is a logical flag, which is true
% if a value is to be assigned to the point.
%% 
% *gd_startendpoints*
% - accept figure to interactively select start and end points on a grid.
%% 
%   points = gd_startendpoints(grid,isdel); %uses dialogue box for input

%%
% *gd_subdomain*
% - accept figure to interactively select a subdomain of a grid.
%%
%   [subdomain,sublimitxt] = gd_subdomain(grid); 

%%
% where _grid_ is a struct containing the gridded data, _npts_ is the
% number of points to be selected, _promptxt_ is the text used for the
% selection of each point, _isdel_ deletes the selection figure if true, 
% _points_ is an x, y struct of the selected points, _subdomain_ is an array of
% [min(x),max(x),min(y),max(y)] of the selected domain and _sublimitxt_ is
% a text description of the subdomain.

%%
% *gd_xy2sn*
% - map grid from cartesian to curvilinear coordinates with option to return 
% elevations on the source cartesian grid, or as a curvilinear grid.
%%
%   sngrid = gd_xy2sn(grid,cline,isxy,isplot);

%%
% *gd_sn2xy*
% - map grid from curvilinear to cartesian coordinates to return 
%   elevations on the source cartesian grid.
%%
%   sngrid = gd_sn2xy(grid,cline,isplot);

%%
% *gd_lineongrid_plot*
% -  plots a defined line onto a countour or surface plot of a grid (e.g a
%   channel centre-line). If z values are included the plot is a 3D plot.
%%
%   hfig = gd_lineongrid_plot(grid,cline,plotxt);

%%
% where _grid_ is the source grid, _cline_ is an x,y struct of centre-line
% path, _isxy_ is true to return cartesian grid and false to return
% curvilinear grid, _isplot_ is true to generate a plot of the grid and
% _plotxt_ is the title text for the plot.

%%
% *gd_user_function*
% - function for user to define bespoke use of grids and grid tools
%%
%   gd_user_functions(obj,mobj);       
%%
% where _obj_ is an instance of a class that inherits _GDinterface_ abstract class
% and _mobj_ is an instance of a model class (i.e. a class that inherits
% _muiModelUI_). _mobj_ enables access to other classes and data sets from the
% user function.
%%
% *getconvergencelength* 
% -  least squares fit using fminsearch to
% find the convergence length of a channel from a distance-width xy data set,
% where yi is the variable with respect to xi (e.g. area as a function of
% distance) and x0 is an initial guess used in fminsearch.
%%
%   convergencelength = getconvergencelength(xi,yi,x0);

%% 
% *getsubgrid*
% - extract a subdomain from a grid and return the extracted
% grid and the source grid indices of the bounding rectangle.
%%
%   [subgrid,ixo,iyo] = getsubgrid(grid,subdomain);
%%
% where _grid_ is an x,y,z struct of the grid and _subdomain_ defines the 
% bounding rectangle using [x0,xN,y0,yN], _subgrid_ is an x,y,z struct with the grid of
% the subdomain, _ixo_ and _iyo_ are the indices of the bounding rectangle in
% the form xo = [ix0,ix0,ixN,ixN,ix0] and iyo = [iyN,iy0,iy0,iyN,iyN].

%%
% *a_star*
% - implements the A* search algorithm to find the shortest path given
% constraints (inaccessible cells) and a cost function (e.g. water depths).
% Author: Alex Ranaldi, 2022, https://github.com/alexranaldi/A_STAR

%%
% *InterX* 
% - intersection of two curves. MATLAB Central File Exchange, 
% Author: NS, 2010, https://www.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections.

%%
% *xy2sn*
% -  Converts the Cartesian coordinates Xi and Yi 
% to the channel fitted coordinate Si and Ni using the centreline defined 
% through the Cartesian coordinates of its vertices Xc and Yc.
% MATLAB Central File Exchange, 
% Author: Bart Vermeulen,2022, Cartesian to Curvilinear coordinate forward and backward transformation.  
% https://www.mathworks.com/matlabcentral/fileexchange/55039-cartesian-to-curvilinear-coordinate-forward-and-backward-transformation 
%%
%   [Si,Ni]=xy2sn(Xc,Yc,Xi,Yi);

%%
% *sn2xy*
% - Converts the channel fitted coordinate Si 
% and Ni to the Cartesian coordinates Xi and Yi using the centreline 
% defined through the Cartesian coordinates of its vertices Xc and Yc.
% MATLAB Central File Exchange, 
% Author: Bart Vermeulen,2022, Cartesian to Curvilinear coordinate forward and backward transformation.  
% https://www.mathworks.com/matlabcentral/fileexchange/55039-cartesian-to-curvilinear-coordinate-forward-and-backward-transformation 
%%
%   [Xi,Yi]=sn2xy(Xc,Yc,Si,Ni);

%% See Also
% The grid classes make use of the <matlab:doc('gdinterface') GDinterface> 
% abstract class or link to classes that use interface properties.