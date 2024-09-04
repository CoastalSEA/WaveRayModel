%% GDinterface
% Abstract class with properties and methods to manipulate gridded datasets 
% for use in applications that import data, or save model ouput. The class 
% inherits <matlab:doc('muidataset') muiDataSet> and together they provide 
% an extensive set of methods to handle datasets of various types (e.g. 
% from models or imported files).

%% Description
% *GDinterface* is used as a superclass to provide data handling 
% functionality for classes that import different types of data or models
% that are gridded, in order to store them in a consistent and documented format.
% The data are stored in a set of <matlab:doc('dstable') dstables>, which
% are assigned as a <matlab:doc('struct') struct> to the <matlab:doc('muicatalogue') muiCatalogue>
% _Data_ property (typically as part of a UI based on the <matlab:doc('muitoolbox') muitoolbox>). 
% The core table is the _'Grid'_ <matlab:doc('dstable') dstable> which holds the
% gridded data and supporting metadata. The metatdata for the dstable content are
% defined in the _*setDSproperties*_ method in *GDinterface*. 

%% 
% <html>
% <table border=1><tr><td><u>Note</u>: In UI's based on <a href="matlab:doc('muitoolbox')">muitoolbox</a>
% the <i>DSproperties</i> can be viewed by using the mouse button callback 
% in the first column of a Case on the Cases/Data/Model tab in the main interface.
% </td></tr></table>
% </html>

%% GDinterface properties
% The class uses the properties defined for <matlab:doc('muidataset') muiDataSet> 

%% GDinterface methods
% The methods in the abstract class include the following:
%%
% Methods to add, delete and update grids:

%%
% *setGrid* - load the grid results into a dstable, where _obj_ is any 
% GDintereface subclass; _gridddata_ is a cell array with a matrix of
% elevations, and any other data to be included in the dstable; _dims_ is a
% struct for the dimensions of griddata, where the fields 'x' and 'y' are the 
% vector dimensions of z and must be unique, and the field 't' can be 
% datetimes, durations or values to assign to rows in dstable; and 
% _meta_ is a stuct with 'source' and 'data' fields to describe file name or model
% used, and additional run details, respectively.
%%
%   obj = setGrid(obj,griddata,dims,meta);

%%
% *getGrid* - retrieve a grid, where _obj_ is any GDintereface subclass;
% and _irow_ is the row index of the dstable to extract grid from
% (optional). Returns grid as a struct with the fields, 'x', 'y', 'z' and
% 't, together with any metadata that may be assigned to the dstable,
% including, 'irow', the row indexin the dstable, 'desc', the dstable Description, and 'cline', a struct 
% with the x and y coordinates of the centre-line if a meander has been added.
%%
%   grid = getGrid(obj,irow);

%%
% *addGrid* - add a grid to an existing <matlab:doc('muitoolbox') muitoolbox> 
% Case table, where _obj_ is any GDintereface subclass; _muicat_ is a 
% handle to the <matlab:doc('muicatalogue') muiCatalogue>; _newgrid_ is the
% grid to add to the dstable and should be of the same dimensions as
% existing grids in the table; _timestep_ is the value to use for RowName
% in the dstable; _sourcetxt_ is the file name, or model name, to identify 
% data source; and _ismsg_ is a flag set to true to use the defualt message, 
% and false to suppress the message.
%%
%   addGrid(obj,muicat,newgrid,timestep,[],sourcetxt,ismsg);

%%
% *deleteGrid* - delete a grid to an existing <matlab:doc('muitoolbox') muitoolbox>
% Case table, where _obj_ is any GDintereface subclass; _classrec_ is the
% class record number; _catrec_ is the catalogue record for the instance; and 
% _muicat_ is a handle to the <matlab:doc('muicatalogue') muiCatalogue>.
%%
%   deleteGrid(obj,classrec,catrec,muicat);

%%
% Methods to manipulate grids:
%%
% *translateGrid* - interactively translate grid x-y coordinates, 
% where _obj_ is any GDintereface subclass; _classsrec_ is the class record
% index of selected instance and _muicat_ is a handle to the
% <matlab:doc('muicatalogue') muiCatalogue>.
%%
%   translateGrid(obj,classrec,muicat);

%%
% *rotateGrid* - interactively flip or rotate grid,
% where _obj_ is any GDintereface subclass; _classsrec_ is the class record
% index of selected instance and _muicat_ is a handle to the
% <matlab:doc('muicatalogue') muiCatalogue>.
%%
%   rotateGrid(obj,classrec,muicat)

%%
% *reGridData* - regrid a gridded dataset to match another grid or to user
% specified dimensions, where _obj_ is any GDintereface subclass; _mobj_ 
% is mui model instance; and _gridclasses_ is a cell array of clases to be
% used for case selection.
%%
%   reGridData(obj,mobj,gridclasses);

%%
% * *subGridData* - interactively define a subgrid and save grid as a new
% Case, where _obj_ is any GDintereface subclass; _classsrec_ is the class record
% index of selected instance and _muicat_ is a handle to the
% <matlab:doc('muicatalogue') muiCatalogue>.
%%
%   subGridData(obj,muicat);

%%
% *addSurface* - add horizontal surface to an extisting grid,
% where _obj_ is any GDintereface subclass; _classsrec_ is the class record
% index of selected instance and _muicat_ is a handle to the
% <matlab:doc('muicatalogue') muiCatalogue>.
%%
%   addSurface(obj,muicat);

%%
% *curvilinear_xy2sn* - map grid from cartesian to curvilinear
% coordinates, where _obj_ is any GDintereface subclass; _classsrec_ is the class record
% index of selected instance and _muicat_ is a handle to the
% <matlab:doc('muicatalogue') muiCatalogue>.
%%
%   curvilinear_xy2sn(cobj,muicat);

%%
% *curvilinear_sn2xy* - map grid from curvilinear to cartesian coordinates,
% where _obj_ is any GDintereface subclass; _classsrec_ is the class record
% index of selected instance and _muicat_ is a handle to the
% <matlab:doc('muicatalogue') muiCatalogue>.
%%
%   curvilinear_sn2xy(cobj,muicat);

%%
% *exportGrid* - select a Case and export grid as xyz tuples,
% where _obj_ is any GDintereface subclass.
%%
%   exportGrid(obj);

%% 
% *displayGridDims* - display the dimensions of a selected grid in a table
% figure, where _obj_ is any GDintereface subclass.
%%
%   displayGridDims(obj)

%% GDinterface static methods
% The static methods in the abstract class include the following:

%%
% *setCombinedGrids* - superimpose one grid on another based on maximum
% or minimum set of values
%%
%   GDinterface.setCombinedGrids(muicat,gridclasses,ismax,prompts); %see below for detail of input variables

%%
% *diffGridsPlot* - generate a plot of the difference between two grids
%%
%   GDinterface.diffGridsPlot(muicat,gridclasses); %see below for detail of input variables

%%
% *plotSections* - display grid and allow user to interactively define start and
% end points of a section line to be plotted in a figure.
%%
%   GDinterface.plotSections(muicat,gridclasses); %see below for detail of input variables

%%
% *getGridLine* - interactively digitise a line and save to a file
%%
%   GDinterface.getGridLine(muicat,gridclasses); %see below for detail of input variables


%%
% where _muicat_ is a handle to the <matlab:doc('muicatalogue') muiCatalogue>; 
% _gridclasses_ is a cell array of clases to be used for case selection; 
% _ismax_ is logical true, then use maximum value at each point, if false 
% then use minimum values (optional); and _prompts_ is a struct, with fields 
% 'first' and 'second' that define prompts for the selection of the two 
% Cases to be combined (optional).

%%
% *gridMenuOptions* - handles the menu calls for the default menu of Grid
% Tools. 
%%
%   GDinterface.gridMenuOptions(mobj,src,gridclasses);
%%
% where _mobj_ is a mui model instance, _src_ is a handle to the calling 
% menu option, and _gridclasses_ is a cell array of classes to use in the case
% selection.

%% See Also
% Some classes that use *GDinterface* and related functions are detailed in
% the help for <matlab:doc('grid_class_fcns') Grid classes and functions>.