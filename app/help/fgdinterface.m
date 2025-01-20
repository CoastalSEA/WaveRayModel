%% FGDinterface
% Abstract class with properties and methods to manipulate gridded datasets 
% for use in applications that import data, or save model ouput. The class 
% inherits <matlab:doc('muidataset') muiDataSet> and <matlab:doc('gdinterface') GDinterface>. Together they provide 
% an extensive set of methods to handle datasets of various types (e.g. 
% from models or imported files).

%% Description
% *FGDinterface* is used as a superclass to provide data handling 
% functionality for classes that import different types of data or models
% that are gridded, in order to store them in a consistent and documented format.
% The data are stored in a set of <matlab:doc('dstable') dstables>, which
% are assigned as a <matlab:doc('struct') struct> to the <matlab:doc('muicatalogue') muiCatalogue>
% _Data_ property (typically as part of a UI based on the <matlab:doc('muitoolbox') muitoolbox>). 
% The core table is the _'Grid'_ <matlab:doc('dstable') dstable> which holds the
% gridded data and supporting metadata - see <matlab:doc('gdinterface')
% GDinterface> documentation.
%%
% A number of additional tables support applications that represent
% inlet, estuary, river and valley formations. These include _'Hypsometry'_ for the
% variation of surface area and volume as a function of elevation,
% _'SectionProps'_ holds the variation of the width and cross-sectional area
% along the x-axis derived from the gridded data, _'Plan'_ holds the widths along 
% the x-axis specified in the model, or based on an exponential fit to the
% gridded data, _'WaterLevels'_ holds the variation in water levels at high
% water, mean tide and low water along the x-axis, and _'GrossProps'_ holds
% various summary properties of the form related to length, area, volume,
% rate of convergence, etc. The metatdata for the dstable content are
% defined in the _*setDSproperties*_ method in *FGDinterface*. 
%% 
% <html>
% <table border=1><tr><td><u>Note</u>: In UI's based on <a href="matlab:doc('muitoolbox')">muitoolbox</a>
% the <i>DSproperties</i> can be viewed by using the mouse button callback 
% in the first column of a Case on the Cases/Data/Model tab in the main interface.
% </td></tr></table>
% </html>

%% FGDinterface properties
% The class uses the properties defined for <matlab:doc('muidataset') muiDataSet> 

%% FGDinterface methods
% The methods in the abstract class include the following:
%%
% Methods to add, delete and update grids and propetry tables:

%%
% *setGrid* - load the grid results into a dstable, where _obj_ is any 
% GDintereface subclass; _gridddata_ is a cell array with a matrix of
% elevations, and any other data to be included in the dstable; _dims_ is a
% struct for the dimensions of griddata, where the fields 'x' and 'y' are the 
% vector dimensions of z and must be unique, and the field 't' can be 
% datetimes, durations or values to assign to rows in dstable, 'ishead' is
% a logical flag defining the x-axis orientation, 'xM' is the distance to 
% the inlet/estuary/valley mouth from the grid origin, and 'cline' is a struct 
% with the x and y coordinates of the centre-line; and 
% _meta_ is a stuct with 'source' and 'data' fields to describe file name or model
% used, and additional run details, respectively.
%%
%   obj = setGrid(obj,griddata,dims,meta);

%%
% *getGrid* - retrieve a grid, where _obj_ is any GDintereface subclass;
% and _irow_ is the row index of the dstable to extract grid from
% (optional). Returns grid as a struct with the fields, 'x', 'y', 'z' and
% 't, together with any metadata that may be assigned to the dstable,
% including, 'irow', the row index in the dstable, 'desc', the dstable Description, 'ishead',
% a logical flag defining the x-axis orientation, 'xM', the distance to 
% the inlet/estuary/valley mouth from the grid origin, and 'cline', a struct 
% with the x and y coordinates of the centre-line.
%%
%   grid = getGrid(obj,irow);

%%
% *addGrid* - add a grid to an existing <matlab:doc('muitoolbox') muitoolbox> 
% Case table, where _obj_ is any GDintereface subclass; _muicat_ is a 
% handle to the <matlab:doc('muicatalogue') muiCatalogue>; _newgrid_ is the
% grid to add to the dstable and should be of the same dimensions as
% existing grids in the table; _timestep_ is the value to use for RowName
% in the dstable; _dims_ is a struct containig xM, Lt and river properties; 
% _sourcetxt_ is the file name, or model name, to identify 
% data source; and _ismsg_ is a flag set to true to use the defualt message, 
% and false to suppress the message.
%%
%   addGrid(obj,muicat,newgrid,timestep,dims,sourcetxt,ismsg);

%%
% *deleteGrid* - delete a grid to an existing <matlab:doc('muitoolbox') muitoolbox>
% Case table, where _obj_ is any FGDintereface subclass; _classrec_ is the
% class record number; _catrec_ is the catalogue record for the instance; and 
% _muicat_ is a handle to the <matlab:doc('muicatalogue') muiCatalogue>.
%%
%   deleteGrid(obj,classrec,catrec,muicat);

%%
% *setProperties* - initialises the tables for derived properties of gridded data set

%%
%    obj = setProperties(obj,zwl,Wz,limits,histint);

%%
% *addProperties* - add properties to an existing set of grid property tables
%%
%   obj = addProperties(obj,irow,zwl,Wz,limits,histint);

%%
% where _obj_ is any class that contains a Grid dstable; _irow_ is the row of 
% Grid table to use; _zwl_ is a struct or 
% table, of water levels at zhw,zmt,zlw; _Wz_ is a struct or table of 
% width plan form data at 3 elevations' _limits_ is an array for the upper 
% and lower limits for vertical range [0=use grid; 1=prompt user; [x,y]=limits to use];         
% histint is the vertical interval for hypsometry (optional user prompted
% if not defined).

%%
% *delProperties* - delete a row of properties from existing grid property
% tables, where _obj_ is any class that contains a Grid dstable; and
% _row2use_ is row(s) to be deleted from ALL property tables
%%
%   [obj,isok] = delProperties(obj,row2use);

%% 
% *addFormProps* - add form properties to a gridded data set (needs to be a
% suitable dtm, ie inlet, estuary, valley), where _obj_ is any FGDintereface subclass; 
% and _muicat_ is the handle to the <matlab:doc('muicatalogue') muiCatalogue>.
%%
%   addFormProps(obj,muicat);

%%
% *delFormProperties* - delete ALL property tables associated with a 
% selected gridded data set, where _obj_ is any FGDintereface subclass; and 
% _muicat_ is a handle to the <matlab:doc('muicatalogue') muiCatalogue>.
%%
%   delFormProperties(obj,muicat);

%%
% *setModelFormProps* -add a set of Hypsometry, Section and Gross form properties 
% (this function is used from models such as CF_TransModel)
% where _obj_ is any FGDintereface subclass.
%%
%   obj = setModelFormProps(obj);

%%
% *editGridInletData* - UI to edit, definition of channel head, x-distance 
% to mouth and definition of centre-line (if used).
%%
%   editGridInletData(obj);


%%
% Methods to manipulate grids use the GDinterface class  - see <matlab:doc('gdinterface')
% GDinterface> documentation.


%% FGDinterface static methods
% The additional static methods in the abstract class include the following:


%%
% *addValleyBase* - use an xyz definition of the channel thalweg to create 
% a valley base this is used to modify an existing valley defined above MHW 
% to include the valley form down to the pre-Holocene surface (i.e. below the
% existing channel bed), where _muicat_ is a handle to the <matlab:doc('muicatalogue') muiCatalogue>; 
% _gridclasses_ is a cell array of clases to be used for case selection.
%%
%   FGDinterface.addValleyBase(muicat,gridclasses);

%%
% *addShore* - add a strip to shore side of a grid by extrapolating from the
% neighbouring strip, where _muicat_ is a handle to the <matlab:doc('muicatalogue') muiCatalogue>; 
% _gridclasses_ is a cell array of clases to be used for case selection.
%%
%   FGDinterface.addShore(muicat,gridclasses);

%%
% *formMenuOptions* - handles the menu calls for the default menu of Form based Grid
% Tools, where _mobj_ is a mui model instance, _src_ is a handle to the calling 
% menu option, and _gridclasses_ is a cell array of classes to use in the case
% selection.
%%
%   FGDinterface.formMenuOptions(mobj,src,gridclasses);


%% See Also
% Some classes that use *FGDinterface* and related functions are detailed in
% the help for <matlab:doc('grid_class_fcns') Grid classes and functions>.