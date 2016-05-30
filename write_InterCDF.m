%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITTEN BY ACHIM MORSCHHAUSER, GFZ POTSDAM, 2016
% mailto: mors/gfz-potsdam.de
%
% Write INTERMAGNET CDF file. The format of the file is based on the
% INTERMAGNET Discussion Document DD22, version 2.3, 23/01/2013.
%
% This function is based on the NASA cdf library for MATLAB which has to
% be added to the search path. This can be done by:
% addpath('Path to NASA CDF library')
% The internal MATLAB cdf library is not able to read 
% InterCDF files as the TT2000 date format is not supported with version
% R2015b or earlier.
%
% TODO: Additional user-specified data
%
% Arguments:
% ==========
%
% The arguments are given in their order.
%
% IN:
% ===
%
% A.) Mandatory fields
% --------------------
%
% NONE
%
% B.) Mandatory argument-value pairs:
% -----------------------------------
%
% The argument name and argument type are given, e.g. GeomagV, STRING, is 
% used as read_InterCDF('GeoMagV','XYZ'). Argument-value pairs can have
% arbitrary order and are not case-sensitive.
%
% GeoMag           ARRAY          Geomagnetic data for writing to CDF.
%                  OF DOUBLES     Different rows correspond to diff.
%                                 data types. Data types are specified in
%                                 the correct order in the field
%                                 'ElementsRecorded'.
%
%                  !!!! WARNING:  UNITS ARE AUTOMATICALLY SET TO 'NT' OR
%                  -------------  'DEGREES OF ARC'.
%
% ElementsRecorded STRING         Describes the recorded components, e.g. 
%                                 'XYZ'. Valid components are any
%                                 combination of 'XYZHDEVIFSG'.
%
% TimeV            ARRAY          Timestamps for geomagnetic vector data.
%                  OF TIME_FORMAT Mandatory only if scalar data are
%                  (see below)    specified in ElementsRecorded.
%
% TimeS            ARRAY          Timestamps for geomagnetic scalar data.
%                  OF TIME_FORMAT Mandatory only if vector data are
%                  (see below)    specified in ELementsRecorded.
%
% TimeFormat       STRING         Format of time. Admissible values are:
%                                 + TT2000    DEFAULT: Nanoseconds since
%                                             midday of 1st January 2000 AD
%                                 + Datenum   Fractional days since 
%                                             December 31,1 BC
%                                 + Datetime  MATLAB datetime object
%                                 + JD2000    Fractional days since
%                                             January 1st, 2000 AD
%                                 + POSIX     Number of seconds since
%                                             1st January, 1970 
%
% G_ATTR           STRUCT         The required Metadata as struct with
%                                 the following fields (case-sensitive!):
%                                 - IagaCode         STRING
%                                 - PublicationLevel STRING
%                                 - PublicationDate  TIME_FORMAT
%                                 - ObservatoryName  STRING
%                                 - Latitude         DOUBLE
%                                 - Longitude        DOUBLE
%                                 - Elevation        DOUBLE
%                                 - Institution      STRING
%                                 - Source           STRING
%                                 - SandardLevel     STRING
%                                 Must contain, if StandardLevel~='NONE':
%                                   - StandardName     STRING
%                                 Must contain, if StandardLevel=='PARTIAL'
%                                   - PartialStandDesc ARRAY OF DOUBLE 
%                                                      OR SINGLE STRING
%
%
% C.) Single optional fields
% --------------------------
%
% If an argument is given which is not listed here, it will be
% interpreted as an global attribute which is requested. The arguments 
% should be passed as strings, e.g. read_InterCDF('GeoMagS').
% Single optional fields can have an arbitrary order and are not
% case-sensitive.
%
% NoCheck            If specified, cdf file will not be checked for 
%                    correctness. Errors might result if cdf is not
%                    properly according to INTERMAGNET standard.
%
% D.) Optional argument-value pairs:
% ----------------------------------
%
% Temp               ARRAY          Temperature(s). Each temperature time
%                    OF DOUBLES     series corresponds to a row in Temp.
%                                   If temperature is specified, then 
%                                   Temp_Loc and Temp_Time need to be
%                                   specified as well.
%
% Temp_Loc           CELL ARRAY     Location where the respective temperature
%                    OF STRINGS     has been recorded.
%
% Temp_Time          ARRAY OF       Time for each temperature.
%                    TIME_FORMAT
%
% FILLVAL            ARRAY OF       The value used to show data is
%                    DOUBLES        missing. Default is 99999. Must be
%                                   specified for all data in the
%                                   order according to elementsrecorded, 
%                                   followed by the values for the
%                                   temperatures.
%
% VALIDMIN           ARRAY OF       As above, but for the smallest valid
%                    DOUBLES        value. Default is -79999, but -360 for 
%                                   'DI'.
%                                 
%
% VALIDMAX          ARRAY OF        As above, but for the largest valid
%                   DOUBLES         value. Default is 79999, but 360 for
%                                   'DI'.
%
% FILENAME          STRING          Filename with *.cdf extension. 
%                                   If not specified, an INTERMAGNET 
%                                   standard filename will be used.
%
% The following attributes can either be given as specified below, or
% may equally be added to G_Attr. If both are specified, the values as
% specified here will be used and warning will be displayed.
%
% VectorSensOrient   STRING
%
% StandardVersion    STRING
%
% TermsofUse         STRING
%
% UniqueIdentifier   STRING
%
% ParentIdentifiers  STRING
%
% ReferenceLinks     STRING
%
%
%
%
% EXAMPLE USAGE:
% ==============
%
% Write simple file without temperature and 1s sampling:
%
% cdf_file='test.cdf';
% M=spdfcdfinfo(cdf_file);
% [T D]=write_InterCDF('elementsrecorded','HDZ','geomag',20.*rand(3,10),...
%                      'timev',(1:10).*1e9,...
%                      'g_attr',M.GlobalAttributes,'timeformat','tt2000');
%
% Read CDF, change ElementsRecorded to XYZ, and save:
%
% cdf_file='test.cdf';
% R=read_InterCDF(cdf_file,'geomagv','XYZ','g_attr');
% D=[R.GeoMagV.X'; R.GeoMagV.Y'; R.GeoMagV.Z'];
% [T D]=write_InterCDF('filename','XYZ.cdf','elementsrecorded','XYZ',...
%                      'geomag',D,'timev',R.GeoMagV.Time,...
%                      'g_attr',R.g_attr,'timeformat','tt2000');
% 
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [M, D, VD]=write_InterCDF(varargin)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Preamble
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %
  % Some valid values definitions
  %

  valid_TimeFormats={'tt2000' 'datenum' 'datetime' 'jd2000' 'posix'};
  
  elements.valid='DEFGHISVXYZ';
  elements.units={'Degrees of arc' 'nT' 'nT' 'nT' 'nT' 'Degrees of arc' ...
      'nT'  'nT'  'nT'  'nT'  'nT'};
  
  %
  % Optional global attributes
  %
  IN.g_attr.opt={'VectorSensOrient' 'StandardVersion' 'TermsofUse' ...
      'UniqueIdentifier' 'ParentIdentifiers' 'ReferenceLinks'};
  
  %
  % Initialize flags for optional arguments
  %
  
  % Data
  IN.geomag.set=0;
  IN.elementsrecorded.set=0;
  IN.timev.set=0;
  IN.times.set=0;
  IN.timeformat.set=0;
  IN.g_attr.set=0;
  
  IN.nocheck.set=0;
  
  IN.temp.set=0;
  IN.temp_loc.set=0;
  IN.temp_time.set=0;
  IN.fillval.set=0;
  IN.validmin.set=0;
  IN.validmax.set=0;  
  IN.filename.set=0;  
  
  IN.vectorsensorient.set=0;
  IN.standardversion.set=0;
  IN.termsofuse.set=0;
  IN.uniqueidentifier.set=0;
  IN.parentidentifiers.set=0;
  IN.referencelinks.set=0;

  % Sampling rate // used for filename generation
  sampling.v=[];
  sampling.s=[];
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Process input arguments
  %  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %
  % Parse optional input arguments. Checking for correct format and
  % consistency is done later.
  %
  arg=1;
  while arg<nargin
      
      % Check if arguments are all strings
      if ~ischar(varargin{arg})
          
          error(['Argument ' num2str(arg+1) ' must be a string, but is '...
              class(varargin{arg}) '.' ])
          
      else

          switch lower(varargin{arg})
              case 'geomag'
                  get_input('geomag','double');
              case 'elementsrecorded'
                  get_input('elementsrecorded','char');
              case 'timev'
                  get_input('timev','numeric')
              case 'times'
                  get_input('times','numeric')
              case 'timeformat'
                  get_input('timeformat','char')
              case 'g_attr'
                   get_input('g_attr','struct');  
              case 'nocheck'
                  if (IN.(item).set)
                      warning(['Warning: ' item ' multiply specified. ' ...
                          'Using first occurence.']);
                  end
                   IN.nocheck.set=1;
                   IN.nocheck.data=varargin{arg};
              case 'temp'
                  get_input('temp','double');
              case 'temp_loc'
                  get_input('temp_loc','cell');
              case 'temp_time'
                  get_input('temp_time','numeric')
              case 'fillval'
                  get_input('fillval','double');
              case 'validmin'
                  get_input('validmin','double');
              case 'validmax'
                  get_input('validmax','double');
              case 'filename'
                  get_input('filename','char');
              otherwise
                  % Check for optional global attributes
                  if (any(strcmpi(IN.g_attr.opt,varargin{arg})))
                    get_input(lower(varargin{arg}),'char');
                  else
                    error(['Argument ' num2str(arg+1) ' not understood: '...
                      varargin{arg}])
                  end
          end
          
      end
      
      
      arg=arg+1;
      
  end
  
  %
  % Check if input was already specified. If not, set IN parameters
  % accordingly.
  %
  function get_input(item,type)
      
      if (IN.(item).set)
          warning(['Warning: ' item ' multiply specified. Using' ...
              ' first occurence.']);
      elseif isa(varargin{arg+1},type)
          IN.(item).data=varargin{arg+1};
          IN.(item).set=1;          
      else
          error(['Warning: ' item ' is of wrong type. Given:' ...
               class(varargin{arg+1}) ', but must be: ' type '.']);
      end
      
      arg=arg+1;
      
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Perform some basic consistency checks. Thorough testing is done later.
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Number of elements recorded
  if (IN.elementsrecorded.set)
    nor=length(IN.elementsrecorded.data);
  else
    error('ElementsRecorded was not specified.');
  end
  
  % Number of temperatures recorded
  if IN.temp.set
      not=size(IN.temp.data,2);
  else
      not=0;
  end
  
  % Number of overall variables
  nov=nor+not;
  
  %
  % 1) Check if geomagnetic data is provided.
  % 2) Check if number of elements recorded is consistent with number of data
  %    provided.
  %
  if ~(IN.geomag.set)
      error('No geomagnetic field data is provided. Nothing to write.');
  elseif size(IN.geomag.data,1)~=nor
      error('Number of ElementsRecorded and  rows in GeoMag must match.');
  end
  
  %
  % Check if recorded elements are valid elements
  %
  for i=1:nor
      if isempty(strfind(elements.valid,IN.elementsrecorded.data(i)))
          error(['Element ' IN.elementsrecorded.data(i) ' is not a '...
              'valid element.']);
      end
  end  
  
  %
  % Check if optional METADATA are uniquely specified
  %
  if IN.g_attr.set
      for i=1:length(IN.g_attr.opt)
          element=lower(IN.g_attr.opt{i});
          if IN.(element).set && ...
                  any(strcmp(fields(IN.g_attr.data),IN.g_attr.opt{i}))
              disp(['Warning: ' element ' explicitly specified and '...
                  'part of G_Attr! Using ' IN.(element).data ...
                  ' as value for ' element]);    
          end
      end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Set some default variable attributes:
  % - FILLVAL       99999
  % - VALIDMIN      DI: -360, T: -273.15, ELSE: -79999
  % - VALIDMAX      DI:  360, ELSE: 79999
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  %
  % Check if number of entries is correctly set.
  %
  function var_attr(item,default,nov)
        
      if (IN.(item).set)
          if (length(IN.(item).data) ~= nov)
              error(['Number of ' item ' entries provided does not ' ...
                  'match the number of elements in ElementsRecorded']);
          end
      else
          IN.(item).data=default*ones(1,nov);
      end
        
  end
  
  %
  % FILLVAL: Entry used for missing value
  %
  var_attr('fillval',99999,nov);
  
  %
  % VALDIMIN: Lowest valid number
  %
  var_attr('validmin',-79999,nov);
  if ~(IN.validmin.set)
    % Different default for DI
    pos=regexp(IN.elementsrecorded.data,'[DI]');
    if ~isempty(pos)
        IN.validmin.data(pos)=-360*ones(1,length(pos));
    end
    % Different default for temperatures
    if (not > 0)
        IN.validmin.data(nor+1:end)=-273.15*ones(1,not);
    end
  end
  
  %
  % VALDIMAX: Maximum valid number
  %
  var_attr('validmax',79999,nov);
  if ~(IN.validmax.set)
    % Different default for DI
    pos=regexp(IN.elementsrecorded.data,'[DI]');
    if ~isempty(pos)
        IN.validmax.data(pos)=360*ones(1,length(pos)); 
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Build data and variable attributes
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %
  % Fill what can be filled
  %
  V.FILLVAL=cell(nov,2);
  V.FILLVAL(:,2) = num2cell(IN.fillval.data);
  
  V.VALIDMIN=cell(nov,2);
  V.VALIDMIN(:,2) = num2cell(IN.validmin.data);
  
  V.VALIDMAX=cell(nov,2);
  V.VALIDMAX(:,2) = num2cell(IN.validmax.data);
  
  V.DISPLAY_TYPE = cell(nov,2);
  V.DISPLAY_TYPE(:,2) = repmat({'time_series'},nov,1);
  
  %
  % Variable Number Counter
  %
  varNUM=1; 
  
  %
  % GeomagneticField data
  %
 
  % Iterate over geomagnetic field elements
  for i=1:nor
            
      % Current element (e.g. X)
      element=IN.elementsrecorded.data(i);
      
      % Variable Name
      vname=['GeomagneticField' element];
      
      % Write data
      D{1,varNUM}=vname;
      D{2,varNUM}=IN.geomag.data(i,:);
      
      % Write variable attributes
      V.FIELDNAM{varNUM,1}=vname;
      V.FIELDNAM{varNUM,2}=['Geomagnetic Field Element ' element];
      
      V.UNITS{varNUM,1}=vname;
      V.UNITS{varNUM,2}=...
              elements.units{strfind(elements.valid,element)};
      
      V.FILLVAL{varNUM,1} = vname;
      
      V.VALIDMIN{varNUM,1} = vname;
      
      V.VALIDMAX{varNUM,1} = vname;
      
      V.DEPEND_0{varNUM,1} = vname;
      if (regexp(element,'[SG]')>0)
        V.DEPEND_0{varNUM,2} = 'GeomagneticScalarTimes';
      else
        V.DEPEND_0{varNUM,2} = 'GeomagneticVectorTimes';
      end
      
      V.DISPLAY_TYPE{varNUM,1} = vname;
      
      V.LABLAXIS{varNUM,1} = vname;
      V.LABLAXIS{varNUM,2} = element;
      
      VD{1,varNUM} = vname;
      VD{2,varNUM} = 'double';
      
      M.Variables{varNUM,1} = vname;
      
      varNUM=varNUM+1;
      
  end 
   
  %
  % Temperatures
  %
  
  % Iterate over given temperature fields
  for i=1:not
      
      % Variable name
      vname=['Temperature' num2str(i)];
      
      % Write data
      D{1,varNUM}=vname;
      D{2,varNUM}=IN.temp.data(i,:);
      
      % Write variable attributes
      V.FIELDNAM{varNUM,1}=vname;
      V.FIELDNAM{varNUM,2}=['Temperature ' IN.temp_loc.data{i}];
      
      V.UNITS{varNUM,1}=vname;
      V.UNITS{varNUM,2}='Celsius';
      
      V.FILLVAL{varNUM,1} = vname;
      
      V.VALIDMIN{varNUM,1} = vname;
      
      V.VALIDMAX{varNUM,1} = vname;
      
      V.DEPEND_0{varNUM,1} = vname;
      V.DEPEND_0{varNUM,2} = ['Temperature' num2str(i) 'Times'];
      
      V.DISPLAY_TYPE{varNUM,1} = vname;
      
      V.LABLAXIS{varNUM,1} = vname;
      V.LABLAXIS{varNUM,2} = '?';
      
      VD{1,varNUM} = vname;
      VD{2,varNUM} = 'double';
      
      M.Variables{varNUM,1} = vname;
      
      varNUM=varNUM+1;
      
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Write timestamps.
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %
  % Set time format
  %
  if ~IN.timeformat.set
      error('Time format must be set!');
  elseif ~any(strcmp(IN.timeformat.data,valid_TimeFormats))
      error('Time Format is not a valid format!');
  end
  
  % Convert to TT2000
  function Tc=convert_to_tt2000(T,TimeFormat)
          
    switch TimeFormat
        
        case 'datenum'
            Tc=spdfdatenumtott2000(T);
            
        case 'datetime'
            if (~isa(T,'datetime'))
               error('Input data was specified as datetime, but is not');
            else  
               tmp1=(T.Second-floor(T.Second))*1000;
               tmp2=(tmp1-floor(tmp1))*1000;
               tmp3=(tmp2-floor(tmp2))*1000;
               milli=floor(tmp1);
               mikro=floor(tmp2);
               nano=floor(tmp3);
               Tc=spdfcomputett2000([T.Year T.Month T.Day T.Hour ...
                   T.Minute floor(T.Second) milli mikro nano]);
            end
            
        case 'jd2000'
            Tc=spdfdatenumtott2000(T+730486);
            
        case 'posix'
            Tc=spdfdatenumtott2000((T+719529)/68400);
            
        case 'tt2000'
            Tc=int64(T);
            
        otherwise
            Tc=T;
            
    end
    
  end
  
  %
  % Geomagnetic vector data
  %
  if ~IN.timev.set && ...
     ~(isempty(regexp(IN.elementsrecorded.data,'[DEFHIVXYZ]','once')))
    error(['Vector time series TimeV mus be set when vector elements '...
        'are specified.']);
  elseif (IN.timev.set)
      vname='GeomagneticVectorTimes';
      D{1,varNUM}=vname;
      D{2,varNUM}=convert_to_tt2000(IN.timev.data,IN.timeformat.data);
      VD{1,varNUM} = vname;
      VD{2,varNUM} = 'tt2000';
      M.Variables{varNUM,1} = vname;
      sampling.v=varNUM;
      varNUM=varNUM+1;
  end
  
  %
  % Geomagnetic scalar data
  %
  if ~IN.times.set && ...
     ~(isempty(regexp(IN.elementsrecorded.data,'[GS]','once')))
    error(['Scalar time series TimeS mus be set when scalar elements '...
        'are specified.']);
  elseif (IN.times.set)
      vname='GeomagneticScalarTimes';
      D{1,varNUM}=vname;
      D{2,varNUM}=convert_to_tt2000(IN.times.data,IN.timeformat.data);
      VD{1,varNUM} = vname;
      VD{2,varNUM} = 'tt2000';
      M.Variables{varNUM,1} = vname;
      sampling.s=varNUM;
      varNUM=varNUM+1;
  end  
  
  %
  % Temperature data
  %
  if (IN.temp.set)
      
      if ~IN.temp_loc.set
          error('Temp_Loc is missing.');
      elseif (length(IN.temp_loc.data)~=size(IN.temp.data,1))
          error('Temperature location for some temperatures missing.');
      end
      
      if ~IN.temp_time.set
          error('Temp_Time is missing.');
      elseif (size(IN.temp.data,1)~=size(IN.temp_time.data,1))
          error('Time series for some temperatures missing.');
      else
          for i=1:not
              vname=['Temperature' num2str(i) 'Times'];
            D{1,varNUM}=vname;
            D{2,varNUM}=convert_to_tt2000(IN.temp.data(i,:),...
                IN.timeformat.data);
            VD{1,varNUM} = vname;
            VD{2,varNUM} = 'tt2000';
            M.Variables{varNUM,1} = vname;
            varNUM=varNUM+1;
          end
      end
      
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  %
  % Build GlobalAttributes
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  if (IN.g_attr.set && ~isempty(IN.g_attr.data))
      M.GlobalAttributes=IN.g_attr.data;
  elseif isempty(IN.g_attr.data)
      error('Global Attributes (G_ATTR) is empty!');
  else
      error('Global Attributes (G_ATTR) must be specified!');
  end
  % Write elementsrecorded
  M.GlobalAttributes.ElementsRecorded{1}=IN.elementsrecorded.data;
  % Write optional elements
  for i=IN.g_attr.opt
      if (IN.(lower(i{1})).set)
          tmp{1}=IN.(lower(i{1})).data;
          M.GlobalAttributes.(i{1})=tmp;
      end
  end
  
  
  % VariableAttributes
  M.VariableAttributes=V;

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Check consistency and validity of data provided.
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  
  
  if (~IN.nocheck.set || IN.nocheck.data)
    M.Variables(:,4)=VD(2,:);
    [valid, desc]=check_InterCDF(M,D(2,:));
    % Display any errors / warnings
    fprintf(desc);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Write data.
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  if valid>0
      %
      % AUTOMATICALLY GENERATE FILENAME
      %
      if ~(IN.filename.set)
          % Determine sampling period
          if ~isempty(sampling.v) && length(D{2,sampling.v})>1
              starttime=spdfbreakdowntt2000(D{2,sampling.v}(1));
              % 1s sampling
              if diff(D{2,sampling.v}(1:2))==1e9
                  datetime=num2str(starttime(1:6),...
                      '%04u%02u%02u_%02u%02u%02u');
                  cad='PT1S';
              % 1m sampling
              elseif diff(D{2,sampling.v}(1:2))==1e9*60
                  datetime=num2str(starttime(1:5),...
                      '%04u%02u%02u_%02u%02u');
                  cad='PT1M';
              % 1h sampling
              elseif diff(D{2,sampling.v}(1:2))==1e9*60*60
                  datetime=num2str(starttime(1:4),'%04u%02u%02u_%02u');
                  cad='PT1H';
              % 1d sampling
              elseif diff(D{2,sampling.v}(1:2))==1e9*60*60*24
                  datetime=num2str(starttime(1:3),'%04u%02u%02u');
                  cad='P1D';
              % 1 month sampling
              elseif diff(D{2,sampling.v}(1:2))<=1e9*60*60*32 && ...
                      diff(D{2,sampling.v}(1:2))>=1e9*60*60*27
                  datetime=num2str(starttime(1:2),'%04u%02u');
                  cad='P1M';
              elseif diff(D{2,sampling.v}(1:2))<=1e9*60*60*368 && ...
                      diff(D{2,sampling.v}(1:2))>=1e9*60*60*365
                  datetime=num2str(starttime(1),'%04u');
                  cad='P1Y';
              else
                  datetime=num2str(starttime(1:6),...
                      '%04u%02u%02u_%02u%02u%02u');
                  cad='UNKNOWN';
              end
          elseif isempty(sampling.v)
              error(['Cannot generate filename as no geomagnetic vector'...
                  ' data was specified!']);
          else
              error(['Cannot generate filename as sampling period'...
                  ' cannot be determined from single datum.']);
          end
          
          IN.filename.data=[M.GlobalAttributes.IagaCode{1} '_' ...
              datetime '_' cad '_' ...
              M.GlobalAttributes.PublicationLevel{1} '.cdf'];
      end
      
      %
      % WRITE DATA
      %
      spdfcdfwrite(IN.filename.data,D,...
          'VariableAttributes',M.VariableAttributes,'Vardatatypes',VD,...
          'GlobalAttributes',M.GlobalAttributes);
      
      disp(['CDF successfully written to ' IN.filename.data]);
      
  else
      
      error('Cannot write file. Incomplete data!\n')
      
  end
  
  
  
end