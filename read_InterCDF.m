%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITTEN BY ACHIM MORSCHHAUSER, GFZ POTSDAM, 2016
% mailto: mors/gfz-potsdam.de
%
%
% This function is based on the NASA cdf library for MATLAB which has to
% be added to the search path. This can be done by:
% addpath('Path to NASA CDF library')
% The internal MATLAB cdf library is not able to read 
% InterCDF files as the TT2000 date format is not supported with version
% R2015b or earlier.
%
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
% STRING                        cdf filename
%
% B.) Optional argument-value pairs:
% ----------------------------------
%
% The argument name and argument type are given, e.g. GeomagV, STRING, is 
% used as read_InterCDF('GeoMagV','XYZ'). Argument-value pairs can have
% arbitrary order and are not case-sensitive.
%
% GeoMagV      STRING           Get the actual data and time of the
%                               specified components which can be
%                               any combination of 'XYZHDEVI'
%                               and 'R' (Recorded)
%
% Temperature  INT ARRAY        Temperature(s) with given number(s) and 
%                               coresponding time
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
% GeoMagS             Get the recorded scalar values including time
%                     information, if scalar values are available.
%
% DataFormat          Get the attributes related to the DataFromat:
%                     - FormatDescription
%                     - FormatVersion
%                     - Title
%
% DataIdent           Get the attributes related to uniquely identifieng 
%                     the data:
%                     - IagaCode           Output IAGA code for observatory
%                     - ElementsRecorded   Recorded Elements, e.g. 'HDZ'
%                     - PublicationLevel   PublicationLevel
%                     - PublicationDate    Date of publication
%
% Observatory        Observatory Attributes
%                    - ObservatoryName    Full name of the observatory
%                    - Latitude           Latitude of the observatory
%                    - Longitude          Longitude of the observatory
%                    - Elevation          Elevation of the observatory
%                    - Institution        Institution
%                    - VectorSensOrient   See INTERMAGNET DD22
%
% DataQuality        Data Quality Attributes
%
% DataPublication    Data Publication Attributes
%
% NoCheck            If specified, cdf file will not be checked for 
%                    correctness. Errors might result if cdf is not properly
%                    according to INTERMAGNET standard
%
% Meta_All           All metadata is returned.
%
% G_Attr             Metadata of global attributes.
%
% TimeFormat         Format of time. Admissible values are:
%                    - TT2000      Nanoseconds since midday of 1st 
%                                  January 2000 AD
%                    - Datenum     Fractional days since December 31,1 BC
%                    - Datetime    MATLAB datetime object
%                    - JD2000      Fractional days since January
%                                  1st, 2000 AD
%                    - POSIX      Number of seconds since 1st January, 1970
%
% OUT:
% ====
%
% The output is arranged in a struct with the requested fields.
%
% EXAMPLE USAGE:
% ==============
%
% cdf_file='test.cdf';
%
% Read geomagnetic field data as stored in CDF. TimeFormat will be TT2000.
% R=read_InterCDF(cdf_file,'geomagv','R');
%
% Read geomagnetic field data as stored in CDF. TimeFormat will be JD2000.
% R=read_InterCDF(cdf_file,'geomagv','R','timeformat','JD2000');
%
% Read first temperature recording
% R=read_InterCDF(cdf_file,'temperature',1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function struct_out=read_InterCDF(filename,varargin)
 
  tic % Time

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Some constants
  % 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %
  % Define attribute classes (INTERMAGNET DD22, Sec. 4)
  %
  
  GAttr.classes={'DataFormat' 'DataIdent' 'Observatory' 'DataQuality'...
      'DataPublication'};
  
  % Elements of attribute classes
  GAttr.DataFormat.class={'FormatDescription' 'FormatVersion' ...
      'Title'};
  GAttr.DataIdent.class={'IagaCode' 'ElementsRecorded' ...
      'PublicationLevel' 'PublicationDate'};
  GAttr.Observatory.class={'ObservatoryName' 'Latitude' 'Longitude' ...
      'Elevation' 'Institution' 'VectorSensOrient'};
  GAttr.DataQuality.class={'StandardLevel' 'StandardName' 'StandardVersion' ...
      'PartialStandDesc'};
  GAttr.DataPublication.class={'Source' 'TermsOfUse' 'UniqueIdentifier' ...
      'ParentIdentifiers' 'ReferenceLinks'};
  
  %
  % Define allowed time formats
  %
  TimeFormats={'tt2000' 'datenum' 'datetime' 'jd2000' 'posix'};
  
  %
  % Initialize flags for optional arguments
  %
  
  % Data
  GeoMagV.set=0;
  GeoMagS.set=0;
  
  % Global Attribute classes
  GAttr.DataFormat.set=0;       % Data Format (Sec. 4.1)
  GAttr.DataIdent.set=0;        % Data Identifiers
  GAttr.Observatory.set=0;      % Observatory Identifiers
  GAttr.DataQuality.set=0;      % Data Quality & Standards
  GAttr.DataPublication.set=0;  % Data Publication
  
  % User Defined Attributes
  UAttr.set=0;
  UAttr.req=cell(0,0);
  
  % Temperature
  Temp.set=0;
  Temp.req=[];
  
  % NoCheck
  NoCheck.set=0;
  
  % TimeFormat
  TimeFormat.req='tt2000';
  TimeFormat.set=0;
  
  % Metadata
  Meta.All.set=0;
  Meta.G.set=0;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Process input arguments
  %  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %
  % Read filename
  %
  if (~ischar(filename))
    error('First argument (filename) must be a string!');
  end
  
  %
  % Parse optional input arguments
  %
  arg=1;
  while arg<nargin

    % Check if argumetns are all strings
    if ~ischar(varargin{arg})
        
        error(['Argument ' num2str(arg+1) ' must be a string, but is '...
               class(varargin{arg}) '.' ])
        
    else
        
      switch lower(varargin{arg})
        case 'geomagv'
          if (GeoMagV.set)
            disp(['Warning: GeoMagV multiply specified. Using ' ...
                GeoMagV.req '.']);
          elseif ~isa(varargin{arg+1},'char')
              error('GeomagV must be of type char.');
          else
            GeoMagV.req=varargin{arg+1};
            GeoMagV.set=1;
          end
          arg=arg+1;
        case 'geomags'
          GeoMagS.set=1;
        case 'temperature'
          Temp.req=floor(varargin{arg+1});
          if ~isa(Temp.req,'numeric')
              error(['Temperature must be of type numeric.']);
          end
          Temp.set=1;
          arg=arg+1;
        case 'dataident'
          GAttr.DataIdent.set=1;
        case 'observatory'
          GAttr.Observatory.set=1;          
        case 'dataquality'
          GAttr.DataQuality.set=1;            
        case 'datapublication'
          GAttr.DataPublication.set=1;
        case 'nocheck'
          NoCheck.set=1;
        case 'timeformat'
          if (TimeFormat.set)
            disp(['Warning: TimeFormat multiply specified. Using ' ...
                TimeFormat.req '.']);
          elseif ~isa(varargin{arg+1},'char')
              error('TimeFormat mus be of type char.')
          else
            % Check if time format is valid
            if (any(strcmpi(varargin{arg+1},TimeFormats)))
                TimeFormat.req=lower(varargin{arg+1});
                TimeFormat.set=1;
            else
                disp('Warning: Specified time format is not valid.');
            end
          end
          arg=arg+1;
        case 'meta_all'
            Meta.All.set=1;
          case 'g_attr'
            Meta.G.set=1;
        otherwise
          UAttr.set=1;
          UAttr.req{end+1}=varargin{arg};
          %error(['Argument ' num2str(arg+1) ' not understood: '...
          %       varargin{arg}])
      end
      
    end
    
    arg=arg+1;
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Read (meta)data and check correctness of file
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
  
  % Read (meta)data
  try
      M=spdfcdfinfo(filename);
      D=spdfcdfread(filename,'KeepEpochAsIs',1);
  catch
      error(['CDF cannot be read. File may be corrupt '...
          'or filename may have invalid characters'])
  end
  
  % Check CDF file
  if (~NoCheck.set)
      [valid, desc]=check_InterCDF(M,D);
      if (valid<0)
          fprintf(desc);
          error(['File ' filename ...
              ' is not in valid INTERMAGNET CDF format!']);
      elseif (valid==0)
          disp(['File ' filename ...
              ' is not in valid INTERMAGNET CDF format!']);
          disp(['However, it will still be read' ...
              ' with the following warnings: ']);
          fprintf(desc);
      end
  end
  
  toc
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Generate output data
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %
  % Meta_All
  %
  if (Meta.All.set)
      struct_out.Meta_All=M;
  end
  
  %
  % Meta_G
  %
  if (Meta.G.set)
      struct_out.g_attr=M.GlobalAttributes;
  end
  
  %
  % GeoMagV
  %
  if (GeoMagV.set)
  
  % Recorded and requested field elements
  if (isfield(M.GlobalAttributes,'ElementsRecorded'))
    Recorded=M.GlobalAttributes.ElementsRecorded;
    Requested=strrep(GeoMagV.req,'R',M.GlobalAttributes.ElementsRecorded{1});
  else
    Recorded='';
    Requested=strrep(GeoMagV.req,'R','');
  end

  
  % Get vector data
  struct_out.GeoMagV=[];
  for i=1:length(Requested)
      % Check if the requested component is recorded
      pos=getPosGF(Requested(i));
      if (pos>0)
          struct_out.GeoMagV.(Requested(i))=D{pos};
      % If the requested component is not recorded, try to calculate it
      elseif ~isfield(struct_out.GeoMagV,upper(Requested(i)))
          switch upper(Requested(i))
              case 'X'
                  posH=getPosGF('H');
                  posD=getPosGF('D');
                  posI=getPosGF('I');
                  posSF=getPosGF('SF');
                  posY=getPosGF('Y');
                  if (posD>0 && posH>0)
                          struct_out.GeoMagV.X=cosd(D{posD}).*D{posH};
                  elseif (posD>0 && posSF>0 && posI>0)
                          struct_out.GeoMagV.X=cosd(D{posD}).*...
                                  (cosd(D{posI}).*D{posSF});
                  elseif (posD>0 && posY>0)
                      struct_out.GeoMagV.X=tand(D{posD}).*D{posY};
                      % For D=+/-90, X will be infinite, but should equal H
                      % However, H is not available
                      struct_out.GeoMagV.X(isinf(struct_out.GeoMagV.X))=NaN;
                  else
                      error(['Vector Element ' element ' cannot be '...
                          ' calculated from ' Recorded '.']);
                  end
              case 'Y'
                  posH=getPosGF('H');
                  posD=getPosGF('D');
                  posI=getPosGF('I');
                  posSF=getPosGF('SF');
                  if (posD>0 && posH>0)
                          struct_out.GeoMagV.Y=sind(D{posD}).*D{posH};
                  elseif (posD>0 && posSF>0 && posI>0)
                          struct_out.GeoMagV.Y=sind(D{posD}).*...
                                  (cosd(D{posI}).*D{posSF});
                  elseif (posD>0 && posX>0)
                      struct_out.GeoMagV.Y=D{posX}./tand(D{posD});
                      % For D=0, Y will be infinite, but should equal H
                      % However, H is not available
                      struct_out.GeoMagV.Y(isinf(struct_out.GeoMagV.Y))=NaN;
                  else
                      error(['Vector Element ' element ' cannot be '...
                          ' calculated from ' Recorded '.']);
                  end
              case 'Z'
                  H=getH;
                  posI=getPosGF('I');
                  posSF=getPosGF('SF');
                  if (posI>0 && posSF>0)
                      struct_out.GeoMagV.Z=D{posF}.*sind(D{posI});
                  elseif (~isempty(H) && posSF>0)
                      struct_out.GeoMagV.Z=sqrt(D{posSF}.^2-H.^2);
                  elseif (posI>0 && ~isempty(H))
                      struct_out.GeoMagV.Z=H.*tand(D{posI});
                      % For I=+/-90, Z will be infinite, but should equal F
                      % However, F is not available
                      struct_out.GeoMagV.Z(isinf(struct_out.GeoMagV.Z))=NaN;
                  else
                      error(['Vector Element ' element ' cannot be '...
                          ' calculated from ' Recorded '.']);
                  end
                  clearvars H;
              case 'H'
                  H=getH;
                  if (~isempty(H))
                      struct_out.GeoMagV.H=H;
                  else
                      error(['Vector Element ' element ' cannot be '...
                          ' calculated from ' Recorded '.']);
                  end
              case 'D'
                  posX=getPosGF('X');
                  posY=getPosGF('Y');
                  H=getH;
                  if (posX>0 && posY>0)
                      struct_out.GeoMagV.D=atand(D{posX}./D{posY});
                  elseif (posX>0 && ~isempty(H))
                      struct_out.GeoMagV.D=asind(D{posX}./H);
                      % If H is zero, D is not defined
                      struct_out.GeoMagV.D(H==0)=NaN;
                  elseif (posY>0 && ~isempty(H))
                      struct_out.GeoMagV.D=acosd(D{posY}./H);
                      % If H is zero, D is not defined
                      struct_out.GeoMagV.D(H==0)=NaN;
                  else
                     error(['Vector Element ' element ' cannot be '...
                         ' calculated from ' Recorded '.']);
                  end 
              case 'I'
                  posZ=getPosGF('Z');
                  posSF=getPosGF('SF');
                  H=getH;
                  if (~isempty(H) && posSF>0)
                      struct_out.GeoMagV.I=acosd(H./D{posSF});
                  elseif (~isempty(H) && posZ>0)
                      struct_out.GeoMagV.I=atand(D{posZ}./H);
                  elseif (posSF>0 && posZ>0)
                      struct_out.GeoMagV.I=asind(D{posZ}./D{posSF});
                  else
                     error(['Vector Element ' element ' cannot be '...
                         ' calculated from ' Recorded '.']);
                  end 
              case 'F'
                  posX=getPosGF('X');
                  posY=getPosGF('Y');
                  posZ=getPosGF('Z');
                  H=getH;
                  posI=getPosGF('I');
                  if (posX>0 && posY>0 && posZ>0)
                      struct_out.GeoMagV.F=sqrt(D{posX}.^2+D{posY}.^2+...
                          D{posZ}.^2);
                  elseif (~isempty(H) && posZ>0)
                      struct_out.GeoMagV.F=sqrt(H.^2+D{posZ}.^2);
                  elseif (~isempty(H) && posI>0)
                      struct_out.GeoMagV.F=H./cosd(D{posI});
                      % For I=+/-90, F will be infinite, but should equal Z
                      % However, Z is not available
                      struct_out.GeoMagV.F(isinf(struct_out.GeoMagV.F))=NaN;
                  else
                     error(['Vector Element ' element ' cannot be '...
                         ' calculated from ' Recorded '.']);
                  end
              otherwise
                   error(['Requested vector element ' Requested(i) ...
              ' is not valid.']);
          end
      end
      
  end
  
  % Get time
  if ~(isempty(fieldnames(struct_out.GeoMagV)))
      pos_T=getPos('GeomagneticVectorTimes');
      struct_out.GeoMagV.Time=convertTime(D{pos_T},TimeFormat);
  end
  
  end
  
  %
  % GeoMagS
  %
  if (GeoMagS.set)
    posS=getPosGF('S');
    if (posS>0)
        struct_out.GeoMagS.S=D{posS};
        pos_T=getPos('GeomagneticScalarTimes');
        struct_out.GeoMagS.Time=converTime(D{pos_T},TimeFormat);
    else
        disp(['Warning: Scalar field is not available.']);
    end
  end
  
  %
  % Global Attributes
  %
  for j=1:length(GAttr.classes)
      
      name=GAttr.classes{j};
      
      if (GAttr.(name).set)
          for i=1:length(GAttr.(name).class)
              try
                  tmp=GAttr.(name).class{i};
                  struct_out.(name).(tmp)=M.GlobalAttributes.(tmp);
              catch
              end
          end
      end
      
  end
  
  %
  % User-specified attributes
  %
  if (UAttr.set)
      for j=1:length(UAttr.req)
          try
              tmp=UAttr.req{j};
              req=M.GlobalAttributes.(tmp);
              if (iscell(req) && length(req)==1)
                  struct_out.(tmp)=req{1};
              else
                  struct_out.(tmp)=req;
              end
          catch
              disp(['Warning: Attribute with name ' tmp ' not found.']);
          end
      end
  end
   
  %
  % Temperature
  %
  if (Temp.set)
      for j=1:length(Temp.req)
          tmp=num2str(Temp.req(j));
          posD=getPos(['Temperature' tmp]);
          posT=getPos(['Temperature' tmp 'Times']);
          if (posD>0 && posT>0)
              struct_out.(['Temperature' tmp]).Data=D{posD};
              FIELDNAM=M.VariableAttributes.FIELDNAM;
              struct_out.(['Temperature' tmp]).FIELDNAM=...
                  FIELDNAM{strcmp(M.Variables{posT,1},FIELDNAM(:,1)),2}
              struct_out.(['Temperature' tmp]).Times=...
                  convertTime(D{posT},TimeFormat);
          else
              disp(['Warning: Temperature ' tmp ' not found.']);
          end
      end
      
   end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                                                                       %
  %                ADDITIONAL INTERNAL FUNCTIONS                          %
  %                                                                       %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Convert time to specified epoch
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function Tc=convertTime(T,TimeFormat)
        switch TimeFormat.req
            case 'tt2000'
                Tc=T;
            case 'datenum'
                Tc=spdftt2000todatenum(T);
            case 'datetime'
                tmp_=spdftt2000todatenum(T);
                Tc=datetime(tmp_,'convertfrom','datenum');
            case 'jd2000'
                Tc=spdftt2000todatenum(T)-730486;
            case 'posix'
                Tc=(spdftt2000todatenum(T)-719529)*68400;
        end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Get horizontal field component from available GeoMagV variables
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function H=getH

  pos_H=getPosGF('H');
  
  % H is readily available
  if (pos_H>0)
      
      H=D{pos_H};
      
  % H has to be calculated
  else
      
      pos_X=getPosGF('X');
      pos_Y=getPosGF('Y');
      pos_I=getPosGF('I');
      pos_D=getPosGF('D');
      pos_SF=getPosGF('SF');
      if (pos_X>0 && pos_Y>0)
          H=sqrt(D{pos_X}.^2+D{pos_Y}.^2);
      elseif (pos_I>0 && pos_SF>0)
          H=cosd(D{pos_I}).*D{pos_F};
      elseif (pos_X>0 && pos_D>0)
          % For D=0, H will be infinite, but should equal Y
          % However, Y is not available
          H=D{pos_X}./sind(D{pos_D});
      elseif (pos_Y>0 && pos_D>0)
          % For D=+/-90, H will be infinite, but should equal X
          % However, X is not available
          H=D{pos_Y}./cosd(D{pos_D});
      else
          H=[];
      end
      
      % Set infinite numbers to NaN
      H(isinf(H))=NaN;
      
  end
  
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Get position of a variable corresponding to the specified element
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function pos=getPos(element)
        
        i_=find(strncmp(element,M.Variables(:,1),inf)==1,1,'first');
        if (i_>0)
            pos=i_;
        else
            pos=-1;
        end
    end

    function pos=getPosGF(element)
        
        if (strcmp(element,'SF'))
            pos_S=getPos('GeomagneticFieldS');
            pos=getPos('GeomagneticFieldF');
            if (pos<=0) pos=pos_S; end
        else
            pos=getPos(['GeomagneticField' element]);
        end
        
    end

    toc

end