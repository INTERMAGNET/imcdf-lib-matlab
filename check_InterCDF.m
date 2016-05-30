%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITTEN BY ACHIM MORSCHHAUSER, GFZ POTSDAM, 2016
% mailto: mors/gfz-potsdam.de
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function checks the compatibility to the INTERMAGNET CDF data
% format. The checks are performed on a struct M containing the
% METADATA and a cell array D containing the data as read from a CDF file.
%
% The checks are in agreement with the format specified in the INTERMAGNET
% Discussion Document DD22, version 2.3, 23/01/2013.%
%
% IN:
% ===
%
% M     STRUCT      Contains the METADATA as read from an INTERMAGNET CDF
%                   file
% D     CELL ARRAY  Contains the data as read from an INTERMAGNET CDF file
%
% OUT:
% ====
%
% valid INTEGER     Describes the validity of M and D
%                   1:    The file completely conforms to the INTERMAGNET
%                         standard
%                   0:    The file does not conform to the INTERMAGNET
%                         standard, but still can be read
%                  -1/-2: The file does not conform to the INTERMAGNET
%                         standard and cannot be further processed as
%                         important data are missing
%
% desc  STRING      A description of the incompatibilities.
%                   Use fprintf(desc) to show a formatted list.
%
% EXAMPLE USAGE:
% ==============
% 
% cdf_file='test.cdf';
% M=spdfcdfinfo(cdf_file);
% D=spdfcdfread(cdf_file,'KeepEpochAsIs',1)
% [valid, desc]=check_InterCDF(M,D);
%
% CHECKS PERFORMED:
% =================
%
% 1) Check if global attributes conform to standard as specified in 
%    sections 4.1 - 4.6 of the INTERMAGNET Discussion Document 22,
%    ver. 2.3 (23/01/2013). 
%    In particular, the presence of the attributem, the correct data type,
%    and the allowed values are checked.
%
%    Some additional checks on attributes:
%    - |Latitude| < 90
%    - 0 <= Longitude <= 360
%
% 2) Check if all mandatory variables are present
%    - At least one geomagnetic field variable
%    - Exactly the elements specified in ElementsRecorded
%    - Temperature for INTERMAGNET_1-Second data
%    - Time Variables
%
% 3) Check if variable attributes are present and in correct format
%    (Section ImagCDF variable attributes on p. 8f.)
%
%    Some additional checks on attributes:
%    - Check if variables specified in DEPEND_0 are available
%    - Vector data msut be present for INTERMAGNET standards
%
% 4) Some consistency checks of the data
%    - double data type for geomagnetic data and temperature data
%    - tt2000 data type for timestamping data
%    - Data range consistent with FILLVAL, VALIDMIN, VALIDMAX
%    - Data range in limits of degrees of arc and degrees of celcius
%    - Data and Times have same length
%    - For INTERMAGNET 1s and 1m data, check that time stamps are
%      consistent
%
% 5) Filename is not checked
%     
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [valid,desc]=check_InterCDF(M,D)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Initialization
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  valid=1;
  desc='';

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Some required values
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  req_FormatDesc      = 'INTERMAGNET CDF Format';
  req_FormatVersion   = {'1.1'};
  req_Title           = 'Geomagnetic time series data';
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Some valid values definitions
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  valid_ElementsRecorded = {'X' 'Y' 'Z' 'H' 'D' 'E' 'V' 'I' 'F' 'S' 'G'};
  valid_PublicationLevel = {'1' '2' '3' '4'};
  valid_StandardLevel    = {'None' 'Partial' 'Full'};
  valid_Source           = {'INTERMAGNET' 'WDC'};
  valid_StandardName     = {'INTERMAGNET_1-Second' ...
                            'INTERMAGNET_1-Minute' ...
                            'INTERMAGNET_1-Minute_QD'};                    
  valid_UNITS_GN         = {'nT'};
  valid_UNITS_GD         = {'Degrees of arc'};
  valid_UNITS_T          = {'Celsius'};
                        
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Process Metadata
  %  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Variable Properties as defined in CDF
  VarProp=M.Variables;
  % Global Attributes
  gloAttr=M.GlobalAttributes;
  % Variable Attributes
  VarAttr=M.VariableAttributes;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Check if global attributes are in correct format
  %  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %
  % Check attributes describing the data
  %        (Sec. 4.1 of DD22)
  %
  
  check_entry(0,1,gloAttr,'FormatDescription','char',[1 1],req_FormatDesc);
  check_entry(0,1,gloAttr,'FormatVersion','char',[3 3],req_FormatVersion);
  check_entry(0,1,gloAttr,'Title','char',[1 1],req_Title);
  
  %
  % Check attributes that uniquely identify the data
  %        (Sec. 4.2 and Sec. 4.6 of DD22)
  %
  
  check_entry(0,1,gloAttr,'IagaCode','char',1);
  
  % ElementsRecorded has special formatting and needs a separate function
  ElementsRecorded=check_entry(0,1,gloAttr,'ElementsRecorded','char',1);
  check_ElementsRecorded(ElementsRecorded,valid_ElementsRecorded,[1 3]);
  
  check_entry(0,1,gloAttr,'PublicationLevel','char',[1 1],...
      valid_PublicationLevel);
  
  check_entry(0,1,gloAttr,'PublicationDate','char',1);
  
  %
  % Check attributes that describe the observatory
  %               (Sec. 4.3 of DD22)
  %
  
  check_entry(0,1,gloAttr,'ObservatoryName','char',1);
  Latitude=check_entry(0,1,gloAttr,'Latitude','double',1);
  if ~isempty(Latitude) && (abs(Latitude)>90)
      err_msg(1,'Latitude is not in range.'); 
  end
  Longitude=check_entry(0,1,gloAttr,'Longitude','double',1);
  if ~isempty(Longitude) && (Longitude>360 || Longitude <0) 
      err_msg(1,'Longitude is not in range.'); 
  end
  check_entry(0,1,gloAttr,'Elevation','double',1);
  
  % The Institution is added as valid source of the data
  Institution=check_entry(0,1,gloAttr,'Institution','char',1);
  if (~isempty(Institution)) valid_Source{end+1}=Institution; end
  
  % VectorSensOrientation has same format as ElementsRecorded
  VectorSensOrient=check_entry(0,0,gloAttr,'VectorSensOrient','char',1);
  check_ElementsRecorded(VectorSensOrient,valid_ElementsRecorded,[1 1]);
  
  %
  % Check attribute that relate to data standards and quality
  %              (Sec. 4.4 and Sec. 4.7 of DD22)
  %
  
  StandardLevel=check_entry(0,1,gloAttr,'StandardLevel','char',[1 1],...
      valid_StandardLevel);
  
  % StandardName is mandatory only if StandardLevel is set to 'Partial' or
  % 'Full' (Se. 4.7)
  ismand=(strcmp(StandardLevel,'Partial') || strcmp(StandardLevel,'Full'));
  StandardName=check_entry(0,ismand,gloAttr,'StandardName','char',[1 1],...
      valid_StandardName);
  
  check_entry(0,0,gloAttr,'StandardVersion','char',1);
  
  % PartialStandDesc is mandatory only if StandardLevel is Partial
  ismand=strcmp(StandardLevel,'Partial');
  PartialStandDesc=check_entry(0,ismand,gloAttr,...
      'PartialStandDesc','char',1);
  % PartialStandDesc has special formatting (Sec. 4.7)
  check_PartialStandDesc(PartialStandDesc,StandardName,1);
  
  %
  % Attributes that relate to publication of the data
  %                (Sec. 4.5 of DD22)
  %
  
  check_entry(0,1,gloAttr,'Source','char',[2 2],valid_Source);
  check_entry(0,0,gloAttr,'TermsOfUse','char',1,0);
  check_entry(0,0,gloAttr,'UniqueIdentifier','char',1);
  check_entry(1,0,gloAttr,'ParentIdentifiers','char',1);
  check_entry(1,0,gloAttr,'ReferenceLinks','char',1);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Create subsets of variable attributes for easier processing
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %
  % Find positions of specific variables in the set of recorded
  % variables. This will be needed for easier processing.
  %
  
  % GeomagneticField variables
  posG=find(strncmp('GeomagneticField',VarProp(:,1),16))';
  
  % GeomagneticField Vector variables with units of nT
  tmp=regexp(VarProp(:,1),'GeomagneticField[XYZHEVF]');
  posV_N=[];
  for i=1:length(tmp)
      if (~isempty(tmp{i}))
        posV_N=[posV_N i];
      end
  end
  
  % GeomagneticField Vector variables with units of degrees
  tmp=regexp(VarProp(:,1),'GeomagneticField[DI]');
  posV_D=[];
  for i=1:length(tmp)
      if (~isempty(tmp{i}))
        posV_D=[posV_D i];
      end
  end
  
  % GeomagneticField Scalar variables
  tmp=regexp(VarProp(:,1),'GeomagneticField[SG]');
  posS=[];
  for i=1:length(tmp)
      if (~isempty(tmp{i}))
        posS=[posS i];
      end
  end
  
  % Temperature variables
  tmp=regexp(VarProp(:,1),'Temperature[1-9]+$','start');
  posT=[];
  for i=1:length(tmp)
      if (~isempty(tmp{i}))
        posT=[posT i];
      end
  end
  
  % Variables with required variable attributes
  posA=[posG posT];
  
  % Vector geomagnetic variables
  posV=[posV_D posV_N];
  
  %
  % Select subsets of variable attributes for easier further processing
  %
  
  fields=fieldnames(VarAttr);
  for i=1:length(fields)
      % Variables that require attributes
      AVarAttr.(fields{i})=VarAttr.(fields{i})([posG posT],:);
      
      % Vector geomagnetic field
      VVarAttr.(fields{i})=VarAttr.(fields{i})(posV,:);
      % Vector geomagnetic field (nT)
      VNVarAttr.(fields{i})=VarAttr.(fields{i})(posV_N,:);
      % Vector geomagnetic field (Degrees)
      VDVarAttr.(fields{i})=VarAttr.(fields{i})(posV_D,:);
      % Scalar geomagnetic field
      SVarAttr.(fields{i})=VarAttr.(fields{i})(posS,:);
      % Geomagnetic (Scalar+Vector)
      GVarAttr.(fields{i})=VarAttr.(fields{i})(posG,:);
      % Temperatures
      TVarAttr.(fields{i})=VarAttr.(fields{i})(posT,:);
  end
  
  %
  % Time variables
  %
  
  V_Times=find(strcmp('GeomagneticVectorTimes',VarProp(:,1)));
  S_Times=find(strcmp('GeomagneticScalarTimes',VarProp(:,1)));
  
  T_Times=[];
  for i=1:length(posT)
      tmp=find(strcmp([VarProp{posT(i),1} 'Times'],VarProp(:,1)));
      if ~isempty(tmp)
          T_Times=[T_Times tmp];
      end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Check if all necessary variables are present
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %
  % Check if geomagnetic data are present
  %
  
  if isempty(posG)
      err_msg(3,'No geomagnetic data are present');
  end
  
  %
  % Vector data must be present if INTERMAGNET standard is specified
  %
  
  if (strncmp(StandardName,'INTERMAGNET_',12) && isempty(posV))
      err_msg(3,['Geomagnetic vector data must be present for '...
          'INTERMAGNET standard ' StandardName]);
  end
  
  
  %
  % Check that provided GeomagmeticField variables match the components
  % listed in the ElementsRecorded field of the METADATA.
  %
  
  if length(VarProp(posG,1))==length(ElementsRecorded)
      for i=1:length(ElementsRecorded)
          ElementsPresent(i)=VarProp{posG(i),1}(17);
      end
      if ~strcmp(sort(ElementsPresent),sort(ElementsRecorded))
        err_msg(3,['ElementsRecorded and present elements are' ...
            ' not consistent']);
      end
  else
     err_msg(3,['Number of ElementsRecorded and number of present' ...
            ' elements differ']);
  end
  
  %
  % At least one temperature is necessary for 1s data!
  %
  if (isfield(gloAttr,'StandardLevel') && ...
          ~strcmp(gloAttr.StandardLevel,'None') && ...
          strcmp(gloAttr.StandardName,'INTERMAGNET_1-Second') && ...
          isempty(posT))
      err_msg(3,'Temperature is mandatary for INTERMAGNET 1s-data');
  end
  
  %
  % Check that time variables are present (correct length is checked below)
  % 
  
  % GeomagneticVectorTimes
  if ~isempty(posV) && isempty(V_Times)
      err_msg(3,'GeomagneticVectorTimes are missing');
  % CDF only allows one instance of variable with same name (this check
  % will never be reached)
  %elseif length(tmp)>1
  %    err_msg(3,'GeomagneticVectorTimes multiply specified');
  end
  
  % Temperature
  for i=posT
      if ~any(strcmp([VarProp{i,1} 'Times'],VarProp(:,1)))
          err_msg(3,[VarProp{i,1} 'Times are missing']);
      end
  end
  
  % GeomagneticScalarTimes
  %if (~isempty(regexp(ElementsRecorded,'(S|G)','start')) && ...
  %tmp=strcmp('GeomagneticScalarTimes',VarProp(:,1));
  if ~isempty(posS) && isempty(S_Times)
      err_msg(3,'GeomagneticScalarTimes are missing');
  %elseif (length(tmp)>1)
  %    err_msg(3,'GeomagneticScalarTimes multiply specified');
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Check if all necessary variables attributes are present and correctly
  % formatted.
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  %
  % Check that every geomagnetic field or temperature variable has a
  % corresponding variable attribute set.
  % Also, create the variable posVF, which describes how Variables are 
  % mapped to VariableAttributes, i.e. the data at position posA(i) has
  % its attributes at position posVF(i).
  %
  
  posVF=[];
  tmp1=fieldnames(VarAttr);
  tmp1=tmp1{1};
  for i=posA
    tmp=strcmp(VarProp{i,1},VarAttr.(tmp1)(:,1));
    if sum(tmp)~=1
        err_msg(3,['VariableAttributes for ' VarProp{i,1} ' are missing']);
    else
        posVF=[posVF find(tmp)];
    end
  end
  
  %
  % Check variable attributes
  %
    
  % Check if FIELDNAM correctly set
  [~, valid_F]=check_entry(0,1,AVarAttr,'FIELDNAM','char',[1 1]);
  if (valid_F)
      for i=posG
          if (~strcmp(strrep(VarAttr.FIELDNAM{i,1},...
                  'GeomagneticField','Geomagnetic Field Element '),...
                  VarAttr.FIELDNAM{i,2}))
              err_msg(3,['Incorrect FIELDNAM for ' VarAttr.FIELDNAM{i,1}]);
          end
      end
      for i=posT
          if(~strncmp(VarAttr.FIELDNAM{i,1},VarAttr.FIELDNAM{i,2},11))
              err_msg(3,['Incorrect FIELDNAM for ' VarAttr.FIELDNAM{i,1}]);
          end
      end
  end
      
%       for i=1:length(FIELDNAM)
%           % Geomagnetic
%           if (strcmp(VarAttr.FIELDNAM{i,1}(1),'G'))
%               if(~strcmp(strrep(VarAttr.FIELDNAM{i,1},...
%                       'GeomagneticField','Geomagnetic Field Element '),...
%                       VarAttr.FIELDNAM{i,2}))
%                   err_msg(3,['Incorrect FIELDNAM for ' ...
%                       VarAttr.FIELDNAM{i,1}]);
%               end
%               % Temperature
%           else
%               if(~strncmp(VarAttr.FIELDNAM{i,1},VarAttr.FIELDNAM{i,2},11))
%                   err_msg(3,['Incorrect FIELDNAM for ' ...
%                       VarAttr.FIELDNAM{i,1}]);
%               end
%           end
%       end
%   end
  
  % Check Units
  if ~isempty(posT)
    check_entry(0,1,TVarAttr,'UNITS','char',[1 1],valid_UNITS_T);
  end
  check_entry(0,1,SVarAttr,'UNITS','char',[1 1],valid_UNITS_GN);
  check_entry(0,1,VNVarAttr,'UNITS','char',[1 1],valid_UNITS_GN);
  check_entry(0,1,VDVarAttr,'UNITS','char',[1 1],valid_UNITS_GD);
  
  % Data ranges and FILLVAL
  FILLVAL  = check_entry(0,1,AVarAttr,'FILLVAL','double',1);
  VALIDMIN = check_entry(0,1,AVarAttr,'VALIDMIN','double',1);
  VALIDMAX = check_entry(0,1,AVarAttr,'VALIDMAX','double',1);
   
  % Check dependencies(DEPEND_0)
  if ~isempty(posV)
    check_entry(0,1,VVarAttr,'DEPEND_0','char',[1 1],...
        'GeomagneticVectorTimes');
  end
  if ~isempty(posS)
    check_entry(0,1,SVarAttr,'DEPEND_0','char',[1 1],...
        'GeomagneticScalarTimes');
  end
  
  % Check if DEPEND_0 variables are available
  DEPEND_0=check_entry(0,0,VarAttr,'DEPEND_0','char',1);
  for i=DEPEND_0'
      if (~any(strcmp(i{1},VarProp(:,1))))
          err_msg(3,['Dependent time variable ' i{1} ' not found']);
      end
  end

  % Check display type
  check_entry(0,1,AVarAttr,'DISPLAY_TYPE','char',1,'time_series');
  
  % Check for consistent element code in LABLAXIS
  %TODO Implement check for temperature
  LABLAXIS=check_entry(0,1,AVarAttr,'LABLAXIS','char',1);
  for i=1:length(LABLAXIS)
      if (strcmp(VarProp{posA(i),1}(1),'G') && ...
              VarProp{posA(i),1}(17)~=LABLAXIS{posVF(i)})
          err_msg(3,['LABLAXIS for ' VarProp{posA(i),1} ' must be '...
              LABLAXIS{posVF(i)}]);
      end
  end
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Check consistency of data
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %
  % Check if data is in correct format
  %
  
  % Geomagnetic data and temperatures
  datavalid=1;
  for i=[posA posT]
      if ~isa(D{i},'double') || ~strcmp(VarProp{i,4},'double')
           err_msg(3,['Data not of type double for variable '...
               VarProp{i,1}]);
           datavalid=0;
      end
  end
  
  % Times
  for i=[T_Times V_Times S_Times]
      if ~isa(D{i},'int64') || ~strcmp(VarProp{i,4},'tt2000')
           err_msg(3,['Data not of type tt2000 for variable '...
               VarProp{i,1}]);
      end
  end
  
  %
  % Check if data is in correct range
  %
  
  if (datavalid)
  for i=1:length(FILLVAL)
      minD=min(D{posA(i)});
      maxD=max(D{posA(i)});
      if ( FILLVAL{posVF(i)}>minD && ...
              FILLVAL{posVF(i)}<maxD )
          err_msg(3,['Data range exceeds FILLVAL for ' ...
              VarProp{posA(i),1}]);
      end
      if ( VALIDMIN{posVF(i)}>minD)
          err_msg(3,['Data range exceeds VALIDMIN for ' ...
              VarProp{posA(i),1}]);
      end
      if ( VALIDMAX{posVF(i)}<minD)
          err_msg(3,['Data range exceeds VALIDMAX for ' ...
              VarProp{posA(i),1}]);
      end
  end
  end
  
  if any(abs(cell2mat(VALIDMAX(posV_D)))>360)
      err_msg(1,['VALIDMAX should be lower than |360| for '...
          VarProp{posV_D,1}]);
  end
  if any(cell2mat(VALIDMIN(posT)))<-273.15
      err_msg(1,['VALIDMIN should be larger than -273.15 for '...
          'temperatures']);
  end
  
  %
  % Check if data has correct length
  %
  
  % Geomagnetic vector times
  GV_times=find(strcmp('GeomagneticVectorTimes',VarProp(:,1)),1,'first');
  if (~isempty(GV_times))
      nod_V=length(D{GV_times});
      for i=posV
          if (length(D{i})~=nod_V)
              err_msg(3,['Number of data inconsistent with time stamps for '...
                  VarProp{i,1}]);
          end
      end
  end
  
  % Geomagnetic Scalar times
  tmp=find(strcmp('GeomagneticScalarTimes',VarProp(:,1)),1,'first');
  if (~isempty(tmp))
      nod_S=length(D{tmp});
      for i=posS
          if (length(D{i})~=nod_S)
              err_msg(3,['Number of data inconsistent with time stamps for '...
                  VarProp{i,1}]);
          end
      end
  end
  
  % Temperatures
  if (~isempty(posT))
      for i=posT
          tmp=find(strcmp([VarProp{i,1} 'Times'],VarProp(:,1)),1,'first');
          if ~isempty(tmp) && (length(D{i})~=length(D{tmp}))
              err_msg(3,['Number of data inconsistent with time stamps for '...
                  VarProp{i,1}]);
          end
      end
  end
  
  %
  % Check the starting time stamp of the data.
  % Check if sampling rate is the same as specified.
  %
  
  if ~isempty(GV_times)
    starttime=spdfbreakdowntt2000(D{GV_times}(1));
    
    if (strcmp(StandardName,'INTERMAGNET_1-Second'))
        % Milli, mikro, nano should be zero
        if sum(sum(starttime(7:9)))~=0
            err_msg(2,'Sampling not at full seconds');
        end
        % Sampling on 1s intervals
        if isempty(find(diff(D{4})~=1e9))
            err_msg(2,'Sampling not at 1s intervals');
        end
        
    elseif (strncmp(StandardName,'INTERMAGNET_1-Minute',20))
        % Seconds, Milli, mikro, nano should be zero
        if sum(sum(starttime(6:9)))~=0
            err_msg(2,'Sampling not at full seconds');
        end
        % Sampling on 1s intervals
        if isempty(find(diff(D{4})~=1e9*60))
            err_msg(2,'Sampling not at 1m intervals');
        end
    end
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Check for correct file name
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % TODO
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                                                                       %
  %                   Additional internal functions                       %
  %                                                                       %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Add error message to desc and set valid level.
  %
  %
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function err_msg(type,msg)
              
        if (type==1)
            desc=[desc 'Warning: ' msg '.\n'];
            valid_=0;
        elseif (type==2)
            desc=[desc 'Error: ' msg '.\n'];
            valid_=-1;
        elseif (type==3)
            desc=[desc 'Fatal Error: ' msg '.\n'];
            valid_=-2;
        end
        
        valid=min(valid,valid_);
        
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Check if CDF attribute entries are correctly formatted. The following
  % checks are performed:
  % - Presence of attribute         msg_type(1), see below
  % - Number of allowed entries     msg_type(1)
  % - Type of entry                 msg_type(1)
  % - Validity of entry             msg_type(2)
  %
  %
  % IN:
  % ===
  %
  % multiple    LOGICAL  Are multiple entries allowed for this attribute?
  % mandatory   LOGICAL  Is this attribute mandatory?
  % attr        STRUCT   Struct with given attributes as obtained from
  %                        CDF file (e.g. M.VariableAttributes)
  % name        STRING   Name of the attribute to check
  % name_type   STRING   Class/Type of the attribute, e.g. double
  % msg_type    DOUBLE   Type of error message as follows:
  %                        - 1: Warning message
  %                        - 2: Error message
  %                        - 3: Severe error message
  % varargin    CELL     Strings with allowed attribute
  %                      values
  %             STRING   Fixed and required attribute name
  %
  % OUT:
  % ====
  %
  % V           CELL     Attribute entries (if multiple==1)
  %             VAR      Attribute entry as respective type
  %                      (if multiple==0)
  % valid       LOGICAL  Was the check valid (1)
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [V, valid_]=check_entry(multiple,mandatory,attr,name,name_type,...
            msg_type,varargin)
        
        V=[];
        valid_=1;
        
        % check if attribute is present
        if (~isfield(attr,name))
            
            if (mandatory)
                err_msg(msg_type(1),[name ' is missing']);
                valid_=0;
            end
            
        else
            
            % Number of entries
            nvar=size(attr.(name),1);
            isvar=size(attr.(name),2)-1;
            V=cell(1,nvar);
            
            for i_=1:nvar
                
                % Variable attributes
                if (isvar)
                    V{i_}=attr.(name)(i_,2);
                % Global attributes
                else
                    V{i_}=attr.(name);
                end
                
                % check number of entries
                if (~multiple && length(V{i_})>1)
                    err_msg(msg_type(1),['Only one entry allowed for ' name]);
                    valid_=0;
                end
                
                
                for j=1:length(V{i_})
                    
                    % check type of entry
                    if (~isa(V{i_}{j},name_type))
                        
                        err_msg(msg_type(1),[name ' is of wrong type']);
                        valid_=0;
                        
                        % if valid_/required argument provided
                    elseif (nargin==7)                       
                        
                        % Check if attribute name is a valid_name
                        if (iscell(varargin{1}) && ...
                                ~any(strcmp(V{i_}{j},varargin{1})))
                            
                            err_msg(msg_type(2),...
                                ['Specified ' name ' not supported']);
                            valid_=0;
                            
                        % Check if the required attribute name is given
                        elseif (ischar(varargin{1}) &&  ...
                                (~strcmp(V{i_}{j},varargin{1})))
                            err_msg(msg_type(2),...
                                ['Wrong entry for ' name ': '...
                                V{i_}{j} ', but one of (' varargin{1} ...
                                ') is required']);
                            valid_=0;
                            
                        end
                        
                    end
                    
                end
                
            end
            
            if (~multiple)
                if (nvar > 1)
                    Vo=V;
                    V=cell(length(Vo),1);
                    for i_=1:length(Vo)
                        V{i_}=Vo{i_}{1};
                    end
                elseif (nvar==1)
                    V=V{1}{1};
                else
                    V=NaN;
                end
            end
            
        end
        
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check validity of ElementsRecorded.
%
% IN:
% ===
% values        STRING       Recorded values, e.g. 'HDZ'
% valid_values  CELL         Valid values
% msg_type      DOUBLE       Type of error message as follows:
%                            - 1: Warning message
%                            - 2: Error message
%                            - 3: Severe error message
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function check_ElementsRecorded(values,valid_values,msg_type)
        
        for i_=1:length(values)
            % Check if valid date elements were specified
            if ~any(strcmp(values(i_),valid_values))
                err_msg(msg_type(1),['Non-standard data records: ' values(i_)])
            end
            %  % Check if variables as specified in metadata are present in
            %  % timeseries data is now done later in the code
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check validity of PartialStandDesc.
%
% IN:
% ===
%
% ParStandDesc     STRING     PartialStandDesc as read from CDF
% StandardName     STRING     StandardName as read from CDF
% msg_type         DOUBLE     Type of error message as follows:
%                             - 1: Warning message
%                             - 2: Error message
%                             - 3: Severe error message
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function check_PartialStandDesc(ParStandDesc,StandardName,msg_type)
        
        pos=strfind(ParStandDesc,'IMO');
        
        if ~isempty(pos)
            if (all(diff(pos)~=8))
                err_msg(msg_type,'PartialStandDesc not correctly formatted');
            elseif ( strcmp(StandardName,'INTERMAGNET_1-Minute') )
                if (isempty(strfind(ParStandDesc(pos+3),'S')))
                    err_msg(msg_type,'Only IMOM is valid');
                else
                    for i_=pos+5
                        if (str2double(ParStandDesc(i_):ParStandDesc(i_+1))>23)
                            err_msg(msg_type,...
                                'Invalid number for IMOM in PartialStandDesc');
                        end
                    end
                end
            elseif ( strcmp(StandardName,'INTERMAGNET_1-Second') )
                if (isempty(strfind(ParStandDesc(pos+3),'M')))
                    err_msg(msg_type,'Only IMOS is valid');
                else
                    for i_=pos+5
                        if (str2double(ParStandDesc(i_):ParStandDesc(i_+1))>42)
                            err_msg(msg_type,...
                                'Invalid number for IMOS in PartialStandDesc');
                        end
                    end
                end
            end
            
        end
        
    end
 
end