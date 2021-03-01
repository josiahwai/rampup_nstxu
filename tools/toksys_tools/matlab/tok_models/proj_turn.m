  function proj= proj_turn(conect,turns)

% proj_turn Generates projection matrix to combine series and anti-series coils.
%           optionally includes turns in output
%
% SYNTAX:   proj= proj_turn(conect,turns)
%
% INPUT:
%       conect=	connection vector containing id of coils to connect
%               same index coils are in series, negative are in anti-series,
%               0 excludes
%               ex: [1,2,2,3,2,0,5,-5]
%                    6th coil is not used, 7th & 8th are anti-series
%       turns=	turns in each coil element (or fraction of turns)
%               can be used to make a one turn system (same size as conect)
%               ex: [9,1/3,1/3,10,1/3,-3,8,8]
%                   coil 2 is a 1 turn system (1/3+1/3+1/3=1)
%               If not present assumes 1-turn in each coil
%
%  Examples:
%
%     Standard NSTX (no PF4, PF5u,l series:)
%             conect= [1,2,3,4,0,5,5,0,6,7,8,9];
%
%     NSTX (PF4u,l, PF5u,l series:)
%             conect= [1,2,3,4,5,6,6,5,7,8,9,10];
%
%     Convert to Neumeyer order: oh, 1au, 1al, ib 2u 2l 3u 3l 4 5:
%            conect= [1 2 5 7 9 10 10 9 8 6 3 4]; 
%               
% OUTPUT:
%       proj=		Projection matrix
%
% USE:
%       mccc=  proj*mcc*proj';
%       mcvv=  proj*mcv;
%       rescc= proj*resc;
%       names:
%         namess=char([]);
%         for ii=1:size(proj,1)
%           nam= char([]);
%           id= find(abs(proj(ii,:)));
%           for jj=1:length(id);
%             if proj(ii,id(jj)) < 0
%              nam=[nam remove_space(names(id(jj),:)) '_']; % '_' = anti-series
%             else    
%              nam=[nam remove_space(names(id(jj),:))];
%             end
%           end
%           namess= strvcat(namess,nam);
%         end    


% Jim Leuer 10-03
% =====================================================================
%  WRITTEN BY:  Jim Leuer 
%                    COPYRIGHT 2001 GENERAL ATOMICS 
%                       RESTRICTED RIGHTS NOTICE
% UNPUBLISHED-RIGHTS RESERVED UNDER THE COPYRIGHT LAWS OF THE UNITED STATES.
% =====================================================================
% input data check/defaults
  if nargin < 1
    disp('ERROR: proj_turn must have at least one argument');
    proj=[];
    return
  elseif nargin == 1
    turns= ones(size(conect));
  else
    if length(turns)~=length(conect)
      disp('ERROR: proj_turn length of turns and conect must be same');
      return
    end
  end;

% make sure inputs are row vectors
  [m,n]= size(conect);
  if m>n
    conect= conect';
  end
  [m,n]= size(turns);
  if m>n
    turns= turns';
  end
  

% ============================================================================
% Main algorithm
  idmin= max(1,min(abs(conect)));
  idmax= max(abs(conect));
  id= idmin:idmax;

% construct projection matrix
  projj= zeros(length(id),length(conect));
  for ii=1:length(conect)
    idd= conect(ii);
    if abs(idd) % if we put in a zero it will not not include in system
       projj(abs(idd),ii)= sign(idd);
    end
  end

% add turns
  if ~isempty(turns)
    tt= ones(length(id),1)*turns;
    projj= projj.*tt;
  end 

% eliminate row if all zeros (some ids from min to max are skipped)
  proj=[];
  count=0;
  for ii=1:length(id)
    if any(projj(ii,:))
      count= count+1;
      proj(count,:)= projj(ii,:);
    end
  end

  if count ~= length(id)
    disp(['%CAUTION: proj_turn eliminated all zero rows in PROJ; # eliminated= ',...
                 int2str(length(id)-count)]);
    disp(['          size of PROJ is: (',int2str(size(proj,1)),', ',...
                                     int2str(size(proj,2)),')'])
  end

  return

% Testing

  conect=[1,2,2,3,2,0,5,-5]        %6th coil is not used, 7th & 8th are anti-series
  turns= [9,1/3,1/3,10,1/3,-3,8,8] % coil 2 is a 1 turn system (1/3+1/3+1/3=1)
  conect=[1,2,2,3,2,4,-4]        %6th coil is not used, 6th & 7th are anti-series
  turns= [9,1/3,1/3,10,1/3,8,8] % 2nd coil set is a 1 turn system (1/3+1/3+1/3=1)
  conect=[1,2,2,3,2,4,+4];        %6th coil is not used, 6th & 7th are anti-series
  turns= [9,1/3,1/3,10,1/3,3,2]; % 2nd coil set is a 1 turn system (1/3+1/3+1/3=1)
  proj= proj_turn(conect,turns);
  a= eye(length(conect));
  b= proj*a*proj'
  c= ones(length(conect));
  d= proj*c*proj'

  conect= [ 2      1     0     -2];
  names=  ['aa1';'bb ';'cc ';'aa2']
  proj= proj_turn(conect);
  namess=char([]);
  for ii=1:size(proj,1)
    nam= char([]);
    id= find(abs(proj(ii,:)));
    for jj=1:length(id);
      if proj(ii,id(jj)) < 0
       nam=[nam remove_space(names(id(jj),:)) '_']; %put '_' for anti-series
      else    
       nam=[nam remove_space(names(id(jj),:))];
      end
    end
    namess= strvcat(namess,nam);
  end    
