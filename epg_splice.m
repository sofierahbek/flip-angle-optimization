
function [s1,s2,phasediag,P] = epg_splice(flipangle,etl,T1,T2,esp)

%	EPG Simulation of SPLICE readout sequence.  First flip angle
%	is 90 (excitation) (Here chosen about y axis)
%   The following flips are refocusing (can be about x or y axis - not important in SPLICE)
%
%   INPUTS:    
%	flipangle: vector of refocusing flip angles (radians)
%	etl: echo train length, if single flipangle is to be repeated.
%	T1,T2,esp: relaxation times and echo spacing (arb. units).
%   OUTPUTS:
%   s1 and s2: the two echo signals (families)
%   phasediag: stored data for the phase states to plot phase diagram. Currently plotting is uncommented. 
%   P: matrix of all states
%
%	Note that if refoc flip angle is constant, less than pi, and etl > 1 the
%	first refocusing flip angle is the average of 90 and the desired
%	refocusing flip angle, as proposed by Hennig.
%
%	All states are kept, for demo purposes, though this 
%	is not really necessary.

%   The function is a further development of an EPG simulation of an CPMG experiment from: 
%   Stanford University Rad229 Class Code: MRI Signals and Sequences. See http://web.stanford.edu/class/rad229   
% 
%   By Sofie Rahbek, June 2022
% -------------------------------------------------------------------------------------

if (nargin < 1) 
    flipangle = pi*ones(1,12); 
end
% -- Default parameters:  ETL = length(flipangle)
if (nargin < 2 || length(etl)==0) 
    etl = length(flipangle); 
end

if (length(flipangle)==1) && (etl > 1) && (abs(flipangle)<pi)	
  % -- 1st flip reduced trick (Hennig)
  flipangle(2)=flipangle(1);
  flipangle(1)=(pi+abs(flipangle(2)))/2/abs(flipangle(2)) * flipangle(2);
end

if (etl > length(flipangle)) 
    flipangle(end+1:etl) = flipangle(end); 
end


if (nargin < 3) 
    T1 = 1; 
end %  sec
if (nargin < 4) 
    T2 = 0.1; 
end % sec
if (nargin < 5) 
    esp = 0.01; 
end	% sec

% First exciation with no displacement gradient:
cons_angle = pi/2*zeros(length(flipangle),1); % pi/2 or 0. It dosent matter for SPLICE

P = zeros(3,2*etl);		% Allocate all known states, 2 per echo.
P(3,1)=1;			% Initial condition/equilibrium.

Pstore = zeros(4*etl,etl);	% Store F,F* states, with 2 gradients per echo
Zstore = zeros(2*etl,etl);	% Store Z states, with 2 gradients per echo

% -- 90 excitation

P = epg_rf(P,pi/2,pi/2);	% Do 90 tip.
s1 = zeros(1,etl);		% Allocate signal vector to store.
s2 = s1;
for ech=1:etl
    
  P = epg_grelax(P,T1,T2,esp/4,1,0,1,1);   % -- Relax. during gradient
  P = epg_rf(P,abs(flipangle(ech)),cons_angle(ech));   % -- Refoc. RF
  P = epg_grelax(P,T1,T2,esp/4,1,0,1,1);   % -- Relax. during gradient
% E1:
  s1(ech) = P(1,1);  	% Signal is F0 state. (first echo) 
  
% E2:  
  P = epg_grelax(P,T1,T2,esp/2,1,0,1,1);   % -- prologation of gradient
  s2(ech) = P(1,1);  % (second echo)
  
  Pstore(2*etl:4*etl-1,ech) = P(2,:).';	% Put in negative states
  Pstore(1:2*etl,ech) = flipud(P(1,:).');  % Put in positive, overwrite center.
  Zstore(:,ech) = P(3,:).';

end

set(0,'defaultAxesFontSize',14);        % Default font sizes
set(0, 'DefaultLineLineWidth', 2);      % Default line width


plotstate = cat(1,Pstore,Zstore);
% subplot(1,2,1);
% plot([1:etl]*esp,abs(s));
% lplot('Echo Time','Signal','CPMG: Signal vs Echo Time');
% 
% subplot(1,2,2);
% dispim(plotstate);
% xlabel('Echo'); ylabel('F(top) and Z(bottom) States');
% title('F and Z states vs time');
phasediag = plotstate;	


