function [c, ceq, GC, GCeq] = mySPLICE_upload(x,T1,T2,esp,ETL,target,sup,profileorder,hs)
% function for optimizing the flip angles in a SPLICE sequence. 
% Method is described in Rahbek et al: Optimized flip angle schemes for the
% diffusion-weighted SPLICE sequence.
%
% INPUT:
% x: the optimization variable (vector of flip angles)
% T1,T2,esp,ETL: tissue and sequence parameters relevant in the optim.
% target: The desired (normalized) k-space signal weighting
% sup: number of start-up-echoes
% profileorder: determining the k-space sampling order. Either 'linear' (default) or 'lowhigh' (same as center-out). 
% hs: the halfscan factor if halfscan, i.e. partial fourier imaging, is
% used. If left empty, halfscan is NOT used. 
% OBS! when using halfscan together with 'lowhigh', it is assumed that the 
% highest value of the target function (max(target)) applies to the center of k-space.  
% 
% Note: a single filter for the sum of the two echo-families is calculated.
% I.e. not taking oscillations in the individual families into account, as these in reality are not present to
% the extend they appaer when simulating with EPG. (Because of imperfect
% sliceprofile effect). 
% 
% By Sofie Rahbek, June 2022
% -------------------------------------------------------------------------
% default values, if no input: 
if nargin > 9
    hs = [];
    if nargin < 8
        profileorder = 'linear'; %default profile order is linear ;
        if nargin < 7
            sup=1; 
        end
    end
end

% Normalizing target (if user input is not normalized):
target=target./sum(target);

% reordering target according to sampling order. E.g. if center-put (low-high) sampling
% is used, the center of the target fct is sampled first. 
N = ETL-sup;
if strcmp(profileorder,'lowhigh') || strcmp(profileorder,'lh')
    order = zeros(numel(target),1);
    Nc = round((N+1)/2);
    order(1) = Nc;
    order(2:2:N) = [Nc-1:-1:1];
    order(3:2:N) = [Nc+1:N];
    
    
    if ~isempty(hs) && hs<1
        order = zeros(numel(target),1);
        [~,Nc] = max(target);
        order(1:2:Nc*2-1) = [Nc:-1:1];
        order(2:2:Nc*2-1)=[Nc+1:Nc*2-1];
        order(Nc*2:end)=[Nc*2:numel(target)];
    end
    
    target = target(order);

end


[s1,s2,~,~] = epg_splice(x,ETL,T1,T2,esp);

S1 = s1(sup+1:end);
S2 = s2(sup+1:end);
%

% Summing echo-signals:
filt = target./(abs(S1)+abs(S2)); % currently using sum (or average)
% (consider using SOS - but usually the 2 echo families are carriyng 
% approximately the same amount of signal, so would not make a differnece)

SNR = 1/sum(filt.^2); 

c = -SNR; % cost function
ceq = [];

%self-defined gradients: (could be implemented to improve computation time)
GC = []; 
GCeq = [];

end