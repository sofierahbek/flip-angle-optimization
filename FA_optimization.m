% Example of optimizing the flipangle scheme for a SPLCIE sequence 
clear, close all;

%% Tissue and seq parameters:
% Several tissues:
T1csf = 2000;
T2csf = 250;
T1gm = 1000;
T2gm = 100;
T1wm = 800;
T2wm = 90;
T1fat = 300;
T2fat = 85;
T1phantom = 370;
T2phantom = 360; 
T1brain = 900;
T2brain = 95;

% Choose a tissue: 
T1choice = T1brain;
T2choice = T2brain;

% acquisition settings:
esp = 4.9; % [ms]
ETL = 55;
startups = 3; % minimum 1
profileorder = 'linear';
hsfactor = 1; % halfscan (partial fourier) factor. Choose 1, if no half-scan. Set to 0.64 for my testing. 
FOV = 20; % Only to define target resolution (not that important)
N = round(ETL/hsfactor)+startups; 

% save ID: 
ID = 'BRAIN_LH_st3_PI';

%% Define target function - here a modified Hanning function is used.
Nx = N-startups;
beta = 1; % Controls only amplitude. Not important, as the targetfunction is normalized to sum to one. 
alpha = 1.5; % Controls shape (width) of the modified Hanning function 
dx = FOV/Nx;
dk = 1/FOV;
Ksamples = ([1:Nx]-(Nx+1)/2)*dk; 
target = beta/2 * (1+cos((2*pi*Ksamples*dx)/alpha)); % the hanning function

% Visualize target (MTF), and the resulting PSF: 
figure('position',[400,400,700,250])
subplot(1,2,1)
plot(target), ylim([0 1])
title('MTF')   
[PSF,PSFtarget,ax,axfine]=MTF2PSF(target,'linear',FOV);
subplot(1,2,2)
plot(axfine,abs(PSFtarget))
title('PSF')
sgtitle(['Alpha=' num2str(alpha)])

% Effectuating possible half-scan: 
if ETL < length(target)
    target = target(end-ETL+1:end);
    N = ETL+startups;
end
[~,center] = max(target);


figure,plot(target), ylim([0 1])
title('Target during acquisition')   

%% Optimization using fmincon:
% Settings:
% lower and upper bounds for the flip angles:
lb = zeros(1,N); 
ub = pi*ones(1,N)+eps;
% intialization: 
%x0 = [120 80 50 45 50 60 80 90*ones(1,N-7)]*pi/180;
x0 = [100*ones(1,N)]*pi/180;

[A,b,Aeq,beq]=deal([]);

opt = optimset('MaxIter',500,'MaxFunEvals',10^6,'TolFun',10^(-9),'AlwaysHonorConstraints','none'); %,'PlotFcns','optimplotfval'

[x1,fval,exitflag,output] = fmincon(@(x)mySPLICE(x,T1choice,T2choice,esp,N,target./sum(target),startups,profileorder,hsfactor),x0,A,b,Aeq,beq,lb,ub,[],opt);

% Alternative use multistart: 
% gs=MultiStart;
% problem = createOptimProblem('fmincon','x0',x0,...
%     'objective',@(x)mySPLICE7(x,T1choice,T2choice,esp,N,0.5*target,startups,profileorder,hsfactor),'lb',lb,'ub',ub,'options',opt);
% [xtest,fval] = run(gs,problem,30);

% OR lsqnonlin: 
%[xx,resnorm,residual,exitflag,output]  = lsqnonlin(@(x)mySPLICE8(x,T1choice,T2choice,esp,N,target,startups,profileorder,hsfactor),x0,lb,ub,opt); 

%% Look at result
flipdeg = x1*180/pi; 

% The flip angles:
figure('position',[100 100 1100 300])
subplot(1,3,1)
plot(flipdeg),
ylim([0 180])
xlim([0 N])
title('FAs'), xlabel('Echo #'), ylabel('Degrees')

% The MTF for both echo families (before filtering):
[s1,s2,phasediag,P] = epg_splice(x1,N,T1choice,T2choice,esp);
MTF = abs(s1)+abs(s2);
SNR_nofilt = abs(sum([s1(startups+1:end) s2(startups+1:end)]));
subplot(1,3,2)
plot(abs(s1)), hold on,
plot(abs(s2))
plot(MTF,'--')
plot([startups startups],[0 0.9],'k--')
legend('E1','E2',['Sum, SNR=' num2str(round(SNR_nofilt,3,'significant'))],'Startup period'), title('MTF (unfiltered)')
ylabel('Magnitude signal (a.u.)'), xlabel('Echo #')

% Remove start-ups: 
S1 = s1(startups+1:end);
S2 = s2(startups+1:end);

% Ordering of samples, if center-out (low-high) readout:
%(for mySPLICE7):
if strcmp(profileorder,'lowhigh') || strcmp(profileorder,'lh')
    order = zeros(numel(target),1);
    
    Nc = round((numel(target)+1)/2);
    order(1) = Nc;
    order(2:2:numel(target)) = [Nc-1:-1:1];
    order(3:2:numel(target)) = [Nc+1:numel(target)];
   
    
    if ~isempty(hsfactor) && hsfactor<1
        order = zeros(numel(target),1);
        [~,Nc] = max(target);
        
        order(1:2:Nc*2-1) = [Nc:-1:1];
        order(2:2:Nc*2-1)=[Nc+1:Nc*2-1];
        order(Nc*2:end)=[Nc*2:numel(target)];
    end
    
    target_order = target(order);
    [~,reverse_order] = sort(order); 
else
    target_order = target;
    reverse_order=[1:numel(target)]; 

end

% The filter: 
filt = abs(target_order)./(abs(S1)+abs(S2)); 
[~, Nc] = max(target);
filt_apply = filt(reverse_order); % the filter to apply to recorded data - ss the sampling order then is irrelevant.
filt_norm = filt_apply./filt_apply(Nc); % the normalized filter to apply to recorded data. 

subplot(1,3,3)
plot([1:numel(filt)],filt_norm), 
%xticks([1 center numel(filt)])
%xticks([1 round(numel(filt)/2) numel(filt)])
%xticklabels({'start','k_{center}','1/2\cdotk_{FOV}'})
ylim([0 max(filt_norm)*1.2])
SNR = 1/sqrt(sum(filt.^2));
title('Correction filter')
legend(['SNR\propto' num2str(round(SNR,3,'significant'))])

% check-up: just to check that target function is obtained after filtering
%(compare with target_order): 
MTFfilt=(abs(S1)+abs(S2)).*filt; % = target function.  
figure,plot(MTFfilt), ylim([0 1])

%% save result:
close all

% save all parameters: 
%save(['FAopt_' ID '.mat']) % saving all parameters used for optimization

% OR save filter and flips: 
%filtername = '..; 
% saving filter and relevant parameters:
%save(['FAfilt_' ID '.mat'],'filt','filt_apply','filt_norm','flipdeg','startups','hsfactor','esp','target','target_order') 

% write flip angles to txt-file: 
fileID = fopen('FILE_FOR_FLIPS.txt','w');
fprintf(fileID,'%f\r\n',flipdeg) % use '%d\n' if it should save/print integers. 
fclose(fileID)



