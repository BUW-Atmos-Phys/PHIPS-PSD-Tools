function [treshist,V,V_err] = generate_histogram(PhipsData, dp, PhipsData_level_0, taxis,tstep,bins,sensitive_area,airspeed,time) 
% Function for Generating histograms for PHIPS SD

%% Make threshist and calculate volume
for i=1:length(taxis)-1
    t1 = taxis(i)-datenum(0,0,0,0,0,1/1000); %1ms tolerance, because sometimes the timestamps are off by 1/100 ms
    t2 = taxis(i+1)+datenum(0,0,0,0,0,1/1000);
    
    %% calculate SD
    % Select time
    % including ForceTriggers for deadtime calculation
    idx= PhipsData_level_0.RealTimeStamp>=t1 & PhipsData_level_0.RealTimeStamp<t2;
    treshist_level_0(i,:)=sum(idx); %histcounts(data,binlimits) sorts data into the number of bins determined by binlimits (intensity)
    % without ForceTriggers for SD
  
    idx= PhipsData.RealTimeStamp>=t1 & PhipsData.RealTimeStamp<t2;
    treshist(i,:)=histcounts(dp(idx),bins); %histcounts(data,binlimits) sorts data into the number of bins determined by binlimits (intensity)
    % Calculate probed volume
    idx = time>=t1 & time<t2;
    qdead = (1-sum(treshist_level_0(i,:),2)./tstep*12*10^-6);% Fraction of missed volumne due to el. dead time of 12*10^-6s
    V(i,:) = sensitive_area.* nanmedian(airspeed(idx)).*100.*tstep*qdead; % The sampling Volume [cm-3]
    V_err(i,:) = tstep.*qdead.*sqrt((0.275.*sensitive_area.*nanmedian(airspeed(idx)).*100).^2 + (sensitive_area.*0.01.*nanmedian(airspeed(idx)).*100).^2); %an error of 27.5% for the Area and 1% for the airspeed.  
    
end


end



