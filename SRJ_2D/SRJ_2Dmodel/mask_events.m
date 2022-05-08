function [events,masked_events] = mask_events( events,mask )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
pad=1e3;
mask(:,2)=mask(:,2)-pad;
mask(:,3)=mask(:,3)+pad;

masked_events=zeros(length(events),1);
for c1=1:length(mask),
    masked_events=masked_events|(events(:,1)==mask(c1,1)&events(:,2)>=mask(c1,2)&events(:,2)<=mask(c1,3))|(events(:,4)==mask(c1,1)&events(:,5)>=mask(c1,2)&events(:,5)<=mask(c1,3));
end
events(masked_events,:)=[];

disp(['masked ', num2str(sum(masked_events)),' events']);
