%Find Mirror Trigger
%M. Hageman 2/10/2019
%code accounts for the possibility that the mirror trigger may be high or
%low at the beginning of the data record.  However, it assumes that the
%first falling edge begins the interferogram.

function[rowstart, rowend] = FindMirrorTrigger(data_full)

trigger = data_full(:,2); %raw mirror trigger signal
dtrigger = gradient(trigger); % Take Derivative of the mirror trigger signal
%figure
%plot(dtrigger)
%dtriggerabs=abs(dtrigger);
%triggerup=0.5; %[volts]
%triggergradientup=0.3; %approx. value of gradient when mirror trigger switches high
%triggergradientdown=-0.3; %approx. value of gradient when mirror trigger switches low
triggerup=0.094; %[volts]
triggergradientup=0.04; %approx. value of gradient when mirror trigger switches high
triggergradientdown=-0.04; %approx. value of gradient when mirror trigger switches low

if trigger(1)<triggerup
    waitedge=1;
else
    waitedge=0;
end

fallingedges=find(dtrigger<triggergradientdown)
risingedges=find(dtrigger>triggergradientup);
rowstart=fallingedges(2);
rowend=risingedges(waitedge+2);
%rowstart=fallingedges(1);
%rowend=risingedges(1);

%code for plotting peaks. Useful for troubleshooting code, but commented now for speed
% [dtriggerpks, x] = findpeaks(dtriggerabs) %
% figure(1)
% plot(trigger)
% hold on
% plot(dtrigger)
% plot(x, dtriggerpks, '^g', 'MarkerFaceColor','g')
% hold off
% grid
% axis([0  1E5    ylim])




