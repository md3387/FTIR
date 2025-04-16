%FTIR External Detector Data conversion

%MDH 02/10/2019
%Converts a single voltage interferogram, mirror trigger, and HeNe Fringe
%into wavenumber vs transmission

%References
%1. Griffiths, P.R., deHaseth, J.A., "Fourier Transform Infrared Spectrometry 2nd Ed." Wiley 2007
%2. Thermo Fisher Scientific "FT-IR Spectrometers Nicolet iS50 Spectrometer User Guide" 2013

%Uses subfunctions:
%.....FindMirrorTrigger
%.....TDMS_readTDMSFile
%..........TDMS_handleGetDataOption
%..........TDMS_preprocessFile
%..........TDMS_getGroupChanNames
%..........TDMS_handleGetDataOption
%.....TDMS_readChannelOrGroup

%Ctrl+R to comment out block of code.
%Ctrl+T to uncomment

clear all
close all
clc
tic

m1 = msgbox('Please select the DATE, you would like to examine.','Warning','warn');
uiwait(m1);
location = uipickfiles('FilterSpec','C:\Users\mitchell.hageman\OneDrive - afacademy.af.edu\Desktop\Course History\AE471 - This Semester\Data','Output','cell')';
[pathstr, date] = fileparts(char(location));
run_no = str2num([char(inputdlg('Which run do you want to examine?: ','Run Number'))]);


%IF USING .TDMS FILES
run_file = [char(location),'\',date,'_FTIR (', num2str(run_no) ').tdms'];
filestruct = TDMS_readTDMSFile(run_file); %load file structure parameters of tdms file
datarows= filestruct.numberDataPointsRaw(6); %sets to read in any number of rows.  However, this assumes the .tdms file structures are identical.
%load interferogram, mirror trigger, and HeNeFringe data from tdms file
[InterferogramData,channelNames(1)] = TDMS_readChannelOrGroup(run_file,'Oscilloscope - Waveform Data','Channel 6 (PXI1Slot2/6)',[1 datarows]);
[MirrorTriggerData,channelNames(2)] = TDMS_readChannelOrGroup(run_file,'Oscilloscope - Waveform Data','Channel 5 (PXI1Slot2/5)',[1 datarows]);
[HeNeFringeData,channelNames(3)] = TDMS_readChannelOrGroup(run_file,'Oscilloscope - Waveform Data','Channel 4 (PXI1Slot2/4)',[1 datarows]);
data_full = [InterferogramData',MirrorTriggerData',HeNeFringeData'];

write_file = [char(location),'\Processed_data',date,'.xlsx'];
data_size = size(data_full);

%Plot full data record
figure(1)
plot (data_full(:,1),'DisplayName','Interferogram')
hold on
plot (data_full(:,2),'DisplayName','Mirror Trigger')
plot (data_full(:,3),'DisplayName','HeNe Fringe')
legend
xlabel('Record Length [row]')
ylabel('Voltage [V]')
xlim([0, length(data_full)])
set(gca, 'FontSize', 18) % Set axis tick labels size
hold off

%filter out unnecessary data, leaving only the portion containing the interferogram
[rowstart, rowend] = FindMirrorTrigger(data_full); %define beginning & end of interferogram
data_parsed=data_full(rowstart:rowend,:); %data parsed for only that corresponding to interferogram,
interferogramvoltage=data_parsed(1:length(data_parsed)-5,1); %column vector of detector voltages - 1st column of Oscope output
HeNeFringevoltage=data_parsed(1:length(data_parsed)-5,3); %column vector of laser output (HeNe Fringe) - 2nd column of Oscope output
x=1:length(HeNeFringevoltage); %row vector - row numbers corresponding to data
[HeNepks, a] = findpeaks(HeNeFringevoltage,'MinPeakHeight',0.055); %find positive LOCAL peaks of HeNe Fringe. Will include some false peaks(MinPeakheight removes negative local maxima, but closely spaced positive local maxima will remain.
%[HeNepks,a]= [peak values, indices]

%Approximate period with false peaks included
for j=1:length(a)-1
    period1(j)=a(j+1)-a(j);
end
MinPeakDist=0.95*mean(period1); %averages length of all periods
figure(2)
plot(period1,'DisplayName','HeNe period1')
hold on

%2nd iteration: Use approximated period to define Min. Peak distance, and perform a refined search of HeNepks.
[HeNepks, a] = findpeaks(HeNeFringevoltage,'MinPeakHeight',0.055,'MinPeakDistance',MinPeakDist);
%Now determine the average period of the sin wave with both negative and closely spaced local maxima removed
localperiod=zeros(1,length(HeNeFringevoltage));
localDistancePerx=zeros(1,length(HeNeFringevoltage));
for j=1:length(a)-1
    period2(j)=a(j+1)-a(j);
end
period=mean(period2); %[index difference]



HeNe_wavelength=632.8E-9; %[m] = 632.8nm Wavelength of a Helium-Neon Laser
%plots to double-check sorting
figure(2)
plot(period2,'DisplayName','HeNe period2')
legend
hold off

DistancePerx=HeNe_wavelength/period; %[m] mirror travel distance represented by each successive data point

%So, every [period = # data points], the mirror travels [DistancePerx] meters.
%FilteredInterferogramVoltage = lowpass(interferogramvoltage(:,i),1,2000000) %2Mhz

%Apodization
centerburst=find(interferogramvoltage==min(interferogramvoltage)); %Returns the index of the centerburst
OpticalPathDifference=(DistancePerx*x)-(DistancePerx*centerburst); %Sets optical path difference using centerburst to define zero path difference.

%Here I'm trying to create and use a "Local period" to determine OPD
%instead of an average period (Distanceperx).
count=1;
addOn=0;
peakindex=1;
DistancePerx(1)=DistancePerx;
CheckPeaks=zeros(1,length(x));
%localperiod(1:a(1))=a(1); 
for j=2:length(x)
    if j<a(length(a)-1) %true for all but the last (partial) periods.
        if j<=a(1)
            localperiod(j)=a(1);%takes care of first (partial) period j=1:a(1)
            localDistancePerx(j)=DistancePerx;
        elseif j<a(peakindex)
            localperiod(j)=a(peakindex)-a(peakindex-1);
            localDistancePerx(j)=HeNe_wavelength/localperiod(j);
        else
            localperiod(j)=localperiod(j-1);
            localDistancePerx(j)=HeNe_wavelength/localperiod(j);
            peakindex=peakindex+1;
        end
    else
        localperiod(j)=length(x)-a(length(a));%final (partial) period
    end
    
    count=count+addOn;
    DeltaOPD=OpticalPathDifference(j)-OpticalPathDifference(1); %Original Delta
    if DeltaOPD<count*HeNe_wavelength
        addOn=0;
        CheckPeaks(j)=0;
    else
        addOn=1;
        CheckPeaks(j)=0.057;
    end
end
figure(3)
plot(HeNeFringevoltage,'DisplayName','HeNe Fringe')
hold on
plot(x(a), HeNepks, '^g', 'MarkerFaceColor','g','DisplayName','Peaks')
grid
plot(x, CheckPeaks, '^b', 'MarkerFaceColor','b','DisplayName','CheckPeaks')
axis([0  500    ylim])
hold off

% %Recalculate using local distance per x
%  OpticalPathDifference=(localDistancePerx.*x)-(mean(localDistancePerx(1:centerburst))*centerburst);
% %re-plot
% count=1;
% addOn=0;
%  for j=2:length(x)
%     count=count+addOn;
%     DeltaOPD=OpticalPathDifference(j)-OpticalPathDifference(1); %Original Delta
%     if DeltaOPD<count*HeNe_wavelength
%         addOn=0;
%         CheckPeaks(j)=0;
%     else
%         addOn=1;
%         CheckPeaks(j)=0.057;
%     end
% end
% figure(4)
% plot(HeNeFringevoltage,'DisplayName','HeNe Fringe')
% hold on
% plot(x(a), HeNepks, '^g', 'MarkerFaceColor','g','DisplayName','Peaks')
% grid
% plot(x, CheckPeaks, '^b', 'MarkerFaceColor','b','DisplayName','CheckPeaks')
% axis([0  500    ylim])
% hold off


d=1-(OpticalPathDifference/max(OpticalPathDifference)).^2; %Ref 1 - Eq.2.23 - this is the bracketed quantity.
ApodizationFunction=(0.152442*(d.^0))-(0.136176*(d.^1))+(0.983734*(d.^2));% Medium Norton-Beer Apodization. Ref 1 - Eq.2.23
ApodizedInterferogramvoltage=ApodizationFunction'.*interferogramvoltage; %Ref 1
%Plot apodization
figure(5)
plot (interferogramvoltage,'DisplayName',' Interferogram')
hold on
plot (ApodizedInterferogramvoltage,'DisplayName','Apodized Interferogram')
set(gca, 'FontSize', 18) % Set axis tick labels size
legend
hold off

rawfft=fft(interferogramvoltage); %discreet fourier transform of detector voltage
%Phase Shift Correction
Re_fft=real(rawfft);% Real part of fft
Im_fft = imag(rawfft); %Imaginary part of fft.
Phase_Angle=atan(Im_fft./Re_fft);
%TrueSpectrum=Re_fft.*cos(Phase_Angle)+Im_fft.*sin(Phase_Angle);
%amp=abs(TrueSpectrum);
amp=2*abs(rawfft/length(x)); % Two-sided fft spectrum -column vector equal in length to "x."
ampreal=amp(1:0.5*length(rawfft)); %Single-sided fft spectrum - column vector half the length of "amp"

%MirrorTravel=DistancePerx*length(x); % [m] original calculation
MirrorTravel=HeNe_wavelength*length(period2)% [m] equivalent to original, but should be more precise.
DeltaWN=(1/MirrorTravel)/100; %[cm^-1] Correlate wavenumber change(deltaWN) to mirror position
WN=DeltaWN.*x'; %[cm^-1] column vector of WN values equal in length to x
WNreal=WN(1:0.5*length(rawfft)); %[cm^-1] column vector of only the first half of WN values - equal in length to "ampreal."
WLreal=10000000./WNreal; %[nm] column vector of wavelengths equal in length to "ampreal" and "WNreal" (nm=10,000,000/cm^-1)

%Plot Transmission vs WN
figure(7)
xllim=2000; %[cm-1] Bottom x-axis left limit
xrlim=8000; %[cm-1] Bottom x-axis right limit
xspace=(xrlim-xllim)/6; %spacing of both x axes
x1ticklabels=xllim:xspace:xrlim;
h2=axes; % Define a plot for the top (secondary) x-axis (waveLENGTH) first
% Wavelength and Wavenumber are inversely proportional (nonlinear relationship) so you have to define the top x-axis tics & spacing.
x2llim=10000000/xrlim; %[cm-1] Top x-axis 2 left limit
x2rlim=10000000/xllim; %[cm-1] Top x-axis 2 right limit
x2tickmagnitude=fliplr(10000000./x1ticklabels); % row vector of non-linear tick values
x2ticks=(x1ticklabels-xllim)/(xrlim-xllim); % define linear spacing of top x-axis (equal to spacing of main x-axis)
set(h2, 'Xdir', 'reverse') %cause x-axis to go from highest number (x2llim) to lowest number (x2rlim)
xticks([x2ticks]);
x2ticklabels=char(split(num2str(x2tickmagnitude)));
xticklabels(x2ticklabels);
%note that no data is being plotted yet.
%Only the top (secondary) x-axis is being established.
set(h2, 'XAxisLocation', 'Top') %move x-axis to top of graph
xlabel('wavelength [nm]')
set(h2, 'Color', 'None') %eliminate color so primary graph can be overlayed
set(h2, 'Ytick', []) %eliminat tick marks for secondary Y-axis, essentially eliminating it.

%now plot the data on the primary axes (waveNUMBER vs intensity).
h1=axes;
plot(WNreal, ampreal) %, 'Color', 'k')
xlabel('wavenumber [cm-1]')
ylabel('single beam power spectrum [Arb Units]')
set(h1, 'Xlim', [xllim, xrlim])
set(gca, 'FontSize', 18) % Set axis tick labels size
hold on
grid

% columnheaders=[{'interferogramvoltageAvg'},{'WLrealAvg'},{'WNrealAvg'},{'amprealAvg'}];
% xlswrite(write_file,columnheaders,1,'A1');
% xlswrite(write_file,interferogramvoltageAvg,1,'A2');
% xlswrite(write_file,WLrealAvg,1,'B2');
% xlswrite(write_file,WNrealAvg,1,'C2');
% xlswrite(write_file,amprealAvg,1,'D2');
toc