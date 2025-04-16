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

m1 = msgbox('Please select the directory, BY DATE, corresponding to the data which you would like to examine.  SELECT ONLY ONE DAY.','Warning','warn');
uiwait(m1);
%location = uipickfiles('FilterSpec','C:\Users\FTIR\Desktop','Output','cell')';
location = uipickfiles('FilterSpec','C:\Users\mitchell.hageman\Desktop\Course History\499 - FTIR - Griffin\FTIR Computer','Output','cell')';
[pathstr, date] = fileparts(char(location));
number_of_runs = str2num([char(inputdlg('Please enter number of runs for this day: ','Number of Runs'))]);
choice=questdlg('Do you want to exclude any runs?','exclude_runs','Yes','No', 'no');
if strcmpi(choice,'Yes')
    skip_runs_txt = inputdlg('List the runs you want to exlude separated by spaces:','skip runs',[1 50]);
    skip_runs_num=str2num(skip_runs_txt{:});
else
    skip_runs_num=0;
end

for i = 1:number_of_runs
    if ismember(i,skip_runs_num)
        run_no = i
    else
        run_no = i
        %IF USING .CSV FILES
        %run_file = [char(location),'\',date,'_FTIR', num2str(run_no)
        %data_full = xlsread(run_file);
        
        %IF USING .TDMS FILES
        run_file = [char(location),'\',date,'_FTIR (', num2str(run_no) ').tdms'];
        filestruct = TDMS_readTDMSFile(run_file); %load file structure parameters of tdms file
        datarows= filestruct.numberDataPointsRaw(6); %sets to read in any number of rows.  However, this assumes the .tdms file structures are identical.
        %load interferogram, mirror trigger, and HeNeFringe data from tdms file
        [InterferogramData,channelNames(1)] = TDMS_readChannelOrGroup(run_file,'Oscilloscope - Waveform Data','Channel 1 (PXI1Slot2/1)',[1 datarows]);
        [MirrorTriggerData,channelNames(2)] = TDMS_readChannelOrGroup(run_file,'Oscilloscope - Waveform Data','Channel 3 (PXI1Slot2/3)',[1 datarows]);
        [HeNeFringeData,channelNames(3)] = TDMS_readChannelOrGroup(run_file,'Oscilloscope - Waveform Data','Channel 4 (PXI1Slot2/4)',[1 datarows]);
        data_full = [InterferogramData',MirrorTriggerData',HeNeFringeData'];
        
        write_file = [char(location),'\Processed_data',date,'.xlsx'];
        data_size = size(data_full);
        
        %Plot full data record
        %         figure(1)
        %         plot (data_full(:,1))
        %         hold on
        %         plot (data_full(:,2))
        %         plot (data_full(:,3))
        %         hold off
        
        %filter out unnecessary data, leaving only the portion containing the interferogram
        [rowstart, rowend] = FindMirrorTrigger(data_full); %define beginning & end of interferogram
        data_parsed=data_full([rowstart:rowend],:); %data parsed for only that corresponding to interferogram,
        HoldRecordLength(i)=length(data_parsed)-5;
        interferogramvoltage(:,i)=data_parsed(1:HoldRecordLength(1),1); %column vector of detector voltages - 1st column of Oscope output
        HeNeFringevoltage(:,i)=data_parsed(1:HoldRecordLength(1),3); %column vector of laser output (HeNe Fringe) - 2nd column of Oscope output
        x=1:length(HeNeFringevoltage(:,i)); %row vector - row numbers corresponding to data
        [HeNepks, a] = findpeaks(HeNeFringevoltage(:,i),'MinPeakHeight',0.06); %find positive LOCAL peaks of HeNe Fringe. Will include some false peaks(MinPeakheight removes negative local maxima, but closely spaced positive local maxima will remain. 
        %Approximate period with false peaks included
        for j=1:length(a)-1
            period1(j)=a(j+1)-a(j);
        end
        MinPeakDist=0.95*mean(period1); %averages length of all periods
        %2nd iteration: Use approximated period to define Min. Peak distance, and perform a refined search of HeNepks. 
        [HeNepks, a] = findpeaks(HeNeFringevoltage(:,i),'MinPeakHeight',0.06,'MinPeakDistance',MinPeakDist);
        %Now determine the average period of the sin wave with both negative and closely spaced local maxima removed
        for j=1:length(a)-1
            period2(j)=a(j+1)-a(j);
        end
        period=mean(period2);
       
        HeNe_frequency=4.47E14; %[Hz] or [1/s] Frequency of a Helium-Neon Laser
        HeNe_wavelength=632.8E-9; %[m] = 632.8nm Wavelength of a Helium-Neon Laser
        
        %plots to double-check sorting
        figure(2)
        plot (interferogramvoltage)
        figure(3)
        plot(period2)
        %figure(4)
        %plot(HeNeFringevoltage)
        hold on
        %plot(x(a), HeNepks, '^g', 'MarkerFaceColor','g')
        %hold off
        grid
        %axis([0  500    ylim])
        
        DistancePerx=HeNe_wavelength/period; %[m] mirror travel distance represented by each successive data point
        %FilteredInterferogramVoltage = lowpass(interferogramvoltage(:,i),1,2000000) %2Mhz
        
        %Apodization
        centerburst(i)=find(interferogramvoltage(:,i)==min(interferogramvoltage(:,i))); %Returns the index of the centerburst
        OpticalPathDifference(:,i)=(DistancePerx*x)-(DistancePerx*centerburst(i)); %Sets optical path difference using centerburts to define zero path difference. 
        d(:,i)=1-(OpticalPathDifference(:,i)/max(OpticalPathDifference(:,i))).^2; %Ref 1 - Eq.2.23 - this is the bracketed quantity.
        ApodizationFunction(:,i)=(0.152442*(d(:,i).^0))-(0.136176*(d(:,i).^1))+(0.983734*(d(:,i).^2));% Medium Norton-Beer Apodization. Ref 1 - Eq.2.23 
        ApodizedInterferogramvoltage(:,i)=ApodizationFunction(:,i).*interferogramvoltage(:,i); %Ref 1 
        %rawfft=fft(interferogramvoltage(:,i)); %discreet fourier transform of detector voltage
        rawfft=fft(ApodizedInterferogramvoltage(:,i)); %discreet fourier transform of detector voltage 
        
        %Phase Shift Correction
        Re_fft=real(rawfft);% Real part of fft
        Im_fft = imag(rawfft); %Imaginary part of fft.
        Phase_Angle=atan(Im_fft./Re_fft);
        %TrueSpectrum=Re_fft.*cos(Phase_Angle)+Im_fft.*sin(Phase_Angle);
        %amp=abs(TrueSpectrum);
        amp=2*abs(rawfft/length(x)); % Two-sided fft spectrum -column vector equal in length to "x." 
        ampreal(:,i)=amp(1:0.5*length(rawfft)); %Single-sided fft spectrum - column vector half the length of "amp"
        
        MirrorTravel=DistancePerx*length(x); % 
        DeltaWN=(1/MirrorTravel)/100; %[cm^-1] Correlate wavenumber change(deltaWN) to mirror position
        WN=DeltaWN.*x'; %[cm^-1] column vector of WN values equal in length to x
        WNreal(:,i)=WN(1:0.5*length(rawfft)); %[cm^-1] column vector of only the first half of WN values - equal in length to "ampreal."
        WLreal(:,i)=10000000./WNreal(:,i); %[nm] column vector of wavelengths equal in length to "ampreal" and "WNreal" (nm=10,000,000/cm^-1)
        %WNreal{i}(:)=WN(1:0.5*length(rawfft)); %[cm^-1] %If array changes
        %length every iteration, lots of variables will need to be changed
        %from "(:,i)" to "{i}(:)"
        
        %Plot individual Scans
        figure(7)
        xllim=2000; %[cm-1] Bottom x-axis left limit
        xrlim=8000; %[cm-1] Bottom x-axis right limit
        plot(WNreal(:,i), ampreal(:,i)) %, 'Color', 'k') !!!TEMP FIX. ASSUMES 1ST SET OF WN IS RIGHT!!!
        xlabel('wavenumber [cm-1]')
        ylabel('signal intensity [Arb Units]')
        hold on
        
    end
    
end
grid
hold off

interferogramvoltageAvg=mean(interferogramvoltage(:,size(interferogramvoltage,2)),2);  %Average interferogram
amprealAvg=mean(ampreal(:,size(ampreal,2)),2);  %Average amplitude
%WNrealAvg=mean(WNreal(:,size(WNreal,2)),2);  %Average WaveNUMBER -This
%will be wrong until you properly fix the offset in WNreal
WNrealAvg=WNreal(:,1);  %!!!TEMP FIX. ASSUMES 1ST SET OF WN IS RIGHT!!!
% WLrealAvg=mean(WLreal(:,size(WLreal,2)),2);  %Average WaveNUMBER - Same
% comment here as for WNrealAvg.
WLrealAvg=WLreal(:,1); %!!!TEMP FIX. ASSUMES 1ST SET OF WN IS RIGHT!!!

figure(8)
xllim=2000; %[cm-1] Bottom x-axis left limit
xrlim=8000; %[cm-1] Bottom x-axis right limit
xspace=(xrlim-xllim)/6; %spacing of both x axes
x1ticklabels=xllim:xspace:xrlim;
h2=axes; % Define a plot for the top (secondary) x-axis (waveLENGTH) first
% Wavelength and Wavenumber are inversely proportional (nonlinear relationship) so you have to define the top x-axis tics & spacing.
x2llim=10000000/xrlim; %[cm-1] Top x-axis 2 left limit
x2rlim=10000000/xllim; %[cm-1] Top x-axis 2 right limit
x2tickmagnitude=fliplr(10000000./x1ticklabels); % row vector of non-linearly tick values
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
plot(WNrealAvg, amprealAvg) %, 'Color', 'k')
xlabel('wavenumber [cm-1]')
ylabel('single beam power spectrum [Arb Units]')
set(h1, 'Xlim', [xllim, xrlim])
grid

columnheaders=[{'interferogramvoltageAvg'},{'WLrealAvg'},{'WNrealAvg'},{'amprealAvg'}];
xlswrite(write_file,columnheaders,1,'A1');
xlswrite(write_file,interferogramvoltageAvg,1,'A2');
xlswrite(write_file,WLrealAvg,1,'B2');
xlswrite(write_file,WNrealAvg,1,'C2');
xlswrite(write_file,amprealAvg,1,'D2');
toc

%xlswrite('C:\Users\mitchell.hageman\Desktop\ApodizationCheck.xlsx',amprealAvg,1,'F2');

%interfer=channelData'
%[channelData',channelNames] = TDMS_readChannelOrGroup('20190122_FTIR1.tdms','Oscilloscope - Waveform Data','Channel 1 (PXI1Slot2/1)',[1 100039])