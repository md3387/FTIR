%Integrated Absorbance Solver
%Uses subfunction "Beads" -
%https://www.mathworks.com/matlabcentral/fileexchange/49974-beads-baseline-estimation-and-denoising-with-sparsity44

%%%%!!!!NOTE!!!!%%%%%
%You built the SpectralSimulationDatabase.mat and LiquidFuelDatabase.mat
%differently.  You'll have to comment/uncomment lines of code if you want
%to switch between them.  At some point, decide which form is better,
%convert both to be the same, and ask ChatGPT how to replace the commands
%LiquidFuelData.(fieldname).wavenumber;  AND
%WaveNumber=SpectralSimulationDatabase(RunNumber).wavenumber;  
%With something generic


%%User Inputs
%Database='SpectralSimulationDatabase.mat';
Database='LiquidFuelDatabase.mat';
RunNumber=25;
WaveNumberStart=2800;%2985  %3000 used for CH4
WaveNumberEnd=3000;%3020  %3020 used for CH4
RunBaselineCorrection='n'; %n= no, don't run it.  y= yes, run it.
%%Load Data
load(Database); %Load Data
fuelNames = fieldnames(LiquidFuelKey);
%WaveNumber=SpectralSimulationDatabase(RunNumber).wavenumber;
%Absorbance=SpectralSimulationDatabase(RunNumber).absorbance;
fieldname = sprintf('LiquidFuel%d', RunNumber);
WaveNumber = LiquidFuelData.(fieldname).wavenumber;
Absorbance = LiquidFuelData.(fieldname).absorbance;

%Print Database metadata
whos %print structure parameters
%SpectralSimulationDatabase(RunNumber).Metadata %print metadata
%Database(RunNumber).Metadata %print metadata
LiquidFuelData.(fieldname).Metadata;


[~, WaveNumberStartIndex] = min(abs(WaveNumber - WaveNumberStart));
[~, WaveNumberEndIndex] = min(abs(WaveNumber - WaveNumberEnd));
%WaveNumberStartIndex=find(WaveNumber==WaveNumberStart);
%WaveNumberEndIndex=find(WaveNumber==WaveNumberEnd);

%BaselineStart=Absorbance(WaveNumberStartIndex);
%BaselineEnd=Absorbance(WaveNumberEndIndex);
Absorbance=Absorbance(WaveNumberStartIndex:WaveNumberEndIndex);
WaveNumber=WaveNumber(WaveNumberStartIndex:WaveNumberEndIndex);
%Baseline=zeros(length(WaveNumber),1);

if RunBaselineCorrection=='n'
    Baseline=zeros(length(WaveNumber),1);
else %if RunBaselineCorrection=='y'
    %% Run the BEADS algorithm

% Filter parameters
fc = 0.00001;     % fc : cut-off frequency (cycles/sample)
d = 1;          % d : filter order parameter (d = 1 or 2)
% Positivity bias (peaks are positive)
r = 6;          % r : asymmetry parameter
% Regularization parameters
amp = 0.8;      
lam0 = 0.5*amp;
lam1 = 5*amp;
lam2 = 4*amp;

tic
%Beads Inputs:
%   d: Filter order (d = 1 or 2)
%   fc: Filter cut-off frequency (cycles/sample) (0 < fc < 0.5)
%   r: Asymmetry ratio
%   lam0, lam1, lam2: Regularization parameters
%Beads Outputs:
%   cost: Cost function history
%   Baseline: Estimated baseline
%   x: Estimated sparse-derivative signal
[x1, Baseline, cost] = beads(Absorbance, d, fc, r, lam0, lam1, lam2); 
toc
%% Display the output of BEADS
ylim1 = [min(Absorbance) max(Absorbance)];
xlim1 = [min(WaveNumber) max(WaveNumber)];

figure(1)
clf

subplot(4, 1, 1)
plot(WaveNumber,Absorbance)
title('Data')
xlim(xlim1)
ylim(ylim1)
set(gca,'ytick', ylim1)

subplot(4, 1, 2)
plot(WaveNumber,Absorbance,'color', [1 1 1]*0.7)
hold on
%plot(WaveNumber,Baseline,'LineWidth', 1)
line(1:length(Absorbance), Baseline, 'LineWidth', 1)
legend('Data', 'Baseline')
legend boxoff
title(['Baseline, as estimated by BEADS', ' (r = ', num2str(r), ', fc = ', num2str(fc), ', d = ', num2str(d),')'])
xlim(xlim1)
ylim(ylim1)

set(gca,'ytick', ylim1)

subplot(4, 1, 3)
plot(WaveNumber,x1)
title('Baseline-corrected data')
xlim(xlim1)
ylim(ylim1)
set(gca,'ytick', ylim1)

subplot(4, 1, 4)
plot(WaveNumber,Absorbance - x1 - Baseline)
title('Residual')
xlim(xlim1)
ylim(ylim1)
set(gca,'ytick', ylim1)

orient tall
print -dpdf example

%% Display cost function history

figure(2)
clf
plot(cost)
xlabel('iteration number')
ylabel('Cost function value')
title('Cost function history')
end

%pre-allocate matrices
dW=zeros(length(Absorbance)-1,1);
dA=zeros(length(Absorbance)-1,1);

for i=1:(length(Absorbance)-1)
    dW(i)=(WaveNumber(i+1)-WaveNumber(i));
    dA(i)=(Absorbance(i)-Baseline(i))*dW(i);
end 
A=sum(dA);

figure
plot(WaveNumber,Absorbance,"DisplayName","Absorbance")
hold on
plot(WaveNumber,Baseline,"DisplayName","Baseline")
%plot(WaveNumber(1:length(dA)),dA,"DisplayName","dA")
set(gca, 'FontSize', 18);
legend('FontSize', 18)
xlabel('Wavenumber [cm-1]','FontSize', 18)
ylabel('Absorbance [A.U.]','FontSize', 18)
%ylim([0 4.5])
hold off




function [x, f, cost] = beads(y, d, fc, r, lam0, lam1, lam2)

% [x, f, cost] = beads(y, d, fc, r, lam0, lam1, lam2)
%
% Baseline estimation and denoising using sparsity (BEADS)
%
% INPUT
%   y: Noisy observation
%   d: Filter order (d = 1 or 2)
%   fc: Filter cut-off frequency (cycles/sample) (0 < fc < 0.5)
%   r: Asymmetry ratio
%   lam0, lam1, lam2: Regularization parameters
%
% OUTPUT
%   x: Estimated sparse-derivative signal
%   f: Estimated baseline
%   cost: Cost function history

% Reference:
% Chromatogram baseline estimation and denoising using sparsity (BEADS)
% Xiaoran Ning, Ivan W. Selesnick, Laurent Duval
% Chemometrics and Intelligent Laboratory Systems (2014)
% doi: 10.1016/j.chemolab.2014.09.014
% Available online 30 September 2014

% The following parameter may be altered.
Nit = 30;       % Nit: Number of iterations
pen = 'L1_v2';  % pen : penalty function for sparse derivative ('L1_v1' or 'L1_v2')
EPS0 = 1e-6;    % cost smoothing parameter for x (small positive value)
EPS1 = 1e-6;    % cost smoothing parameter for derivatives (small positive value)

switch pen
    case 'L1_v1'
        phi = @(x) sqrt(abs(x).^2 + EPS1);
        wfun = @(x) 1./(sqrt(abs(x).^2 + EPS1));
    case 'L1_v2'
        phi = @(x) abs(x) - EPS1 * log(abs(x) + EPS1);
        wfun = @(x) 1./( abs(x) + EPS1);
    otherwise
        disp('penalty must be L1_v1, L1_v2')
        x = []; cost = []; f = [];
        return
end

theta = @(x) sum(x(x>EPS0)) - r * sum(x(x<-EPS0)) ...
    + sum( (1+r)/(4*EPS0)*x(abs(x)<=EPS0).^2 ...
    + (1-r)/2 * x(abs(x)<=EPS0) + EPS0*(1+r)/4 );

y = y(:);
x = y;
cost = zeros(1, Nit);
N = length(y);
[A, B] = BAfilt(d, fc, N);
H = @(x) B*(A\x);
e = ones(N-1, 1);
D1 = spdiags([-e e], [0 1], N-1, N);
D2 = spdiags([e -2*e e], 0:2, N-2, N);
D = [D1;  D2];
BTB = B'*B;

w = [lam1 * ones(N-1, 1); lam2 * ones(N-2, 1)];
b = (1-r)/2 * ones(N, 1);
d = BTB * (A\y) - lam0 * A' * b;

gamma = ones(N, 1);

for i = 1:Nit
    
    Lambda = spdiags(w.*wfun(D*x), 0, 2*N-3, 2*N-3);
    
    k = abs(x) > EPS0;
    gamma(~k) = ((1 + r)/4) / abs(EPS0);
    gamma(k) = ((1 + r)/4) ./  abs(x(k));
    Gamma = spdiags(gamma, 0, N, N);
    
    M = 2 * lam0 * Gamma + D' * Lambda * D;
    x = A * ((BTB + A'*M*A)\d);
    
    cost(i) = 0.5 * sum(abs(H(y - x)).^2) + lam0 * theta(x) ...
        + lam1 * sum(phi(diff(x))) + lam2 * sum(phi(diff(x, 2)));
end

f = y - x - H(y-x);

end


% --- local function ----

function [A, B] = BAfilt(d, fc, N)
% [A, B] = BAfilt(d, fc, N)
%
% Banded matrices for zero-phase high-pass filter.
% The matrices are 'sparse' data type in MATLAB.
%
% INPUT
%   d  : degree of filter is 2d (use d = 1 or 2)
%   fc : cut-off frequency (normalized frequency, 0 < fc < 0.5)
%   N  : length of signal

b1 = [1 -1];
for i = 1:d-1
    b1 = conv(b1, [-1 2 -1]);
end
b = conv(b1, [-1 1]);

omc = 2*pi*fc;
t = ((1-cos(omc))/(1+cos(omc)))^d;

a = 1;
for i = 1:d
    a = conv(a,[1 2 1]);
end
a = b + t*a;

A = spdiags( a(ones(N, 1), :), -d:d, N, N);   % A: Symmetric banded matrix
B = spdiags(b(ones(N, 1), :), -d:d, N, N);    % B: banded matrix

end