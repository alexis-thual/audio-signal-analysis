%Script for multi-pitch detection (reference: Klapuri)
% but with pitch detection with spectral product 
% 
% Author: G. Richard, Janv. 2005 - MAJ:2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;


%Lecture du signal
[x,Fs]=audioread('A3_piano.wav');
soundsc(x,Fs);

spect_smooth=2 % options: 0: no spectral smoothness
               %          1: spectral smoothness principle with triangle 
               %           2 spectral smoothness principle with mean of three harmonics 

N=floor(0.7*Fs);      % Window size of analysed signal (only one window of signal is analysed)
Fmin=100;             % Minimal F0 frequency that can be detected
Fmax=900;             % Maximal F0 frequency that can be detected
H=4;                  % H = nombre de versions compress�es
prod_spec = 1;        %  m�thode pour la d�tection de pitch est produit spectral sinon par autocorrelation
freq_fond=[];         % tableau contenant la valeur des fr�quences fondamentales
seuil_F0 = 0.005;  %seuil pour l'acceptation d'un nouveau pitch

%Minimal frequency resolution
dF_min=Fs/N;             

%beta coef d'inharmonicit�
beta=0;

%alpha: coefficient for harmonic location around the true theoric harmony
alpha=0.02;

%Window
w=hamming(N);

%Beginning of signal (e.g. attack) is discarded
offset=floor(0.1*Fs);
xw=x(offset+1:offset+N).*w;    %xw est la fenetre de signal analys�

%Minimal number of data points to satisfy the minimal frequency resolution
Nfft_min=Fs/dF_min;

%compute the smallest power of two that satisfies the minimal frequency resolution for FFT
p = nextpow2(Nfft_min);
Nfft=2^p;

%calcul FFT
Xk=fft(xw,Nfft);

%frequency resolution of FFT
df=Fs/Nfft;

%normalisation
Xk=Xk/max(abs(Xk)+eps);

%"Reduced" frequency
f=[0:Nfft-1]/Nfft+eps;

% fr�quence maximale
Rmax = floor((Nfft-1)/(2*H));



        figure(1);
        plot([1:Nfft/2]*(Fs/Nfft), 20*log10(abs(Xk(1:Nfft/2))),'b');
        axis([0 5000 -80 0]);
        xlabel('Fr�quences (Hz)','fontsize',14);
        ylabel('Spectre d''Amplitude (dB)','fontsize',14);      
        title('Spectre d''amplitude du signal original','fontsize',14);
        pause

       
%loop on the number of pitches
%example criterion could use an energy ratio
while criterion > seuil_F0
    
    %detection of main F0
        %Compute spectral product
        % locate maximum
        %store value of estimated F0
criterion =seuil_F0+1
    %Subtraction of main note (Main F0 with its harmonics)
        
        %localisation of harmonics around theoretical values (with or without inharmonicy coeeficient) 
        % beta: harmonicity coefficient ;  alpha: coefficient of tolerance
                
        % Harmonic suppression (wideness of an harmonic to be suppressed
        % depends on the main lob of the TF of the analysis window)
            % suppression of harmonics is done on abs(Xk) on forcing all values
            % of a harmonic peak to the minimum value of the peak (e.g. the
            % level of noise).
         
end
%end of loop