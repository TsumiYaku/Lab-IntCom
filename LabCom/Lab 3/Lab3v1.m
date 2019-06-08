clear all
close all

%% Parametri
MPAM = 2; 
nbits = 1e4; % Numero di bit trasmessi
Mbps = 100; % Velocità di trasmissione in Mbps
SpS = 10; % Campioni per simbolo

BpS = log2(MPAM); % Numero di bits per simbolo
Rs = Mbps/BpS; % Velocità di trasmissione simboli in Baud
Ts = 1/Rs; % Tempo di trasmissione di un simbolo
fs = SpS*Rs; % Banda di simulazione (freq di campionamento)
sym2alpha = [-1;1]; % Tabella di conversione da simboli a coefficienti PAM
alpha2sym = [0;1]; % Tabella di conversione da coefficienti PAM a simboli
CarrierFreq = 2e3; % Frequenza carrier sinusoidale

%% Calcolo bits

% Matricole affiancate in vettore
matr = [2,1,7,2,0,4,2,3,5,2,0,6,2,2,6,6,1,9];

Bits = [];
for i=1:length(matr)
    Bits = [Bits, de2bi(matr(i),8)];
end

%% Generazione segnale elettrico
alphas = 2*Bits-1; % Conversione da bits ad ampiezze del sengnale elettrico
sig = zeros(nbits*SpS,1);

for ii=1:nbits
    for jj=1:SpS
        sig((ii-1)*SpS+jj) = alphas(ii);
    end
end
sig = sig';

%% Calcolo e riproduzione del segnale in uscita con carrier 

samples = length(sig); % Campioni del segnale
t = linspace(0, length(Bits)/Mbps, length(sig));

sig_out = real(sig.*cos(2*pi*CarrierFreq*t));
%sound(sig_out, CarrierFreq+fs);
%% Plot delle PSD (test)

f1 = linspace(-Bsim, Bsim,length(sig));
f2 = linspace(-Bsim-CarrierFreq, Bsim+CarrierFreq,length(sig_out));

figure(1), plot(f1,abs(fftshift(fft(sig))).^2);
figure(2), plot(f2,abs(fftshift(fft(sig_out))).^2);

