clear all
close all

%% Definizione segnale tramite rand

V = 1; % Ampiezza massima del segnale (da -V a +V)
fc = 10e5; % Frequenza di campionamento (scelta arbitraria)
Tc = 1/fc; % Tempo di campionamento
samples = 1e6; % Numero di campioni del segnale (arbitrario)
t = (0:1:samples-1)*Tc; % Valori asse dei tempi dei campioni del segnale

sig = rand(1,samples); % Generazione del segnale
sig = (sig-0.5)*2*V; % Normalizzazione del segnale (media nulla e varianza pari a V)

%% Calcolo e plot spettro di potenza

psd = abs(fft(sig)).^2;
f = linspace(-fc/2, +fc/2, samples); % Valori trasformata asse delle frequenze

figure(1)

subplot(1,2,1)
grid on
plot(f, fftshift(psd))
xlabel('f')
ylabel('Sx')
title('Spettro di potenza')

%% Plot densità di probabilità

subplot(1,2,2)
histogram(sig, 30)
title('Densità di probabilità')

%% Calcolo delle SNR

nbit = [4 6 8];
peA = logspace(-9,-1, 1e3);
peB = logspace(-9, -1, 9);

figure (2)

for i = 1:length(nbit)
    
    % Definizione partizioni quantizzazione
    M = 2^nbit(i); % Numero intervalli di quantizzazione
    DV = 2*V/M; % Passo di quantizzazione
    partition = -V+DV:DV:V-DV; % Partizione asse delle ampiezze
    codebook = -V+DV/2:DV:V-DV/2; % Valori quantizzati
    
    % Quantizzazione
    [index, quants] = quantiz(sig,partition,codebook);
    
    % Calcolo SNR
    SNRt = M^2./(1+4*(M^2-1)*peA); % SNR teorica
    SNR = fullSNR(sig, index, codebook, peB); % SNR segnale
    
    subplot(3,1,i)
    line = ['SNR (', num2str(nbit(i)), ' bit)'];
    semilogx(peA, 10*log10(SNRt));
    hold on
    grid on
    semilogx(peB, SNR, 'o');
    title(line);
    legend('SNR teorica', 'SNR segnale');
    xlabel('P_e')
    
end

%% Funzioni di utilità

function SNR=fullSNR(sig, index, codebook, pe) % funzione per il calcolo di SNR su tutte le probabilità

SNR = zeros(1, length(pe));
indata = de2bi(index); % Codifica

for i = 1:length(pe)
    outdata = bsc(indata, pe(i)); % Simulazione trasmissione
    outidx = bi2de(outdata);
    vout = codebook(outidx+1);
    e = sig - vout; % Segnale d'errore/rumore
    SNR(i) = snr(sig, e);
end

end
