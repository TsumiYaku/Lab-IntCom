clear all

nbit = [4 6 8];
peA = logspace(-9,-1, 1e3);
peB = logspace(-9, -1, 9);

% Definizione segnale tramite rand

V = 1; % Ampiezza massima del segnale (da -V a +V)
fc = 10e5; % Frequenza di campionamento (scelta arbitraria)
Tc = 1/fc; % Tempo di campionamento
samples = 1e6; % Numero di campioni del segnale (arbitrario)
t = (0:1:samples-1)*Tc; % Valori asse dei tempi dei campioni del segnale
sig = rand(1,samples); % Generazione del segnale
sig = (sig-0.5)*2*V; % Normalizzazione del segnale (media nulla e varianza pari a V)

% Quantizzazione e analisi del segnale

for i = 1:length(nbit)
    
    % Definizione partizioni quantizzazione
    
    M = 2^nbit(i); % Numero intervalli di quantizzazione
    DV = 2*V/M; % Passo di quantizzazione
    partition = -V+DV:DV:V-DV; % Partizione asse delle ampiezze
    codebook = -V+DV/2:DV:V-DV/2; % Valori quantizzati
    
    % Quantizzazione
    [index, quants] = quantiz(sig,partition,codebook);
    
    figure(1)
    
    % Spettro di potenza
    psd = abs(fft(quants)).^2;
    f = linspace(-fc/2, +fc/2, samples); % Valori trasformata asse delle frequenze
    
    subplot(3,2,2*(i-1)+1)
    grid on
    plot(f, fftshift(psd))
    xlabel('f')
    ylabel('Sx')
    line = ['Spettro di potenza (', num2str(nbit(i)),' bit)'];
    title(line)
    
    % Densità di probabilità
    subplot(3,2,2*(i-1)+2)
    histogram(quants, M)
    line = ['Densità di probabilità (', num2str(nbit(i)), ' bit)'];
    title(line)
    
    % Calcolo SNR
    SNR = zeros(1,length(peB));
    
    SNRt = M^2./(1+4*(M^2-1)*peA);
    words = de2bi(index, nbit(i));
    outquants = zeros(1, length(words));    
    for j = 1:length(peB)
        outdata = bsc(words, peB(j));
        outidx = bi2de(outdata);
        outidx = outidx';
        vout = codebook(outidx + 1);
        e = sig - vout;
        SNR(j) = snr(sig, e);
    end
    
    figure(2)
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

