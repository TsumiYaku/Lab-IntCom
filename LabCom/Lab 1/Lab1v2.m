V = 1; % Ampiezza massima dei segnali

figure(1)

%% Analisi voce registrata

% Registrazione audio
fsr = 8000; % Frequenza di campionamento
Tr = 5; % Tempo di registrazione

recorder = audiorecorder(fsr, 16, 1);
recordblocking(recorder, Tr);
sigrec = getaudiodata(recorder); % Segnale registrato
sigrec = sigrec'; % Trasposto per compatibilità
samples = length(sigrec); % Numero di campioni

% Spettro di potenza
psd = abs(fft(sigrec)).^2;
fr = linspace(-fsr/2, fsr/2, samples); % Valori asse delle frequenze

subplot(2,2,1);
plot(fr, fftshift(psd));
title('Spettro di potenza (voce registrata)');
xlabel('f')
grid on

% Densità di probabilità
subplot(2,2,2)
histogram(sigrec, 30);
title('Densità di probabilità (voce registrata)');

%% Analisi file audio

% Apertura file
[sigfile, fs] = audioread('record.wav');
sigfile = sigfile'; % Trasposto per compatibilità
fsamples = length(sigfile); % Numero di campioni

% Spettro di potenza
psd = abs(fft(sigfile)).^2;
f = linspace(-fs/2, fs/2, fsamples); % Valori asse delle frequenze

subplot(2,2,3)
plot(f, fftshift(psd));
title('Spettro di potenza (file audio)');
xlabel('f')
grid on

% Densità di probabilità
subplot(2,2,4)
histogram(sigfile, 30)
title('Densità di probabilità (file audio)')

%% Calcolo delle SNR

nbit = [4 6 8]; % Bit di quantizzazione
pe_theory = logspace(-9,-1,1e3); % Valori di probabilità (SNR teorica)
pe_sig = logspace(-9,-1,9); % Valori di probabilità (SNR segnale)

figure(2)

for i = 1:length(nbit)
    
    % Definizione partizioni quantizzazione
    M = 2^nbit(i); % Numero intervalli di quantizzazione
    DV = 2*V/M; % Passo di quantizzazione
    partition = -V+DV:DV:V-DV; % Partizione asse delle ampiezze
    codebook = -V+DV/2:DV:V-DV/2; % Valori quantizzati
    
    % SNR teorica
    SNRt = M^2./(1+4*(M^2-1)*pe_theory);
    
    % SNR registrazione vocale
    % Quantizzazione
    [index, quants] = quantiz(sigrec,partition,codebook);
    
    % Calcolo SNR
    SNR = fullSNR(sigrec, index, codebook, pe_sig); % SNR segnale
    
    subplot(3,2,2*(i-1)+1)
    line = ['SNR (voce registrata, ', num2str(nbit(i)), ' bit)'];
    semilogx(pe_theory, 10*log10(SNRt));
    hold on
    grid on
    semilogx(pe_sig, SNR, 'o');
    title(line);
    legend('SNR teorica', 'SNR segnale');
    xlabel('P_e')
    
    % Analisi file audio
    % Quantizzazione
    [index, quants] = quantiz(sigfile,partition,codebook);
    
    % Calcolo SNR
    SNR = fullSNR(sigfile, index, codebook, pe_sig); % SNR segnale
    
    subplot(3,2,2*(i-1)+2)
    line = ['SNR (file audio, ', num2str(nbit(i)), ' bit)'];
    semilogx(pe_theory, 10*log10(SNRt));
    hold on
    grid on
    semilogx(pe_sig, SNR, 'o');
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