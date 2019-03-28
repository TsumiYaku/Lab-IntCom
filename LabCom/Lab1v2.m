nbit = [4 6 8]; % Bit di quantizzazione
cases = length(nbit);
V = 1; % Ampiezza massima dei segnali
peA = logspace(-9,-1,1e3);
peB = logspace(-9,-1,9);
SNR = zeros(1,length(peB));

% Registrazione audio
fsr = 8000; % Frequenza di campionamento
Tr = 3; % Tempo di registrazione

recorder = audiorecorder(fsr, 16, 1);
recordblocking(recorder, Tr);
sigrec = getaudiodata(recorder); % Segnale registrato

rsamples = length(sigrec); % Numero di campioni
fr = linspace(-fsr/2, fsr/2, rsamples); % Valori asse delle frequenze
    

% Apertura file audio
[sigfile, fsf] = audioread('record.wav');
fsamples = length(sigfile); % Numero di campioni
Tf = fsamples * 1/fsf; % Durata del segnale
ff = linspace(-fsf/2, fsf/2, fsamples); % Valori asse delle frequenze

for i = 1:length(nbit)
    
    % Definizione partizioni quantizzazione
    M = 2^nbit(i); % Numero intervalli di quantizzazione
    DV = 2*V/M; % Passo di quantizzazione
    partition = -V+DV:DV:V-DV; % Partizione asse delle ampiezze
    codebook = -V+DV/2:DV:V-DV/2; % Valori quantizzati
    SNRt = M^2./(1+4*(M^2-1)*peA);
    
    % Analisi registrazione vocale
    % Quantizzazione
    [index, quants] = quantiz(sigrec,partition,codebook);
    
    figure(1)
    
    % Spettro di potenza
    psd = abs(fft(quants)).^2;
    
    subplot(cases,2,2*(i-1)+1);
    line = ['Spettro di potenza (voce registrata ', num2str(nbit(i)), ' bit)'];
    title(line);
    grid on
    hold on
    plot(fr, fftshift(psd));
    xlabel('f')
    
    % Densità di probabilità
    subplot(cases,2,2*(i-1)+2)
    histogram(quants, M)
    line = ['Densità di probabilità (voce registrata ', num2str(nbit(i)), ' bit)'];
    title(line)
    hold on
    
    % Calcolo SNR
    words = de2bi(index, nbit(i));    
    for j = 1:length(peB)
        outdata = bsc(words, peB(j));
        outidx = bi2de(outdata);
        outidx = outidx';
        vout = codebook(outidx + 1);
        e = sigrec - vout;
        SNR(j) = snr(sigrec, e);
    end
    
    figure(2)
    subplot(3,2,2*(i-1)+1)
    line = ['SNR (voce registrata, ', num2str(nbit(i)), ' bit)'];
    semilogx(peA, 10*log10(SNRt));
    hold on
    grid on
    semilogx(peB, SNR, 'o');
    title(line);
    legend('SNR teorica', 'SNR segnale');
    xlabel('P_e')
    
    % Analisi file audio
    % Quantizzazione
    [index, quants] = quantiz(sigfile,partition,codebook);
        
    figure(2);
    
    % Spettro di potenza
    psd = abs(fft(quants)).^2;
    
    subplot(cases,2,2*(i-1)+1)
    grid on
    hold on
    plot(ff, fftshift(psd));
    xlabel('f')
    line = ['Spettro di potenza (file audio ', num2str(nbit(i)), ' bit)'];
    title(line);
    
    % Densità di probabilità
    subplot(cases,2,2*(i-1)+2)
    histogram(quants, M)
    line = ['Densità di probabilità (file audio ', num2str(nbit(i)), ' bit)'];
    title(line)
    hold on
    
    % Calcolo SNR
    words = de2bi(index, nbit(i));    
    for j = 1:length(peB)
        outdata = bsc(words, peB(j));
        outidx = bi2de(outdata);
        outidx = outidx';
        vout = codebook(outidx + 1);
        e = sig - vout;
        SNR(j) = snr(sigfile, e);
    end
    
    figure(2)
    subplot(3,2,2*(i-1)+2)
    line = ['SNR (file audio, ', num2str(nbit(i)), ' bit)'];
    semilogx(peA, 10*log10(SNRt));
    hold on
    grid on
    semilogx(peB, SNR, 'o');
    title(line);
    legend('SNR teorica', 'SNR segnale');
    xlabel('P_e')
end
