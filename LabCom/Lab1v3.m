% Registrazione audio
fs = 8000; % Frequenza di campionamento
Tr = 5; % Tempo di registrazione

recorder = audiorecorder(fs, 16, 1);
recordblocking(recorder, Tr);
sigrec = getaudiodata(recorder); % Segnale registrato
sigrec = sigrec'; % Trasposto per compatibilità

% Apertura file audio
[sigfile, fsf] = audioread('record.wav');
sigfile = sigfile'; % Trasposto per compatibilità

% Calcolo delle SNR
nbit = [4 6 8]; % Bit di quantizzazione
V = 1; % Ampiezza massima dei segnali
mu = 255; % Coefficiente per compand
peA = logspace(-9, -1, 1e3); % Valori di probabilità (SNR teorica)
peB = logspace(-9, -1, 9); % Valori di probabilità (SNR segnali)

for i = 1:length(nbit)
    
    M = 2^nbit(i); % Numero intervalli di quantizzazione
    
    % SNR teorica
    SNRt = M^2./(1+4*(M^2-1)*peA);
    
    % Intervalli di quantizzazione uniforme (per companding)
    DV = 2*V/M; % Passo di quantizzazione
    partition = -V+DV:DV:V-DV; % Partizione asse delle ampiezze
    codebook = -V+DV/2:DV:V-DV/2; % Valori quantizzati
    
    % Quantizzazione registrazione vocale
    [lpartition, lcodebook] = lloyds(sigrec, M); % Partizioni algoritmo di Lloyd
    [lindex, lquants] = quantiz(sigrec,lpartition,lcodebook); % Quantizzazione di Lloyd
    comprec = compand(sigrec, mu, V); % Segnale compresso con compand
    [cindex, cquants] = quantiz(comprec,partition,codebook); % Quantizzazione con compand
    
    % Calcolo SNR
    lrSNR = fullSNR(sigrec, lindex, lcodebook, peB); % SNR con quantizzazione di Lloyd
    crSNR = fullSNR(comprec, cindex, codebook, peB); % SNR con compand
    
    % Quantizzazione file audio
    [lpartition, lcodebook] = lloyds(sigfile, M); % Partizioni algoritmo di Lloyd
    [lindex, lquants] = quantiz(sigfile,lpartition,lcodebook); % Quantizzazione di Lloyd
    compfile = compand(sigfile, mu, V); % Segnale compresso con compand
    [cindex, cquants] = quantiz(compfile,partition,codebook); % Quantizzazione con compand
    
    % Calcolo SNR
    lfSNR = fullSNR(sigfile, lindex, lcodebook, peB); % SNR con quantizzazione di Lloyd
    cfSNR = fullSNR(compfile, cindex, codebook, peB); % SNR con compand
    
    % Plot dei risultati
    figure(1)
    subplot(3,2,2*(i-1)+1)
    semilogx(peA,10*log10(SNRt));
    hold on
    grid on
    semilogx(peB, lrSNR, 'o')
    semilogx(peB, crSNR, 'o')
    line = ['SNR (voce registrata, ', num2str(nbit(i)), ' bit)'];
    title(line);
    legend('SNR teorica', 'SNR Lloyd', 'SNR Companding');
    xlabel('P_e');
    
    subplot(3,2,2*(i-1)+2)
    semilogx(peA,10*log10(SNRt));
    hold on
    grid on
    semilogx(peB, lfSNR, 'o')
    semilogx(peB, cfSNR, 'o')
    line = ['SNR (file audio, ', num2str(nbit(i)), ' bit)'];
    title(line);
    legend('SNR teorica', 'SNR Lloyd', 'SNR Companding');
    xlabel('P_e');
end

function SNR=fullSNR(sig, index, codebook, pe) % funzione per il calcolo di SNR su tutte le probabilità

SNR = zeros(1, length(pe));
indata = de2bi(index);

for i = 1:length(pe)
    outdata = bsc(indata, pe(i));
    outidx = bi2de(outdata);
    vout = codebook(outidx+1);
    e = sig - vout;
    SNR(i) = snr(sig, e);
end

end