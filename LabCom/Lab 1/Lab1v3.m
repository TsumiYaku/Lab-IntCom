%% Registrazione audio

fs = 8000; % Frequenza di campionamento
Tr = 5; % Tempo di registrazione

recorder = audiorecorder(fs, 16, 1);
recordblocking(recorder, Tr);
sigrec = getaudiodata(recorder); % Segnale registrato
sigrec = sigrec'; % Trasposto per compatibilità

%% Apertura file audio

[sigfile, fsf] = audioread('record.wav');
sigfile = sigfile'; % Trasposto per compatibilità

%% Calcolo delle SNR
nbit = [4 6 8]; % Bit di quantizzazione
V = 1; % Ampiezza massima dei segnali
mu = 255; % Coefficiente per compand
pe_theory = logspace(-9, -1, 1e3); % Valori di probabilità (SNR teorica)
pe_sig = logspace(-9, -1, 9); % Valori di probabilità (SNR segnali)

for i = 1:length(nbit)
    
    M = 2^nbit(i); % Numero intervalli di quantizzazione
    
    % SNR teorica
    SNRt = M^2./(1+4*(M^2-1)*pe_theory);
    
    % Intervalli di quantizzazione uniforme (per companding)
    DV = 2*V/M; % Passo di quantizzazione
    partition = -V+DV:DV:V-DV; % Partizione asse delle ampiezze
    codebook = -V+DV/2:DV:V-DV/2; % Valori quantizzati
    
    % Quantizzazione registrazione vocale
    [index, quants] = quantiz(sigrec, partition, codebook);
    [lpartition, lcodebook] = lloyds(sigrec, M); % Partizioni algoritmo di Lloyd
    [lindex, lquants] = quantiz(sigrec, lpartition, lcodebook); % Quantizzazione di Lloyd
    comprec = compand(sigrec, mu, V); % Segnale compresso con compand
    [cindex, cquants] = quantiz(comprec, partition, codebook); % Quantizzazione con compand
    
    % Calcolo SNR
    rSNR = fullSNR(sigrec, index, codebook, pe_sig);
    lrSNR = fullSNR(sigrec, lindex, lcodebook, pe_sig); % SNR con quantizzazione di Lloyd
    crSNR = zeros(1, length(pe_sig)); % SNR con compand
    indata = de2bi(cindex);
    
    for j = 1:length(pe_sig)
        outdata = bsc(indata, pe_sig(j));
        outidx = bi2de(outdata);
        vout = codebook(outidx+1);
        vout = compand(vout, mu, V, 'mu/expander');
        e = sigrec - vout;
        crSNR(j) = snr(sigrec, e);
    end
    
    % Quantizzazione file audio
    [index, quants] = quantiz(sigfile, partition, codebook);
    [lpartition, lcodebook] = lloyds(sigfile, M); % Partizioni algoritmo di Lloyd
    [lindex, lquants] = quantiz(sigfile, lpartition, lcodebook); % Quantizzazione di Lloyd
    compfile = compand(sigfile, mu, V); % Segnale compresso con compand
    [cindex, cquants] = quantiz(compfile, partition, codebook); % Quantizzazione con compand
    
    % Calcolo SNR
    fSNR = fullSNR(sigfile, index, codebook, pe_sig);
    lfSNR = fullSNR(sigfile, lindex, lcodebook, pe_sig); % SNR con quantizzazione di Lloyd
    cfSNR = zeros(1, length(pe_sig)); % SNR con compand
    indata = de2bi(cindex);
    
    for j = 1:length(pe_sig)
        outdata = bsc(indata, pe_sig(j));
        outidx = bi2de(outdata);
        vout = codebook(outidx+1);
        vout = compand(vout,mu,V, 'mu/expander');
        e = sigfile - vout;
        cfSNR(j) = snr(sigfile, e);
    end
    
    % Plot dei risultati
    figure(1)
    subplot(3,1,i)
    semilogx(pe_theory,10*log10(SNRt));
    hold on
    grid on
    semilogx(pe_sig, rSNR, 'o')
    semilogx(pe_sig, lrSNR, 'o')
    semilogx(pe_sig, crSNR, 'o')
    line = ['SNR (voce registrata, ', num2str(nbit(i)), ' bit)'];
    title(line);
    legend('SNR teorica', 'SNR Uniforme', 'SNR Lloyd', 'SNR Companding');
    xlabel('P_e');
    
    figure(2)
    subplot(3,1,i)
    semilogx(pe_theory,10*log10(SNRt));
    hold on
    grid on
    semilogx(pe_sig, fSNR, 'o')
    semilogx(pe_sig, lfSNR, 'o')
    semilogx(pe_sig, cfSNR, 'o')
    line = ['SNR (file audio, ', num2str(nbit(i)), ' bit)'];
    title(line);
    legend('SNR teorica', 'SNR Uniforme', 'SNR Lloyd', 'SNR Companding');
    xlabel('P_e');
end

%% Funzioni di utilità

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