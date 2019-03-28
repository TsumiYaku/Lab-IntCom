nbit = [4 6 8]; % Bit di quantizzazione
cases = length(nbit);
V = 1; % Ampiezza massima dei segnali
mu = 255;

% Registrazione audio
fsr = 8000; % Frequenza di campionamento
Tr = 5; % Tempo di registrazione

recorder = audiorecorder(fsr, 16, 1);
recordblocking(recorder, Tr);
sigrec = getaudiodata(recorder); % Segnale registrato

rsamples = length(sigrec); % Numero di campioni
fr = linspace(-fsr/2, fsr/2, rsamples); % Valori asse delle frequenze
    

% Apertura file audio
[sigfile, fsf] = audioread('/home/yaku/Desktop/record.wav');
fsamples = length(sigfile); % Numero di campioni
Tf = fsamples * 1/fsf; % Durata del segnale
ff = linspace(-fsf/2, fsf/2, fsamples); % Valori asse delle frequenze

for i = 1:length(nbit)
    
    M = 2^nbit(i); % Numero intervalli di quantizzazione
    
    % Intervalli di quantizzazione uniforme (per companding)
    DV = 2*V/M; % Passo di quantizzazione
    partition = -V+DV:DV:V-DV; % Partizione asse delle ampiezze
    codebook = -V+DV/2:DV:V-DV/2; % Valori quantizzati
    
    % Quantizzazione registrazione vocale
    
    [lpartition, lcodebook] = lloyds(sigrec, M);
    [lrindex, lrquants] = quantiz(sigrec,lpartition,lcodebook);
    comprec = compand(sigrec, mu, V);
    [crindex, crquants] = quantiz(comprec,partition,codebook);
    
    % Quantizzazione file audio
    [lpartition, lcodebook] = lloyds(sigfile, M);
    [lfindex, lfquants] = quantiz(sigfile,lpartition,lcodebook);
    compfile = compand(sigfile, mu, V);
    [cfindex, cfquants] = quantiz(compfile,partition,codebook);
    
    % Calcolo delle SNR
    SNRt = M^2./(1+4*(M^2-1)*pe);
    lrSNR = zeros(1, length(pe));
    crSNR = zeros(1, length(pe));
    lfSNR = zeros(1, length(pe));
    cfSNR = zeros(1, length(pe));
    
    
    % Plot dei risultati
    figure(1)
    subplot(3,2,2*(i-1)+1)
    plot(pe,[SNRt, lrSNR, crSNR]);
    line = ['SNR (voce registrata, ', num2str(nbit(i)), ' bit)'];
    title(line);
    legend('SNR teorica', 'SNR Lloyd', 'SNR Companding');
    xlabel('P_e');
    
    subplot(3,2,2*(i-1)+2)
    plot(pe,[SNRt, lfSNR, cfSNR]);
    line = ['SNR (file audio, ', num2str(nbit(i)), ' bit)'];
    title(line);
    legend('SNR teorica', 'SNR Lloyd', 'SNR Companding');
    xlabel('P_e');
end