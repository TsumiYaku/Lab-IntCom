nbit = [4 6 8]; % Bit di quantizzazione
V = 1; % Ampiezza massima dei segnali

% Registrazione audio
fsr = 8000; % Frequenza di campionamento
Tr = 3; % Tempo di registrazione

% Registrazione
recorder = audiorecorder(fsr, 8, 1);
recordblocking(recorder, Tr);
sigrec = getaudiodata(recorder);

rsamples = length(sigrec); % Numero di campioni
fr = linspace(-fsr, fsr, rsamples); % Valori asse delle frequenze
    

% Apertura file audio
[sigfile, fsf] = audioread('/home/yaku/Desktop/record.wav');
fsamples = length(sigfile); % Numero di campioni
Tf = fsamples * 1/fsf; % Durata del segnale
ff = linspace(-fsf, fsf, fsamples); % Valori asse delle frequenze

for i = 1:length(nbit)
    
    % Definizione partizioni quantizzazione
    M = 2^nbit(i); % Numero intervalli di quantizzazione
    DV = 2*V/M; % Passo di quantizzazione
    partition = -V+DV:DV:V-DV; % Partizione asse delle ampiezze
    codebook = -V+DV/2:DV:V-DV/2; % Valori quantizzati
    
    % Analisi registrazione vocale
    % Quantizzazione
    [rindex, rquants] = quantiz(sigrec,partition,codebook);
    
    figure(1)
    
    % Spettro di potenza
    rpsd = abs(fft(sigrec)).^2;
    
    subplot(1,2,1);
    title('Spettro di potenza (voce registrata)');
    grid on
    hold on
    plot(fr, fftshift(rpsd));
    xlabel('f')
    
    % Densità di probabilità
    subplot(1,2,2)
    title('Densità di probabilità (voce registrata)')
    histogram(rquants, M)
    hold on
    
    % Analisi file audio
    % Quantizzazione
    [findex, fquants] = quantiz(sigfile,partition,codebook);
        
    figure(2);
    
    % Spettro di potenza
    fpsd = abs(fft(fquants)).^2;
    
    subplot(1,2,1)
    grid on
    hold on
    plot(ff, fftshift(fpsd));
    xlabel('f')
    title('Spettro di potenza (file audio)');
    
    % Densità di probabilità
    subplot(1,2,2)
    histogram(fquants, M)
    hold on
    title('Densità di probabilità (file audio)')
end
