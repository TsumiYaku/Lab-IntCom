nbit = [4 6 8]; % Bit di quantizzazione
V = 1;

% Registrazione audio
fsr = 8000; % Frequenza di campionamento
Tr = 3; % Tempo di registrazione
recorder = audiorecorder(fs, 8, 1);
recordblocking(recorder, Tr);
sigrec = getaudiodata(recorder);
rsamples = length(sigrec); % Numero di campioni
fr = linspace(-fsr, fsr, rsamples); % Valori asse delle frequenze
    

% Apertura file audio
[sigfile, fsf] = audioread('/home/yaku/Desktop/record.wav');
fsamples = length(sigfile);
Tf = samples * 1/fsf;
ff = linspace(-fsf, fsf, fsamples);

for i = 1:length(nbit)
    
    % Definizione partizioni quantizzazione
    M = 2^nbit(i); % Numero intervalli di quantizzazione
    DV = 2*V/M; % Passo di quantizzazione
    partition = -V+DV:DV:V-DV; % Partizione asse delle ampiezze
    codebook = -V+DV/2:DV:V-DV/2; % Valori quantizzati
    
    % Analisi registrazione vocale
    % Quantizzazione
    [index, rquants] = quantiz(sigrec,partition,codebook);
    
    % Spettro di potenza
    rpsd = abs(fft(sigrec)).^2;
    
    figure(1)
    title('Spettro di potenza (voce registrata)');
    grid on
    hold on
    plot(fr, fftshift(rpsd));
    xlabel('f')
    
    % Densità di probabilità
    figure(2)
    title('Densità di probabilità (voce registrata)')
    histogram(rquants, M)
    hold on
    
    % Analisi file audio
    % Quantizzazione
    [index, fquants] = quantiz(sigfile,partition,codebook);
        
    % Spettro di potenza
    fpsd = abs(fft(fquants)).^2;
    
    figure(3);
    grid on
    hold on
    plot(f, fftshift(fpsd));
    xlabel('f')
    title('Spettro di potenza (file audio)');
    
    % Densità di probabilità
    figure(4)
    histogram(fquants, M)
    hold on
    title('Densità di probabilità (file audio)')
end
