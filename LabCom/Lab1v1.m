clear all

% Calcolo grafico SMR

nbit = [4 6 8];
pe = logspace(-9,-1,1e5);

figure(1);

for i = 1:length(nbit)
    M = 2^nbit(i);
    SNR = M^2./(1+4*(M^2-1)*pe);
    semilogx(pe,10*log10(SNR));
    grid on;
    hold on;
end

title('SNR')
xlabel('P_e')
ylabel('SNR')

% Definizione segnale tramite rand

V = 1; % Ampiezza massima del segnale (da -V a +V)
fc = 10e4 % Frequenza di campionamento (scelta arbitraria)
Tc = 1/fc % Tempo di campionamento
samples = 1e4; % Numero di campioni del segnale (arbitrario)
t = (0:1:samples-1)*Tc; % Valori asse dei tempi dei campioni del segnale
sig = rand(1,samples); % Generazione del segnale
sig = (sig-0.5)*2*V; % Normalizzazione del segnale (media nulla e varianza pari a V)

% Quantizzazione e analisi del segnale

for i = 1:length(nbit)
    
    % Quantizzazione
    
    M = 2^nbit(i); % Numero intervalli di quantizzazione
    DV = 2*V/M; % Passo di quantizzazione
    partition = -V+DV:DV:V-DV; % Partizione asse delle ampiezze
    codebook = -V+DV/2:DV:V-DV/2; % Valori quantizzati
    
    % Quantizzazione
    [index, quants] = quantiz(sig,partition,codebook);
    
    % Codifica
    words = de2bi(index, nbit(i));
    
    % Spettro di potenza
    psd = abs(fft(quants)).^2;
    f = (linspace(-1,1,samples))*fc % Valori trasformata asse delle frequenze
    
    figure(2)
    hold on
    grid on
    plot(f, fftshift(psd))
    xlabel('f')
    ylabel('Sx')
    title('Spettro di potenza')
    
    % Densità di probabilità
    figure(3)
    hold on
    title('Densità di probabilità')
    histogram(quants, M)
end

