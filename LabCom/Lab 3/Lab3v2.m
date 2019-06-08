clear all
close all

%% Parametri

SampleRate=2e6;
NsBlock=1024; % Bits per blocco
Nblocks=1000; % Numero di blocchi
Div = 1; % Numero di blocchi per ogni PSD parziale
StartFreq = 210e6; % Frequenza iniziale di cattura
EndFreq = 226e6; % Frequenza finale di cattura

% Centri di frequenza
FreqCenters = StartFreq:SampleRate:EndFreq;

%% Cattura segnale con scheda RTL

failed = false; % Utilizzata per determinare se la chiavetta funziona

Sigs = NaN*ones(Nblocks*NsBlock, length(FreqCenters));

for i=1:length(FreqCenters)
    CenterFrequency=FreqCenters(i);
    
    hRadio = comm.SDRRTLReceiver('CenterFrequency', CenterFrequency, ...
        'SampleRate', SampleRate, ...
        'EnableTunerAGC', true, ...
        'SamplesPerFrame', NsBlock, ...
        'OutputDataType', 'single');
    
    if ~isempty(sdrinfo(hRadio.RadioAddress))
        x_saved=NaN*ones(Nblocks*NsBlock,1);
        lost = NaN*ones(Nblocks);
        for Counter=1:Nblocks
            [x, len, lost(Counter)] = step(hRadio);
            x_saved((Counter-1)*NsBlock+1:(Counter-1)*NsBlock+NsBlock)=x;
        end
        Sigs(:,i) = x_saved;
    else
        warning('SDR Device not connected')
        failed = true;
        break;
    end
    
    release(hRadio)
end

% Utilizza l'ultima cattura da file come segnale se la chiavetta non è connessa 
% (utile per lavorare quando non si ha la chiavetta a disposizione)
if(failed)
    load('Sigs.mat');
else
    save('Sigs.mat', Sigs);
end

%% Calcolo delle medie delle PSD di tutti i segnali

PSDs = NaN*ones(NsBlock*Div,length(FreqCenters)); % PSD dei vari segnali

for i = 1:length(FreqCenters)
    
    Sig = Sigs(:,i);
    
    % Divisione e calcolo delle PSD parziali
    PartPSDs = NaN*ones(NsBlock*Div, Nblocks/Div);
    for j = 1:Nblocks/Div
        portion = Sig((j-1)*NsBlock*Div+1:(j-1)*NsBlock*Div+NsBlock*Div);
        PartPSDs(:,j) = abs(fftshift(fft(portion)).^2);
    end
    
    % Media delle PSD parziali
    for j = 1:NsBlock*Div
        PSDs(j,i) = mean(PartPSDs(j,:));
    end
end

%% Unione delle PSD dei vari segnali

PSD = []; % PSD finale
ftot = []; % asse delle frequenze finale
Df = SampleRate/(NsBlock*Div);
f = -SampleRate/2:Df:SampleRate/2-Df;

for i = 1:length(FreqCenters)
    ftot = [ftot f+FreqCenters(i)];
    PSD = [PSD PSDs(:,i)'];
end

%% Plot delle PSD

semilogy(ftot,PSD);

grid on
xlim([FreqCenters(1)-SampleRate FreqCenters(length(FreqCenters))+SampleRate])
xticks(FreqCenters(1)-SampleRate:SampleRate:FreqCenters(length(FreqCenters))+SampleRate)
ylim([0 max(PSD)+0.1*max(PSD)])