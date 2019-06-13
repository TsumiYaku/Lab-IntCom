clear all
close all

%% Parametri

MPAM = 2; 
nbits = 1e4; % Numero di bit trasmessi
Mbps = 100; % Velocità di trasmissione in Mbps
fc = 16000; % Banda di simulazione (freq di campionamento)

BpS = log2(MPAM); % Numero di bits per simbolo
Rs = Mbps/BpS; % Velocità di trasmissione simboli in Baud
SpS = fc/Rs; % Campioni per simbolo
Ts = 1/Rs; % Tempo di trasmissione di un simbolo
Tc = 1/fc;
sym2alpha = [0;1]; % Tabella di conversione da simboli a coefficienti PAM
alpha2sym = [0;1]; % Tabella di conversione da coefficienti PAM a simboli
CarrierFreq = 2e3; % Frequenza carrier sinusoidale

%% Calcolo bits

% Matricole affiancate in vettore
matr = [2,1,7,2,0,4,2,3,5,2,0,6,2,2,6,6,1,9];

Bits = [];
for i=1:length(matr)
    Bits = [Bits, de2bi(matr(i),8)];
end

%% Generazione segnale elettrico

sig = toSig(Bits, sym2alpha, alpha2sym, SpS, BpS);
sig = sig';

%% Calcolo e riproduzione del segnale in uscita con carrier 

samples = length(sig); % Campioni del segnale
t = [0:1:samples-1]*Tc;

sig_out = sig.*cos(2*pi*CarrierFreq*t);
audiowrite('audio.wav', sig_out, fc);

%% Plot delle PSD (test)

df = fc/samples;
f = -fc/2:df:fc/2-df;

subplot(2,2,1), plot(f, abs(fftshift(fft(sig))).^2);
subplot(2,2,2), plot(f, abs(fftshift(fft(sig_out))).^2);


sig_out_rx = [zeros(1,1e4), sig_out, zeros(1,2e4)];
samples1 = length(sig_out_rx);
t1 = [0:1:samples1-1]*Tc;
sig_out_rx = sig_out_rx.*cos(2*pi*CarrierFreq*t1);
df1 = fc/samples1;
f1 = -fc/2:df1:fc/2-df1;
subplot(2,2,3), plot(f1, abs(fftshift(fft(sig_out_rx))).^2);

H_matched = sinc(f1*Ts);
sig_out_rx = sigFilter(sig_out_rx, H_matched);
subplot(2,2,4), plot(f1, abs(fftshift(fft(sig_out_rx))).^2);

H_matched = sinc(f*Ts);
signal_ref = sigFilter(sig, H_matched);
sig_out_rx = sigAlign(signal_ref, sig_out_rx);
sig_out_rx = sigNorm(sig_out_rx, signal_ref);

bits_out = toBits(sig_out_rx, sym2alpha, alpha2sym, SpS);

BER = ber(Bits, bits_out);

eyediagram(sig_out_rx,SpS*2,SpS*2);
eyediagram(signal_ref,SpS*2,SpS*2);

%% Funzioni di utilità

function sig = toSig(bits, sym2alpha, alpha2sym, SpS, BpS)
    sig = zeros(length(bits)*SpS/BpS,1);
    sym = zeros(1, BpS);
    
    for i=1:BpS:length(bits)
        % Salva il simbolo singolo per convertirlo
        for j = 1:BpS
            sym(j) = bits(i+j-1);
        end
        
        % Converte il simbolo in coefficiente del segnale
        for j=1:length(sym2alpha)
            if(isequal(sym, alpha2sym(j,:)))
                for k=1:SpS
                    sig(((i-1)*SpS/BpS)+k) = sym2alpha(j);
                end
                break;
            end
        end
    end
end

% Converte un segnale elettrico in bits
function bits = toBits(sig, sym2alpha, alpha2sym, SpS)
bits = [];

% Creazione delle soglie 
th = [];
for i=1:length(sym2alpha)-1
    th = [th, (sym2alpha(i+1)+sym2alpha(i))/2]; 
end

% Comprime il segnale a 1 SpS
alphas = []; % Coefficienti
for i=0:length(sig)/SpS-1
    alphas = [alphas sig(i*SpS+ceil(SpS/2))];
end

% Converte i coefficienti in simboli
for i=1:length(alphas)
    for j=1:length(th)
        sym = alpha2sym(j+1, :);
        if(abs(alphas(i)) <= th(j)) 
            sym = alpha2sym(j, :);
            break;
        end
    end
    bits = [bits, sym];
end
end

% Calcola la BER tra segnale entrante e segnale uscente
function BER = ber(bits_in, bits_out)
t = abs(bits_in-bits_out); % Calcola i bits di differenza tra i due segnali
BER = sum(t)/length(t);
end

% Filtra il segnale con il filtro passato come parametro
function sig_out = sigFilter(sig, filter) 
Sig_in = fft(sig);
Sig_out = Sig_in.*fftshift(filter);
sig_out = real(ifft(Sig_out));
end

function signal_out = sigAlign(signal_ref,signal_rx) 
[signal_ref_aligned,signal_rx_aligned,Delay] = alignsignals(signal_ref,signal_rx);
signal_out=signal_rx(Delay+1:Delay+length(signal_ref));
end

function sig_out = sigNorm(sig, ref)

pow_sig = mean(abs(sig.^2));
pow_ref = mean(abs(ref.^2));

sig_out = sig.*sqrt(pow_ref)./sqrt(pow_sig);

end