clear all
close all

%% Parametri
MPAM = 4; 
nbits = 1e4; % Numero di bit trasmessi
Mbps = 100; % Velocità di trasmissione in Mbps
SpS = 16; % Campioni per simbolo

a = 2; b = 15; % dB iniziali e finali Eb/N0
EbN0_sim_db = a:1:b; % Eb/N0 per simulazione (in dB)
EbN0_sim = 10.^(EbN0_sim_db./10); % Eb/N0 per simulazione (in valore assoluto)
EbN0_theory_db = a:1e-5:b; % Eb/N0 per calcolo teorico (in dB)
EbN0_theory = 10.^(EbN0_theory_db./10); % Eb/N0 per calcolo teorico (in valore assoluto)

BpS = log2(MPAM); % Numero di bits per simbolo
Rs = Mbps/BpS; % Velocità di trasmissione simboli in Baud
Ts = 1/Rs; % Tempo di trasmissione di un simbolo
fs = SpS*Rs; % Banda di simulazione (freq di campionamento)
sym2alpha = [-3;-1;1;3]; % Tabella di conversione da simboli a coefficienti PAM
alpha2sym = [ % Tabella di conversione da coefficienti PAM a simboli
    0,0;
    0,1;
    1,1;
    1,0];

%% Generazione segnale elettrico

Bits = randi([0 1], [nbits, 1]); % Bits casuali
Bits = Bits';

sig_in = toSig(Bits, sym2alpha, alpha2sym, SpS, BpS);
sig_in = sig_in';

%% Creazione filtri

% Campioni del segnale
samples = length(sig_in); 

% Valore asse delle frequenze
f = linspace(-fs/2, fs/2, samples); 

% Filtro adattato
H_matched = sinc(f*Ts);

% Filtri RC
poles = [50, 100, 150];
H_pole = [
    1./(1+(1i*f/poles(1)));
    1./(1+(1i*f/poles(2)));
    1./(1+(1i*f/poles(3)));
    ];

%% Calcolo delle BER per diversi Eb/N0

% Curva BER teorica
BER_theory = 3/8*erfc(sqrt(2/5*EbN0_theory)); 

% Potenza del segnale
psig = mean(abs(sig_in).^2); 

BER_matched = zeros(1,length(EbN0_sim));
BER_poles = zeros(length(poles), length(EbN0_sim));

for i=1:length(EbN0_sim)
    
    % Calcolo AWGN
    sigma = (psig*SpS)/(2*EbN0_sim(i)*BpS);
    sigma = sqrt(sigma);
    noise = sigma*randn(1,samples);
    
    sign = sig_in+noise;
    
    % Filtraggio segnali
    sig_out_matched = sigFilter(sign, H_matched);
    sig_out_pole = [
        sigFilter(sign, H_pole(1,:));
        sigFilter(sign, H_pole(2,:));
        sigFilter(sign, H_pole(3,:));
        ];
    
    % Conversione da segnali a bits
    bits_matched = toBits(sig_out_matched, sym2alpha, alpha2sym, SpS);
    bits_pole = [
        toBits(sig_out_pole(1,:), sym2alpha, alpha2sym, SpS);
        toBits(sig_out_pole(2,:), sym2alpha, alpha2sym, SpS);
        toBits(sig_out_pole(3,:), sym2alpha, alpha2sym, SpS);
    ];
    
    % Calcolo BER segnali
    BER_matched(i) = ber(Bits, bits_matched);
    BER_poles(1,i) = ber(Bits, bits_pole(1,:));
    BER_poles(2,i) = ber(Bits, bits_pole(2,:));
    BER_poles(3,i) = ber(Bits, bits_pole(3,:));
    
end

%% Plot delle BER

figure(1)
semilogy(EbN0_theory_db, BER_theory);
hold on
semilogy(EbN0_sim_db, BER_matched, 'o');
semilogy(EbN0_sim_db, BER_poles(1,:), 'o');
semilogy(EbN0_sim_db, BER_poles(2,:), 'o');
semilogy(EbN0_sim_db, BER_poles(3,:), 'o');
legend('SNR teorica', 'SNR f. adattato', 'SNR f. 50MHz', 'SNR f. 100 Mhz', 'SNR f. 150 MHz')
grid on

%% Plot diagrammi ad occhio

% Segnali filtrati senza rumore
sig_out_matched = sigFilter(sig_in, H_matched);
sig_out_pole = [
    sigFilter(sig_in, H_pole(1,:));
    sigFilter(sig_in, H_pole(2,:));
    sigFilter(sig_in, H_pole(3,:));
    ];

eyediagram(sig_out_matched,SpS*2,SpS*2);
eyediagram(sig_out_pole(1,:),SpS*2,SpS*2);
eyediagram(sig_out_pole(2,:),SpS*2,SpS*2);
eyediagram(sig_out_pole(3,:),SpS*2,SpS*2);

%% Funzioni di utilità

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
        if(alphas(i) <= th(j)) 
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

% Converte bits in segnale elettrico
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