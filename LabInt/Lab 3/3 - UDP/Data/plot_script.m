filename = 'Data/Data.txt';
delimiterIn = '\t';
headerlinesIn = 1;
Data = importdata(filename,delimiterIn,headerlinesIn);

t = Data.data(:,1);
ports = Data.data(:,2);

figure; plot(t, ports, '.'), xlabel('Time [s]'), ylabel('Port'), title('Ports scanned over time')
