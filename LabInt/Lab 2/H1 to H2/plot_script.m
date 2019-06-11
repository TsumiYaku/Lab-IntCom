filename = 'Data.txt';
delimiterIn = '\t';
headerlinesIn = 1;
Data = importdata(filename,delimiterIn,headerlinesIn);

t = Data.data(:,1);
Seq = Data.data(:,2);
Ack = Data.data(:,3);
Win = Data.data(:,4);
Len = Data.data(:,5);

fig = figure; plot(t, Seq, '.'), xlabel('Time [s]'), ylabel('Seq N°'), title('Sequence number'), saveas(fig,'Plots/Seq.png');
fig = figure; plot(t, Ack, '.'), xlabel('Time [s]'), ylabel('Ack N°'), title('Ack Number'), saveas(fig,'Plots/Ack.png');
fig = figure; plot(t, Win, '.'), xlabel('Time [s]'), ylabel('Win size'), title('Congestion window'), saveas(fig,'Plots/Win.png');
fig = figure; plot(t, Len, '.'), xlabel('Time [s]'), ylabel('Packet len'), title('Packet length'), saveas(fig,'Plots/Len.png');
