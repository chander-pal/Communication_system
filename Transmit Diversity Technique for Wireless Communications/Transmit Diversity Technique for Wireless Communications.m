clear all
clc
% constants
N = 1e7;% No. of bits Transmitted
snr= 0:1:50; % Eb/N0 vector in db from 0 to 50
Message_Sequence= sign(randi(2,1,N)-1.5);%randi generates sequence of 1's and 2's the -1.5 and sign maps it to 1 and -1 
N0 = -snr ;% It is assumed that P = 1 watt, then for Eb/N0 to vary from 0 to 10
                         % N0 varies from 0 to -40 dB, the power of the
                         % Gaussian Pdf  is N0 dB.
 format long
 %% traditional MRRC with one transmitter and two reciever
ber_no_diversity= zeros(1,length(snr)); %Initialising BER Vector

h = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N));

for K = 1:length(N0)       
       faded_sequence1 = h(1,:).*Message_Sequence; 
       n = 1/sqrt(2)*(randn(2,N) + 1i*randn(2,N));
       Received_Sequence_rayl = faded_sequence1+ 10^(N0(K)/20).*n(1,:);
       
       combined_seq =(conj(h(1,:)).*Received_Sequence_rayl);%./(conj(h(1,:)).*h(1,:) + conj(h(2,:)).*h(2,:));

       decoded_seq = sign(real(combined_seq));
       
       ber_no_diversity(K) = sum(abs(decoded_seq-Message_Sequence))/(2*N);
       
end 
%% traditional MRRC with one transmitter and two reciever
ber_1_tx_2_rx= zeros(1,length(snr)); %Initialising BER Vector

h = 1/sqrt(2)*(randn(2,N) + 1i*randn(2,N));

for K = 1:length(N0)       
       faded_sequence1 = h(1,:).*Message_Sequence; 
       n = 1/sqrt(2)*(randn(2,N) + 1i*randn(2,N));
       Received_Sequence_rayl1 = faded_sequence1+ 10^(N0(K)/20).*n(1,:);
       faded_sequence2 = h(2,:).*Message_Sequence; 
       Received_Sequence_rayl2 = faded_sequence2+ 10^(N0(K)/20).*n(2,:);
       
       combined_seq =(conj(h(1,:)).*Received_Sequence_rayl1 + conj(h(2,:)).*Received_Sequence_rayl2);%./(conj(h(1,:)).*h(1,:) + conj(h(2,:)).*h(2,:));

       decoded_seq = sign(real(combined_seq));
       
       ber_1_tx_2_rx(K) = sum(abs(decoded_seq-Message_Sequence))/(2*N);
       
end 
%% the new transmit diversity scheme with two transmitter and one reciever

Tx1 = zeros(1,N);
Tx2 = Tx1;
k=1;
while k<N  
    Tx1(k)= Message_Sequence(k);
    Tx1(k+1) = -conj(Message_Sequence(k+1));
    Tx2(k) = -conj(Tx1(k+1));
    Tx2(k+1) = Tx1(k);
    k= k+2;
end 
ber_2_tx = zeros(1,length(N0));
h_general = 1/sqrt(2)*(randn(2,N/2) + 1i*randn(2,N/2));
h = zeros(2,N);
     o=1;
     while o<=N/2
         h(:,2*o-1) = h_general(:,o);
         h(:,2*o) = h(:,2*o-1);    
         o=o+1;
     end  
for k =1:length(N0)  
     n = 10^(N0(k)/20).*1/sqrt(2)*(randn(2,N) + 1i*randn(2,N));
     Rx1 = h(1,:).*Tx1 + n(1,:);
     Rx2 = h(2,:).*Tx2 + n(2,:);
     R = (Rx1 + Rx2);
     decide = zeros(1,N);
     for t=1:N/2
         decide(2*t-1)= conj(h_general(1,t))*R(2*t-1) + h_general(2,t)*conj(R(2*t));
         decide(2*t)= conj(h_general(2,t))*(R(2*t-1)) - conj(R(2*t))*h_general(1,t);
     end
     decoded_seq = sign(real(decide));
     ber_2_tx(k) = sum(abs(Message_Sequence - decoded_seq))./(2*N);
end



%% tx1_rx4
ber_1_tx_4_rx = zeros(1,length(N0));
h = 1/sqrt(2)*(randn(4,N) + 1i*randn(4,N));   
for k =1:length(N0)  
     n = 10^(N0(k)/20).*1/sqrt(2)*(randn(4,N) + 1i*randn(4,N));
     Rx11 = h(1,:).*Message_Sequence + n(1,:);
     Rx12 = h(2,:).*Message_Sequence + n(2,:);
     Rx13 = h(3,:).*Message_Sequence + n(3,:);
     Rx14 = h(4,:).*Message_Sequence + n(4,:);
     R = [Rx11;Rx12;Rx13;Rx14];
     combined_seq = sum(conj(h).*R,1);
     decoded_seq = sign(real(combined_seq));
     ber_1_tx_4_rx(k) = sum(abs(Message_Sequence - decoded_seq))./(2*N);
end

%% tx2_rx2
%Transmit Diversity 2Tx 2Rx

Tx1 = zeros(1,N);
Tx2 = Tx1;
k=1;
while k<N
    
    Tx1(k)= Message_Sequence(k);
    Tx1(k+1) = -conj(Message_Sequence(k+1));
    Tx2(k) = -conj(Tx1(k+1));
    Tx2(k+1) = Tx1(k);
    k= k+2;
end 

ber_2_tx_2_rx = zeros(1,length(N0));
h_general = 1/sqrt(2)*(randn(4,N/2) + 1i*randn(4,N/2));
     h = zeros(4,N);
     o=1;
     while o<=N/2
         h(:,2*o-1) = h_general(:,o);
         h(:,2*o) = h(:,2*o-1);    
         o=o+1;
     end  
for k =1:length(N0) 
     n = 10^(N0(k)/20).*1/sqrt(2)*(randn(4,N) + 1i*randn(4,N));
     Rx11 = h(1,:).*Tx1 + n(1,:);
     Rx21 = h(2,:).*Tx2 + n(2,:);
     Rx12 = h(3,:).*Tx1 + n(3,:);
     Rx22 = h(4,:).*Tx2 + n(4,:);
     R1 = Rx11 + Rx21;
     R2 = Rx12 + Rx22;
     decide = zeros(1,N);
     for t=1:N/2
         decide(2*t-1)= conj(h_general(1,t))*R1(2*t-1) + h_general(2,t)*conj(R1(2*t)) +conj(h_general(3,t))*R2(2*t-1) + h_general(4,t)*conj(R2(2*t));
         decide(2*t)= conj(h_general(2,t))*(R1(2*t-1)) - conj(R1(2*t))*h_general(1,t) + conj(h_general(4,t))*(R2(2*t-1)) - conj(R2(2*t))*h_general(3,t);
     end
     decoded_seq = sign(real(decide));
     ber_2_tx_2_rx(k) = sum(abs(Message_Sequence - decoded_seq))./(2*N);
end

%% plots
%Plot
snr = 0:1:50;
semilogy(snr,ber_no_diversity,'bo-')
hold on
semilogy(snr,ber_1_tx_2_rx,'mv-')
semilogy(snr,ber_2_tx,'rd-')
semilogy(snr,ber_1_tx_4_rx,'ks-')
semilogy(snr,ber_2_tx_2_rx,'-^g')
 axis([0 50 10^-6 0.5])
legend('No Diversity','1 Tx,2 Rx','2 Tx, 1 Rx','1 Tx, 4  Rx','2 Tx, 2 Rx')

grid on
title('Error Performance Analysis')
xlabel('Eb/N0 , dB')
ylabel('BER')