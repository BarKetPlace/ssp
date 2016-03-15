clear all
close all
clc

load assignment2.mat

in = male8;%[-6:.01:6]'; %Signal must be a columns

nbits_tab = 3;%:3; %Loop on the element of this array


k_ = [0.1:0.5:2 2:0.1:5 5:.5:11];
m_tab=[0];% 1.5]; %Offset of the quantizer reconstruction levels =0 ->midrise quantizer =delta/2 ->midtread
D = zeros(length(k_),length(m_tab));%Distorsion
PSNR = zeros(length(nbits_tab),length(m_tab)); 
IN = ones(length(k_),1)*in';%contains the input in rows
OUT = zeros(length(k_), length(in));%Output of USQ with nbits specified in nbits_tab

i_bits =1;
n_bits = nbits_tab(i_bits);

en_plots = 0; %Enables plot in functions
for ik=1:length(k_)
    xmax = k_(ik)*sqrt(var(in));
    m = m_tab(1);
   
    fprintf('%d bits\n', n_bits);

    %Find the index of the samples
    [ idx ] = sq_enc(in, n_bits, xmax, m, en_plots);
    %Convert the index into actual values
    OUT(ik,:) = sq_dec(idx, n_bits, xmax, m);
end


%Computes the distorsion
D = 1/length(in)*sum( (IN-OUT).^2,2);
SNR = var(in-mean(in))./D;


[M, imax] = max(SNR);

figure, %Plot Rate vs distorsion curve with m_tab of dimension 2
plot(k_, SNR,'LineWidth',2);hold on;
plot(k_(imax)*[1 1], [min(SNR) M]);
% plot(nbits_tab, D(:,2),'LineWidth',2);
% legend(['m= ' num2str(m_tab(1))] ,['m= ' num2str(m_tab(2))]);
xlabel('Value of k. xmax = k*sqrt(var(in))'); ylabel('SNR');
 title({'USQ'; ['SNR = f(k) with R = ' num2str(n_bits)]; ['>Optimal k  = ' num2str(k_(imax))  '<']});
 
% soundsc(in,8000);
% soundsc(OUT(3,:),8000);

% figure, plot(in,outq);
% title(['Uniform Scalar Quantizer. ' num2str(n_bits) ' bits, m= ' num2str(m) ', xmax= ' num2str(xmax) ]);
% xlabel('Input');ylabel('Output');

%% part2



clear all
close all
clc

load assignment2.mat

in = male8;%[-6:.01:6]'; %Signal must be a columns

nbits_tab = 1:11; %Loop on the element of this array
IN = ones(length(nbits_tab),1)*in';%contains the input in rows
OUT = zeros(length(nbits_tab), length(in));%Output of USQ with nbits specified in nbits_tab


k_ = [0.95, 2.1, 3.4, 4.95, 6.3, 7.65, 8.85, 9.95, 10.6, 11.0, 11.1, 11.2, 11.15, 11.2,...
11.15, 11.15];
m_tab=[0];% 1.5]; %Offset of the quantizer reconstruction levels =0 ->midrise quantizer =delta/2 ->midtread
D = zeros(length(nbits_tab),length(m_tab));%Distorsion


en_plots = 1; %Enables plot in functions

for im = 1:length(m_tab)%For all offset
 
    
    for i_bits = 1:length(nbits_tab)%For all bit rate we want to try
        xmax = k_(i_bits)*sqrt(var(in));
        n_bits = nbits_tab(i_bits);
        fprintf('%d bits\n', n_bits);
        m = xmax/2.^n_bits;
        
        % n_bits = 2;
        %Find the index of the samples
        [ idx ] = sq_enc(in, n_bits, xmax, m, en_plots);
        
        %Convert the index into actual values
%         OUT(i_bits,:) = sq_dec(idx, n_bits, xmax, m);
           OUT(i_bits,:) = sq_dec(idx, n_bits, xmax, m);
    end
    
    %Computes the distorsion
    D(:,im) = 1/length(in)*sum( (IN-OUT).^2,2);
    
end
ERROR=IN-OUT;
SNR = 10*log10(var(in)./D);
figure, %Plot Rate vs distorsion curve with m_tab of dimension 2
plot(nbits_tab, SNR,'LineWidth',2);hold on;
% plot(nbits_tab, D(:,2),'LineWidth',2);
% legend(['m= ' num2str(m_tab(1))] ,['m= ' num2str(m_tab(2))]);
xlabel('Rate (bits)'); ylabel('SNR (dB)');
title({'USQ'; 'Rate vs SNR'});
 
% soundsc(in,8000);
%% Listen to quantized
soundsc(OUT(1,:),8000);
%% Listen to Error
soundsc(ERROR(11,:),8000);

% figure, plot(in,outq);
% title(['Uniform Scalar Quantizer. ' num2str(n_bits) ' bits, m= ' num2str(m) ', xmax= ' num2str(xmax) ]);
% xlabel('Input');ylabel('Output');

