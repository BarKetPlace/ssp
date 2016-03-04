clear all
close all
clc

load assignment2.mat

in = male8;%[-6:.01:6]'; %Signal must be a columns

nbits_tab = 2:10; %Loop on the element of this array
IN = ones(length(nbits_tab),1)*in';%contains the input in rows
OUT = zeros(length(nbits_tab), length(in));%Output of USQ with nbits specified in nbits_tab

xmax = max(in);
m_tab=[0];% 1.5]; %Offset of the quantizer reconstruction levels =0 ->midrise quantizer =delta/2 ->midtread
D = zeros(length(nbits_tab),length(m_tab));%Distorsion
PSNR = zeros(length(nbits_tab),length(m_tab)); 




en_plots = 0; %Enables plot in functions

for im = 1:length(m_tab)%For all offset
    m = m_tab(im);
    
    for i_bits = 1:length(nbits_tab)%For all bit rate we want to try
        n_bits = nbits_tab(i_bits);
        fprintf('%d bits\n', n_bits);
        
        % n_bits = 2;
        %Find the index of the samples
        [ idx ] = sq_enc(in, n_bits, xmax, m, en_plots);
        
        %Convert the index into actual values
        OUT(i_bits,:) = sq_dec(idx, n_bits, xmax, m);
    end
    
    %Computes the distorsion
    D(:,im) = 1/length(in)*sum( (IN-OUT).^2,2);
    
end

figure, %Plot Rate vs distorsion curve with m_tab of dimension 2
plot(nbits_tab, D(:,1),'LineWidth',2);hold on;
% plot(nbits_tab, D(:,2),'LineWidth',2);
% legend(['m= ' num2str(m_tab(1))] ,['m= ' num2str(m_tab(2))]);
xlabel('Rate (bits)'); ylabel('Distorsion');
 title({'USQ'; 'Rate vs Distorsion'});
 
% soundsc(in,8000);
soundsc(OUT(3,:),8000);

% figure, plot(in,outq);
% title(['Uniform Scalar Quantizer. ' num2str(n_bits) ' bits, m= ' num2str(m) ', xmax= ' num2str(xmax) ]);
% xlabel('Input');ylabel('Output');
