function s = synthesis2(E, ZC, V, A, P, ulen)
% We have included all the estimated parameters as input arguments
% but here we only use A!
n_frames = size(A,1); % Assuming filter coefficients are stored row-wise
% Create a pulse train excitation:
cp = 800; % Constant pitch period in samples
pexc = zeros(n_frames*ulen, 1);
pexc(1:cp:end) = 1;
% Create noise excitation:
nexc = 2*rand(n_frames*ulen,1) - 1;
n1 = 1;
n2 = ulen;
Z = [];
s = zeros(n_frames*ulen,1);
for n=1:n_frames
    % Filter the excitation through the production (vocal tract) filter:
    [s(n1:n2), Z] = varifilter(1, A(n,:), nexc(n1:n2), Z);
    
    Es(n) = 1/(n2-n1+1)*sum(s(n1:n2).^2);
    %To normalize we resolve sum((s(n1:n2)*corr).^2) = E(n)
    
    corr = sqrt(E(n)/Es(n));
    s(n1:n2) = s(n1:n2)*corr;
    
    Es(n) = 1/(n2-n1+1)*sum(s(n1:n2).^2);

    n1 = n1+ulen;
    n2 = n2+ulen;
    
end
figure, plot(E/max(E)); hold on; plot(Es/max(Es));

end
