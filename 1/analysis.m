function [E, ZC, V, A, P] = analysis(x, alen, ulen,voicedthreshold,M)
%function [E, ZC, V, A, P] = analysis(x, alen, ulen,M)

% Initialization
xlen = length(x);
naf = ceil(xlen/alen);
ZC = zeros(naf, 1);
V = zeros(naf, 1);
E = zeros(naf, 1);
A = zeros(naf, M+1);
P = zeros(naf,1);
% A(:,1) = ones(naf,1);


for n = 1:naf
    iframe = (n-1)*alen+1 :min(xlen,n*alen);
    xf = x(iframe);
    flen = length(iframe);
    % Inside loop
    E(n) = 1/alen*sum(x( iframe ).^2);
    
    %ZERO counting in a frame
    for i =iframe(1):iframe(end)
        if(x(i)*x(min(i+1,xlen)) <=0 ) ZC(n) = ZC(n) +1; end
    end
    ZC(n) = ZC(n)/flen;
    V(n) = ZC(n)<=voicedthreshold;%Decision making
    
    %LP
    c = xcorr(xf,xf,M);%, M);
    [a, e]= levinson(c(M+1:2*M+1)); % a is a vector always starting
    %with a 1.
    A(n,:) = a;
    %Find the pitch
    %DFT 
%     X = abs(fft(x(iframe)));
    corrframe = xcorr(x(iframe),x(iframe));
    
    [maxcorr] = max(corrframe);
    thresh = .1*maxcorr;
    for i = 1:length(corrframe)
        if corrframe(i)>=thresh
            index = i;
            break;
        end
    end
    P(n) = index;
    
%     i=1;
%     dX = diff(X);
%     dX=dX(1:round(end/2));
%     for i =1:length(dX)
%         if(dX(i)*dX(min(i+1,length(dX))) <=0 ) P(n) = i; break; end
%     end
%     
%     k = find(~dX);
%     P(n) = X(k);
end



end