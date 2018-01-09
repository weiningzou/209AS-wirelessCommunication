function ser = project_test_v2_1(SNR_dB,switch1)

K = 64;
M = 64;
CP = 16;
NoSym = 500;
NoBits = log2(M)*NoSym;
alpha = 4*pi*50*5*(10^-8);

%% Transmitter BITS
% b1 =  randi([0,1],NoBits,1);
% stem(b1);
% axis([0 60 0 1]);
B = randi([0,1],NoBits,K);
%% QAM
X = qammod(B,M,'UnitAveragePower',true,'InputType','bit');
%scatterplot(X(:,10));
%% PILOT IFFT AND PREFIX
for s = 1:NoSym
     samples = X(s,:);
    pilot=(1+1j)/abs(1+1j);
    for i = 8:16:K
        samples(i) = pilot;
    end
    samples_IFFT = ifft(samples,K);
    samples_CP = [samples_IFFT(K-CP+1:K) samples_IFFT];
    OFDMsymbols(s,:) = samples_CP;
    
end

OFDMSeries = reshape(OFDMsymbols,1,[]);

OFDMSeries = OFDMSeries./sqrt(mean(abs(OFDMSeries).^2));
%%  CHANNEL MULTIPATH
h = 1;
%h=[-0.7298-0.2884i  -0.1235+0.0411i  -0.0108-0.0381i  -0.0654-0.0165i];
%h = [-1.4896+0.7091i 0.4849+0.1602i 0.0339-0.0344i -0.2158+0.0961i -0.0675+0.0593i 0.1587-0.0241i -0.0049-0.0335i -0.0232+0.0434i];
channel_output= conv(OFDMSeries,h);
output_temp = reshape(channel_output,NoSym,[]);


for s = 1:NoSym
    w= sqrt(alpha).*randn(1,80)...
    +1i*sqrt(alpha).*randn(1,80);
    phi(1) = w(1);
    for n = 1:79
        phi(n+1) = phi(n)+w(n);  
    end
    a = output_temp(s,:);
    for i = 1:(K+CP)
        b(s,i) = a(i)*exp(phi(i)*1j);
    end
end
channel_real_output = reshape(b,1,[]);
%AWGN
%SNR_dB=50;
SNR_lin=10^(SNR_dB/10);
noisePower=sum(abs(h).^2)./SNR_lin;
wNoise= sqrt(0.5*noisePower).*randn(1,length(channel_output))...
    +1j*sqrt(0.5*noisePower).*randn(1,length(channel_output));
if switch1 == 1
    y = channel_real_output+ wNoise;
else
%y= channel_output.*exp(phi*1i) + wNoise;
    y = channel_output + wNoise;
end

    %y = channel_output+ wNoise;

%outputSymbol = reshape(channel_real_output,NoSym,[]);


%% Rx
%REMOVE CYCLIC PREFIX
samples_parallel= reshape(y(1:NoSym*(K+CP)),NoSym,[]);
samples_noCP= samples_parallel(:,CP+1:end);

%FFT
for s = 1:NoSym
    Y=fft(samples_noCP(s,:),K);
    for i=8:16:K
        channelEstimate(i)=pilot'*Y(i)./power(abs(pilot),2);
        channelEstimate(i-7:i+8)=channelEstimate(i);
        
    end
    for i=1:K
        X_hat(s,i) = channelEstimate(i)'*Y(i)...
        ./power(abs(channelEstimate(i)),2);
        %X_hat(s,i) = Y(i);
        
    end
end
%% correction
% 
% for k = 1: 64
%     for l = 1:64
%         E(k,l) = exp(-1*abs(k-l)*alpha/2);
%     end
% end
% %get Rjmjm
% Rj = fft2(E)/(K+CP)^2;
% %get Remem
% %Re = 0
% % for l = 1:(K+CP)
% %     Re = Re + sum(diag(Rj(:,l))) + alpha;
% % end
% Re = sum(diag(Rj)) + alpha;
% %MMSE
% Mm = Rj* X(1,:)'*(X(1,:) * Rj * X(1,:)' + Re)^-1;
% Jm = Mm * X_hat;   
% % for s = 1:NoSym
% %     for n = 1:(K+CP)
% %         for p = 1:(K+CP)
% %             R = fft2(E,80,80);
% %         end
% %     end
% % end
% X_linear = reshape(X_hat,1,[]);
% X_estimate = conv(X_linear,Jm);
%% Estiamted bits
B_hat=qamdemod(X_hat,M,'UnitAveragePower',true,'OutputType', 'bit');
%scatterplot(X_hat(:,30));
% h1 = scatterplot(X_hat(:,10));
% title('Subcarrier: 10');
% h2 = scatterplot(X_hat(:,25));
% title('Subcarrier: 25');
% h3 = scatterplot(X_hat(:,30));
% title('Subcarrier: 30');
% 
% h4 = figure;
% u1 = uipanel('position',[0,0,1/3,1]);
% u2 = uipanel('position',[1/3,0,1/3,1]);
% u3 = uipanel('position',[2/3,0,1/3,1]);
% 
% set(get(h1,'Children'),'parent',u1);
% set(get(h2,'Children'),'parent',u2);
% set(get(h3,'Children'),'parent',u3);
% close(h1,h2,h3);
%%SER vs SNR plot 
total_symbol = 500*64;
error_count = 0 ;
is_error = 0;
for sys = 1:500
       for sc = 1:64
            bt = 1;
            is_error=0;
            while (bt<=6)&&(is_error == 0) 
                if  (B(bt+(sys-1)*6 , sc) ~= B_hat(bt+(sys-1)*6 , sc))&&sc~=8 &&sc~=24 &&sc~=40 &&sc~=56 
                    
                    is_error = 1;
                    error_count = error_count + 1; 
                else
                    bt = bt + 1;
                end    
            end
       end     
end

ser = error_count/total_symbol;
