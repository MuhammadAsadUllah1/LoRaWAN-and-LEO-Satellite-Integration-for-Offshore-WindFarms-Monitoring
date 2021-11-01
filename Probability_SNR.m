function [P_SNR] = Probability_SNR (Pt,Gt,Gr,k,D_SNR,Distance,MonteCarlo,wavelength,eta, std)

%% LOS : rician fading channel

%R = (sigma*G + mu) + j*(sigma*G + mu);

%A. Gaysin, V. Fadeev and M. Hennhöfer, "Survey of modulation and coding schemes for application in CubeSat systems," 
%2017 Systems of Signal Synchronization, Generating and Processing in Telecommunications (SINKHROINFO), Kazan, 2017,

%https://ieeexplore.ieee.org/abstract/document/7997514

%k                         % Rician factor K 
% mu = sqrt(k/2*(k+1));    % Mean 
% sigma = sqrt(1/2*(k+1)); % Variance 

% hr=sigma*randn(1,MonteCarlo)+mu;
% hi=1j*sigma*randn(1,MonteCarlo)+mu;

% h1=(abs(hr+hi)).^2;

%% PSNR
Noise_Figure=6;           
Band_Width=125e3; %Band Width 

NoisePower=(10.^((-174+Noise_Figure+10*log10(Band_Width))/10))/1000;

P_SNR = zeros(length(D_SNR),length(Distance));

do= 1000;   %1 km

for i=1:3 % SF 7, SF10 and SF12
    
    for pointer=1:length(Distance) % Distance of satellite from the ground user
       
        mu = sqrt(k(pointer)./(2.*(k(pointer)+1)));    % Mean 
        sigma = sqrt(1./(2.*(k(pointer)+1)));          % Variance 

        hr=sigma*randn(1,MonteCarlo)+mu;
        hi=1j*(sigma*randn(1,MonteCarlo)+mu);
        
        %h1 = (sqrt(0.5)*abs((sigma*randn(length(d1),samples) + mu) + (sigma*1j*randn(length(d1),samples))+mu).^2);

        h1=(abs(hr+hi)).^2;
        
        %% Log Distance Path loss
        Lo = eta * 10*log10(4*pi*do*(1/wavelength));
        L(pointer,:)  = Lo + 10*eta*log10(Distance(pointer)/do) + std*randn(1,MonteCarlo);
        
        %% Pr = pt*Gt*Gr*path loss * h  
        L_Lin(pointer)=10.^((L(pointer))./10);
        
        pr = Pt.*Gt.*Gr.*h1.*(1./L_Lin(pointer));


        %% P_SNR = Probability(SNR > D_SNR)
        
        P_SNR(i,pointer) = (sum((pr./NoisePower) >= D_SNR(i)))/MonteCarlo; 
    end
end
end