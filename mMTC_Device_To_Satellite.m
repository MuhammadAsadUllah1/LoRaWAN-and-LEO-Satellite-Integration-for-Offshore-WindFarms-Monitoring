%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
clc;close all;clear; % reset all
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%% Gains and Pt are converted into linear form
Gr=(10.^((22.6)/10));      %22.6 dBi
%https://lora.readthedocs.io/en/latest/
Gt=(10.^((2.15)/10));      %2.15 dBi 

Freq_Band = 868e6;         % 868 MHz (frequency band Europe)
SpeedLight  = 3e8;         % Speed of light

wavelength = SpeedLight/Freq_Band; %Wavelength
Channels = 8;              % Eight frequency channels
Pt = 10^(14/10)/1000;      % Transmit Power of LoRa 14 dBm
MonteCarlo = 1e2;          % 1e5 for the results in paper 
std = 0.1;                 % Shadowing

D_SNR = 10.^([-6,-15,-20]./10); % SF specific LoRa demodulator thresholds
do= 1000;                  %1km

Grid=csvread('Grid.mat');  %Coordinates of Wind turbines
Total_Wind_Turbines = length(Grid);

Sensors_per_Turbines = 50;
Total_Wind_Turbines = Total_Wind_Turbines * Sensors_per_Turbines;

R_Time = 600;              % Report period
Period = 60.*60.*1000;     % 1 hour in milliseconds
Transmissions = (Period/(R_Time*1000))*Total_Wind_Turbines; 

std = 0.1;                  % STD of shadowing in dB

D_SNR = 10.^([-6,-15,-20]./10); % SF specific LoRa demodulator thresholds
E = 10:1:90;                    % Elevation Angles
R = 6378e3;                     % Radius of earth
H = 780e3;                      % Orbital height

Satellite_subpoint=54;                    %Center of satellite beam at 54th Wind tirbine (in mid of farms)

[Satellite_Link_Farms,Ground_distance,Difference_from_GW_Slant_Range] = Simulations_Distance_Points(R,H,Satellite_subpoint);

[Elevation_Angles, Distance] = Satellite_Geometry(sort(Ground_distance),H);


Angles = [10 20 30 40 50 60 70 80 90];
K_Factor = [1.24 3.07 3.24 3.6 3.89 5.63 9.77 17.06 25.11];
K = sort(interp1(Angles,K_Factor,Elevation_Angles),'descend');
eta = 2;                   % path loss exponent for free space (LOS)

tic
[P_SNR] = Probability_SNR (Pt,Gt,Gr,K,D_SNR,Distance,MonteCarlo,wavelength,eta,std);

figure
plot(Distance/1000,P_SNR,'-','LineWidth',2)
ylabel('Probability ($P_{SNR}$)','Interpreter','Latex','FontSize', 12)
xlabel('Distance from node to satellite (km)','Interpreter','Latex','FontSize', 12)
axis([Distance(1)/1000 Distance(end)/1000 0 1])
legend('P_{SNR} SF7','P_{SNR} SF10','P_{SNR} SF12')

toc

%% PSIR: Considering interfering signals 
ToA = [61.696 370.688 1482.752]; %Time on Air for 10 bytes in milliseconds
P_SIR=zeros(length(D_SNR),length(Distance));

% Count : runs for SF7, SF10, SF12
% c : it picks the a ground distance (increases with step)
Count_capture=zeros(length(Distance),MonteCarlo);
for count=1:3
    
    for c=1: length(Distance)   
        
        discarded = 0;     % To calculate destructive collision at a given distance(as function of elevation angle)

        for m=1:MonteCarlo
        simultaneous=0;
        %% calculating the number of interfering signals
        
        Transmission_per_channel=round(Transmissions/Channels);  % Transmission per channel

    
        clear TX_stamp
        TX_stamp  = rand(1,round(Transmission_per_channel))*3600;   %Generating transmissions
        %TX_stamp  = rand(1,round(Total_Wind_Turbines/8)*6)*3600;
        
        Tsort=sort(TX_stamp);
        TimeStamp=Tsort;

        TimeStampEnd = TimeStamp + (ToA(count)/1000);       %End time of transmission of devices
    
        Transmit = randi(length(TimeStamp));                % Selecting one random device
        % Ts: transmission started
        Ts= TimeStamp(Transmit);
        % Tend: transmission ended
        Tend= Ts + (ToA(count)/1000);
        
        % Finding simultaneous transmission (collisions) in that time window 
        
        Starting = find(TimeStamp >= Ts & TimeStamp <= Tend);
        
        %Starting: Previous transmissions in that time window?
        Ending = find(TimeStampEnd >= Ts & TimeStampEnd <= Tend);
        %unique : find nodes other than the desired one in that window 
        
        simultaneous = length(unique([Starting Ending]))-1;          % colliding signals

                if simultaneous > 0                                  % if there is any simultaneous transmission (if collision)
                  %% Generating location of interfering signals
                  
                  index = randi(length(Grid),simultaneous,1);
                  
                  Interfering_Nodes_X = Satellite_Link_Farms(c,index,1);
                  Interfering_Nodes_Y = Satellite_Link_Farms(c,index,2);
                  
                  Interfering_Nodes =[Interfering_Nodes_X', Interfering_Nodes_Y'];
                  Location_Nodes_Int = sqrt((Interfering_Nodes(:,1)-Grid(Satellite_subpoint,1)).^2 + (Interfering_Nodes(:,2)-Grid(Satellite_subpoint,2)).^2)';
                  
   
                clear dPropogation
                clear E_dpro
                clear kC
                clear E_AngPro
                dPropogation=zeros(1,simultaneous);
                E_dpro=zeros(1,simultaneous);
                
                for track=1:1:length(Location_Nodes_Int)
                dPropogation(1,track) = sqrt(H^2 + Location_Nodes_Int(track).^2);
                end
                
                for Fin_ang=1:length(Location_Nodes_Int)
                E_dpro(Fin_ang) = (H*((H+2*R)) - dPropogation(Fin_ang).^2)./(2.*dPropogation(Fin_ang).*R);
                end

                E_AngPro=asind(E_dpro);
                
                kC = interp1(Angles,K_Factor,E_AngPro);
                
                %% Rician fading for interfering signals
               
                muC = sqrt(kC./(2.*(kC+1)));    % Mean 
                sigmaC = sqrt(1./(2.*(kC+1)));  % Variance 

                hrC=sigmaC.*randn(1,simultaneous)+muC;
                hiC=1j.*(sigmaC.*randn(1,simultaneous)+muC);
                
                h1C=(abs(hrC+hiC)).^2;
                
                
                %% Log Distance Path loss
                Lo = eta * 10*log10(4*pi*do*(1/wavelength));
                clear L_int
                clear L_Lin_int
                clear pr_h_g_I
                L_int=zeros(1,simultaneous);
                L_Lin_int=zeros(1,simultaneous);
                
                L_int  = Lo + 10*eta*log10(dPropogation/do) + std*randn(1,1);
                L_Lin_int=10.^((L_int)./10);
                
                %% total interferance = sum of all the interfering signals 
                pr_h_g_I = sum(Pt.*h1C.*Gt.*Gr.*(1./L_Lin_int));

                %% Rician fading for desired signals                
                kD = K(c);
                muD = sqrt(kD./(2*(kD+1)));     % Mean 
                sigmaD = sqrt(1./(2*(kD+1)));   % Variance 
                
                hrD=sigmaD*randn(1,1)+muD;
                hiD=1j.*(sigmaD*randn(1,1)+muD);
                
                h1D=(abs(hrD+hiD)).^2;
                
                mean(h1D);

                %% Received power of desired signal
                L_des  = Lo + 10*eta*log10(Distance(c)/do) + std*randn(1,1);
                L_Lin_des=10.^((L_des)./10);
                
                pr_h_g_D = Pt.*h1D.*Gt.*Gr.*(1./L_Lin_des);

                %% Data Collisions and Spreading Factor Orthogonality
                
                 %"Two packets with the same spreading factor arriving at the same time... 
                 ...on the same channel might result in a collision. 
                  ...However if one of the two packets is stronger by six dB, it will survive."
                  
                % https://lora-developers.semtech.com/library/tech-papers-and-guides/lora-and-lorawan/
                % 6 dB  = 3.9811 (approx 4)
                
                      if  pr_h_g_D  < (pr_h_g_I*4)
                      discarded = discarded + 1; %% Destructive collision
                      end
                     
               end
        end
     P_SIR(count, c) = 1- (discarded/MonteCarlo); %% Success in presence of interfering nodes
    end
end
%% Plotting PSIR: Considers interference 
figure
plot(Distance/1000,P_SIR(1,:),'-r','LineWidth',2);
hold on
plot(Distance/1000,P_SIR(2,:),'-b','LineWidth',2);
plot(Distance/1000,P_SIR(3,:),'-k','LineWidth',2);
grid on
ylabel('Probability ($P_{SIR}$)','Interpreter','Latex','FontSize', 12)
xlabel('Distance from user to satellite (km)','Interpreter','Latex','FontSize', 12)
axis([Distance(1)/1000 Distance(end)/1000 0 1])
legend('PSIR SF7','PSIR SF10','PSIR SF12')
toc
