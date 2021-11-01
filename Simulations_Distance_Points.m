function [Satellite_Link_Farms,Distance_from_sub_point,Difference_from_GW_Slant_Range] = Simulations_Distance_Points(Earth_Radius,Orbit_Height,Subpoint);

Elevation_Simulation = [90:-1:10.01];
Slant_Range=zeros(1,length(Elevation_Simulation));
Base =zeros(1,length(Elevation_Simulation));
    
Location_X_Y=csvread('Grid.mat');
    
for WhichAngle=1:length(Elevation_Simulation)

    % Calculating the equation (5) of the paper    
    Eq_pt_1 = cosd(Elevation_Simulation(WhichAngle)).*cosd(Elevation_Simulation(WhichAngle));  
    Eq_pt_2 = ((Orbit_Height + Earth_Radius)./Earth_Radius)^2;
    
    Slant_Range(WhichAngle) = Earth_Radius.*(sqrt(Eq_pt_2- Eq_pt_1) - sind(Elevation_Simulation(WhichAngle)));  % Slant_Range(E): distance from user to satellite
        
        if(Elevation_Simulation(WhichAngle)==90)
            Base (WhichAngle) = abs(sqrt(Slant_Range(WhichAngle)^2-Orbit_Height^2));
        else
            Base (WhichAngle) = sqrt(Slant_Range(WhichAngle)^2-Orbit_Height^2);
        end
    
%% simulation points (distance): Satellite mobility
% Here we are increasing the distance from Wind turbine to the Satellite as
% shown on X-axis of Figure 8-10
    
    Satellite_Link_Farms(WhichAngle,:,:)=[(Location_X_Y(:,1)+ Base(WhichAngle)), Location_X_Y(: ,2)];
    Distance_from_sub_point(WhichAngle) = sqrt((Location_X_Y(Subpoint,1)-Satellite_Link_Farms(WhichAngle,Subpoint,1))^2 + (Location_X_Y(Subpoint,2) - Satellite_Link_Farms(WhichAngle,Subpoint,2))^2);
    
    Distance_from_sub_point_4_all(WhichAngle,:) = sqrt((Location_X_Y(Subpoint,1)-Satellite_Link_Farms(WhichAngle,:,1)).^2 + (Location_X_Y(Subpoint,2) - Satellite_Link_Farms(WhichAngle,:,2)).^2);
    
    Slant_Range_4_All(WhichAngle,:) =  sqrt(Orbit_Height^2 + Distance_from_sub_point_4_all(WhichAngle,:).^2);
    
    Difference_from_GW_Slant_Range(WhichAngle,:) = abs(Slant_Range_4_All(WhichAngle,Subpoint) - Slant_Range_4_All(WhichAngle,:));
end

end