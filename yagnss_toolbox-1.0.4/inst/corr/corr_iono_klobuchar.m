function [diono] = corr_iono_klobuchar(X,Y,Z,Xs,Ys,Zs,alpha_ion,beta_ion,tGPS,Eph,freq);
%% function [diono] = corr_iono_klobuchar(X,Y,Z,Xs,Ys,Zs,alpha_ion,beta_ion,tGPS,Eph,freq);
%% Klobuchar model computation (source : NAVIPEDIA)
%%
%% Jacques BEILIN - 2008-11-22
%% Clement FONTAINE - 2013-11-18
%% 
%% Input : 
%% - X,Y,Z : station approx. coordinates (m)
%% - Xs, Ys, Zs : satellite coordinates (m)
%% - alpha_ion, beta_ion : klobuchar parameters (present in NAV_header)
%% - tGPS : mjd
%% - Eph : ephemeris (from get_ephemeris())
%% - used frequency ('F1' or 'F2')
%%
%% Output : 
%% - diono : ionospheric correction (m)
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if X == 0
	diono = 0;
	return;
end

% Klobuchar model computation for GPS-L1

c = 299792458.0;
% Geographic coordinates
[lon,lat,h] = tool_cartgeo_GRS80(X,Y,Z);
 
% Azimuth and Elevation computation
[A,E] = tool_az_ele_h(X,Y,Z,Xs,Ys,Zs);

lon = lon / pi;
lat = lat / pi;

A = A / pi;
E = E / pi;

% 1. Calculate the earth-centred angle   
psi = 0.0137/(E+0.11) - 0.022;

% 2. Compute the latitude of the Ionospheric Pierce Point 
phi_i = lat + psi * cos(A);
if phi_i <= -0.416 
	phi_i=-0.416;
end
if phi_i >= 0.416 
	phi_i=0.416;
end		

% 3. Compute the longitude of the IPP. 
lambda_i = lon + psi * sin(A) / cos(phi_i);

% 4. Find the geomagnetic latitude of the IPP.
phi_m = phi_i + 0.064 * cos(lambda_i - 1.617);

% 5. Find the local time at the IPP.
tGPS = tGPS - floor(tGPS);
t = 43200*lambda_i+tGPS*86400;

if t>=86400
	t = t - 86400;
elseif t < 0
	t = t + 86400;
end

% 6. Compute the amplitude of ionospheric delay
AMP = alpha_ion(1) + alpha_ion(2) * phi_m^1  + alpha_ion(3) * phi_m^2  + alpha_ion(4) * phi_m^3;
if AMP < 0 
	AMP = 0;
end

% 7. Compute the period of ionospheric delay.
PER = beta_ion(1) + beta_ion(2) * phi_m^1  + beta_ion(3) * phi_m^2  + beta_ion(4) * phi_m^3;
if PER < 72.000 
	PER = 72.000;
end	

% 8. Compute the phase of ionospheric delay.
x = (2 * pi *(t -50400))/PER;


% 9. Compute the slant factor 
F = 1.0 + 16.0 * (0.53-E)^3;

% 10. Compute the ionospheric time delay.
diono = F * 5.0*10^-9;
if abs(x) <=	pi/2 
	diono = diono + F * (1-(x^2)/2 + (x^4)/24) * AMP;
end

% 11. ionospheric time delay converted to distance
diono = c * diono;

% Transformation for others frequencies

if (isfield(Eph,'const') && isfield(Eph,'PRN'))
	const = Eph.const;
	PRN = Eph.PRN;
	
	if(strcmp(const,'G')) % GPS
	
		if(strcmp(freq,'F1'))
			% diono already set
			return;
		elseif(strcmp(freq,'F2'))
			diono = diono * ( 1575.42 / 1227.60 ) ^ 2;
			return;
		else
			tool_print_info(sprintf('Frequency %s is not defined : change it for F1 or F2. diono set to 0',freq),2);
			diono = 0;
			return;
		end
		
	elseif(strcmp(const,'E')) % Galileo
	
		if(strcmp(freq,'F1'))
			% diono already set
			return;
		elseif(strcmp(freq,'F2')) % E5a
			diono = diono * ( 1575.42 / 1176.45 ) ^ 2;
			return;
		else
			tool_print_info(sprintf('Frequency %s is not defined : change it for F1 or F2. diono set to 0',freq),2);
			diono = 0;
			return;
		end
	
	elseif(strcmp(const,'R')) % Glonass
	
		if(isfield(Eph,'freq_num'))
		
			if(strcmp(freq,'F1'))
				diono = diono * ( 1575.42 / (1602.00 + Eph.freq_num*0.5625) ) ^ 2;
				return;
			elseif(strcmp(freq,'F2')) % 
				diono = diono * ( 1575.42 / (1246.00 + Eph.freq_num*0.4375) ) ^ 2;
				return;
			else
				tool_print_info(sprintf('Frequency %s is not defined : change it for F1 or F2. diono set to 0',freq),2);
				diono = 0;
				return;
			end
		
		else
			tool_print_info('Glonass constellation with sp3 file : freq_num is missing, diono = 0',2);
			diono = 0;
			return;
		end
	
	else
		tool_print_info('Constellation not recognized : diono set to 0',2);
		diono = 0;
		return;
	end
		
	
	
else
	tool_print_info('Eph is empty : diono set to 0',2);
	diono = 0;
	return
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
