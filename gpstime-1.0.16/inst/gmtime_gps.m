function TM_STRUCT = gmtime_gps(T)

% Test octave/matlab
if exist('OCTAVE_VERSION') % octave
    
    TM_STRUCT = gmtime(T);

else % matlab
    
    T = T/86400+719529; % sec to day , 719529 = day number between 1900 and 1970
    
    Tvec = datevec(T);
    TM_STRUCT.year = Tvec(1) - 1900;
    TM_STRUCT.mon = Tvec(2)-1;
    
    yearbegin = datenum(Tvec(1),1,1,0,0,0);
    
    TM_STRUCT.mday = Tvec(3);
    TM_STRUCT.yday = floor(T - yearbegin);
    TM_STRUCT.wday = weekday(T)-1;
    
	TM_STRUCT.hour = Tvec(4);
	TM_STRUCT.min = Tvec(5);
	TM_STRUCT.sec = fix(Tvec(6));
	TM_STRUCT.usec = (Tvec(6) - TM_STRUCT.sec)*1e6;
    
end
   

end
