function [Station] = get_station(RNX_header)

    Station.X = RNX_header.X;
    Station.Y = RNX_header.Y;
    Station.Z = RNX_header.Z;
    
    
   
    