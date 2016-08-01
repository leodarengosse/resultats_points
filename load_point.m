function [Points_data] = load_point(filename)
    
    [Num, Time, Lat, Long, Alt, Acc, List_PRN] = textread(filename, '%s %f %f %f %f %f %s');    
     
    l = length(Num);
    
    Points_data = cell(l,7);
    
    List_PRN{1}(1); 
    Temp_PRN = zeros(l,1);
    size(List_PRN{1});

    for i=1:1:l
        Points_data{i,1} = Num{i};
        Points_data{i,2} = Time(i)/1000; %%conversion en secondes
        Points_data{i,3} = Lat(i);
        Points_data{i,4} = Long(i);
        Points_data{i,5} = Alt(i);
        Points_data{i,6} = Acc(i);
        Points_data{i,7} = List_PRN{i};
        
    end
end




