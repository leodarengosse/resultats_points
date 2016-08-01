function [Point] = get_data_point(Points_data, Num)
   
    tab = Points_data;
    l = size(tab,1);
    
    for i=1:1:l+1
        if (i==l+1)
                disp('Unknown Point');
        else
            if (tab{i,1} == Num)
           Point = struct('Num', tab{i,1}, 'Time', tab{i,2}, ...
               'Acc', tab{i,6}, 'List_PRN', convertVect(tab{i,7}), ...
               'Lat',tab{i,3}, 'Long', tab{i,4}, 'Alt', tab{i,5});
            break
       
            else
                disp('Research...');
                Point = 0;
            end
        end
    end
end

%%
%Function who convert the String of PRN number into a vector whith each
%numbers (so we have a PRN for an indice given)

    function [Vec_PRN] = convertVect(List_PRN)
    
        str = strsplit(List_PRN, ','); %découpage de la chaine
    
        l = length(str);
    
        Vec_PRN = zeros(l,1); %initialisation du vecteur
    
        for i=1:1:length(str)
        
            Vec_PRN(i) = str2num(str{i});
    
        end
    end

            
          
            

 

    