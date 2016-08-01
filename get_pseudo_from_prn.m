function Pseudo_range = get_pseudo_from_prn(RNX_header,RNX_data,PRN,epoch)

Pseudo_range = {};
for i=1:1:length(List_PRN)
    [Obs] = get_obs(RNX_header,RNX_data,'G',PRN,epoch);
    
    if (isempty(fieldnames(Obs)) || Obs.C1 == 0 )
        sprintf('PRN %s introuvable',num2str(List_PRN(i)))
    else
        Pseudo_range{i-1,1} = List_PRN(i);
        Pseudo_range{i-1,2} = Obs.C1;
    end
end