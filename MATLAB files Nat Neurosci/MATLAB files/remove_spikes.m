function [Reas]=remove_spikes(txt2)

if isempty(txt2)
    disp('No sheaths were removed')
else
Reas=txt2(2:end,4);
end

end
