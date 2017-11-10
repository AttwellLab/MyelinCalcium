function [Fig, Spk]=read_num(num2)

if isempty(num2)
    Fig=0;
    Spk=0;
else
    Fig=num2(1:end,1);
    Spk=num2(1:end,2);
end

end
