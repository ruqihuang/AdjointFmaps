function [present, mappedverts] = ismapped(map)
    nv = length(map);
    mappedverts = unique(map);
    present = zeros(nv,1);
    present(mappedverts) = 1;
end