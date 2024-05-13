function [ enable ] = tx_counter( hit,impulse )

persistent counter first;

    if isempty(first)
        % initialisation
        counter = 0;
        first = 1;
    end
    counter = counter + 1;
    
    if counter == hit
        enable = impulse;
        counter = 0;
    else
        enable = false;
    end
end