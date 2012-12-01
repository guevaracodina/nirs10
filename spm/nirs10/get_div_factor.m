function [div_factor div_factor0] = get_div_factor(d,fs,bl_m,f,win,div_factor0)
try
    % Normalization factor
    switch bl_m
        case 0 %0: Median;
            div_factor = median(d(:,win),2);
        case 1 %1: Initial value - not too sensible - better last value
            div_factor = d(:,win(end));
        case 2 %mean
            div_factor = mean(d(:,win),2);
        case 3 %median - session 1
            %do nothing
            if f == 1
                div_factor0 = median(d,2);
            end
            div_factor = div_factor0;
        case 4 %mean of first 60 seconds
            div_factor = mean(d(:,win(1:(fs*60))),2);
        case 5 %mean of first 10 seconds
            div_factor = mean(d(:,win(1:(fs*10))),2);
        case 6 %mean of first 60 seconds of first session only
            if f == 1
                div_factor0 = mean(d(:,win(1:(fs*60))),2);
            end
            div_factor = div_factor0;
        otherwise %take median
            div_factor = median(d(:,win),2);
    end   
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Could not compute div_factor');
    div_factor = 1;
    div_factor0 = 1;
end