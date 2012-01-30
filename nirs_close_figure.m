function nirs_close_figure(h)
try 
    close(h)
catch exception
    disp(exception.identifier);
    disp(exception.stack(1));
    disp('Problem closing figure -- perhaps user closed it!');
end