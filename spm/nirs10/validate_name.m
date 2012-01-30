function name = validate_name(name)
%check if there are unallowed characters in file name part and remove them
name = regexprep(name,'>','g');
name = regexprep(name,'<','l');
name = regexprep(name,' ','_');
name = regexprep(name,'?','u');
name = regexprep(name,'!','e');
name = regexprep(name,'|','v');
name = regexprep(name,'#','n');
name = regexprep(name,'%','p');
name = regexprep(name,'&','a');
name = regexprep(name,'+','s');
name = regexprep(name,'@','t');
name = regexprep(name,'$','d');
name = regexprep(name,'^','c');
name = regexprep(name,'(','_');
name = regexprep(name,')','_');