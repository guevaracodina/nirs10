function a = nirs_rep_array(a,n)
a = str2num(a);
if ~isempty(a)
    if n > 1 && length(a) == 1
        a = repmat(a,[1 n]);
    end
end