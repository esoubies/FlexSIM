function DispMsg(verb,msg,level)

if nargin <3, level=1; end

if verb>0
    if level
        disp(msg);
    else
        fprintf('%s',msg);
    end
end

end