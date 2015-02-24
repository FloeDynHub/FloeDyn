% Convert, recursively, an object to a structure
function str = object2struct_rec( obj )
    if length(obj) > 1
        if isnumeric(obj) || ischar(obj) || islogical(obj)
            str = obj;
        else
            cnt = prod(size(obj));
            str(cnt) = object2struct_rec(obj(cnt));
            str = reshape(str, size(obj));
            for i = 1:(cnt-1)
                str(i) = object2struct_rec( obj(i) );
            end
        end
    else
        if isnumeric(obj) || ischar(obj) || islogical(obj)
            str = obj;
        elseif      iscell(obj)
            str = { object2struct_rec(obj{1}) };
        elseif  isobject(obj)
            str = object2struct_rec( struct(obj) );
        elseif  isstruct(obj)
            str = struct();
            props = fieldnames(obj);
            for i = 1:length(props)
                str = setfield(str, props{i}, object2struct_rec(getfield(obj, props{i})));
            end
        else
            str = obj;
        end;
    end
end