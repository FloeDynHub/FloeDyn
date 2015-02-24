% Convert, recursively, a structure to an object
function obj = struct2object_rec( ref_obj, str )
    obj = ref_obj;
    if length(str) > 1
        if isnumeric(ref_obj) || ischar(ref_obj) || islogical(ref_obj)
            obj = str;
        else
            cnt = prod(size(str));
            obj(cnt) = struct2object_rec(ref_obj(1), str(cnt));
            obj = reshape(obj, size(str));
            for i = 1:cnt
                obj(i) = struct2object_rec(ref_obj(i), str(i) );
            end
        end
    else
        if isnumeric(ref_obj) || isnumeric(ref_obj) || islogical(ref_obj)
            obj = str;
        elseif iscell(ref_obj)
            obj = { struct2object_rec( ref_obj{1}, str{1} ) };
        elseif isobject(ref_obj) || isstruct(ref_obj)
            props = fieldnames(str);
            for i = 1:length(props)
                obj = setfield(obj, props{i}, ...
                        struct2object_rec( ...
                            getfield(ref_obj, props{i}), ...
                            getfield(str, props{i}) ...
                        ) ...
                );
            end
        else
            obj = str;
        end
    end
            
end