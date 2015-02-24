% Convert, recursively, a structure to an object
function obj = struct2object( ref_obj, str )
    obj = struct2object_rec( ref_obj, str );
end