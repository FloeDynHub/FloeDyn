% Convert, recursively, an object to a structure 
function str = object2struct( obj )
    warning('off', 'MATLAB:structOnObject');
    str = object2struct_rec(obj);
    warning('on', 'MATLAB:structOnObject');
    
end