function value = get_param(params, fieldName, defaultValue)
%GET_PARAM Read a struct field with a default fallback.

    if isstruct(params) && isfield(params, fieldName) && ~isempty(params.(fieldName))
        value = params.(fieldName);
    else
        value = defaultValue;
    end
end
