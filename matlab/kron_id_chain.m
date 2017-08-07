function res = kron_id_chain( id_dims )
    %KRON_ID_CHAIN Creates a chain of tensor products of identity matrices of
    %specified dimension.

    res = eye(id_dims(1));

    for k = 2:length(id_dims)
        res = kron(res, eye(id_dims(k)));
    end
end