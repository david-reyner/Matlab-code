function t = tangent_vector(A)
    % A is 2x3 matrix
    t = zeros(length(A),1);
    for k = 1:length(t)
        aux = A;
        aux(:,k) = [];
        t(k) = (-1)^(k-1)*det(aux);
    end
    t = t/norm(t);
end

