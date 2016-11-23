function x = solve_system(A,b)
    [U,S,V] = svd(A);
    [n,N] = size(A);

    sT = sum(diag(S));

    % size(V)

    x = zeros(N,1);
    for i = 1:N
        if (S(i,i))/sT>0.1
            x = x+(1/(S(i,i)))*V(:,i)*(U(:,i)'*b);
        else
            break;
        end
    end
end