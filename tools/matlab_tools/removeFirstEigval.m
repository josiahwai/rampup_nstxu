function Mhat = removeFirstEigval(M)

[V,D] = eig(M);
eig_max = max(real(diag((D))));

if eig_max > 0
    [row, col] = find(real(D) == eig_max);
    D(row,col) = 0;    
end
Mhat = V * D / V;

end
