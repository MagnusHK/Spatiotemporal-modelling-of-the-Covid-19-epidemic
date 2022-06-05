function AU = AmultBoundary(U, int, m, n, h, D)

% D has to be a matrix

U = reshape(U, [m,n]);
AU = zeros(size(U));


for i = 1:m         %Loop for rows
    for j = 1:n     %Loop for coloumns


        if int(i,j)   %Checking points

            %Check if point to the left is not in domain
            if (j==1 || ~int(i,j-1)) && j<n
                AU(i,j) = AU(i,j) + 2*U(i,j+1);

            %Check if point to the right is not in domain
            elseif (j==n || ~int(i,j+1)) && j>1
                AU(i,j) = AU(i,j) + 2*U(i,j-1);
            
            else
                AU(i,j) = AU(i,j) + U(i,j-1) + U(i,j+1);

            end

            %Check if point above is not in domain
            if (i==1 || ~int(i-1,j)) && i<m
                AU(i,j) = AU(i,j) + 2*U(i+1,j);

            %Check if point below is not in domain
            elseif (i==m || ~int(i+1,j)) && i>1
                AU(i,j) = AU(i,j) + 2*U(i-1,j);

            else
                AU(i,j) = AU(i,j) + U(i-1,j) + U(i+1,j);
            end

            AU(i,j) = D*(AU(i,j) - 4*U(i,j));


        end

    end


end

AU = AU(:)/(h^2);

end