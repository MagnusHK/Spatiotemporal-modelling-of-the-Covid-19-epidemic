function I0 = FindICS(im, s, N)

I0 = zeros(size(im));

if s*N > 0.5
    %Finds a random number of points that starts with virus
    num = randi([2*ceil(s*N), sum(im, "all")]);
    %Determine the dosage in each point
    dose = s*N/num;
    counter = 0;

    %Find starting places
    while counter < num
        ind = RandomIndex(im);

        %Check if index is unique
        if I0(ind) == 0
            I0(ind) = dose;
            counter = counter + 1;
        end

    end
else
    ind = RandomIndex(im);
    I0(ind) = s*N;
end


end