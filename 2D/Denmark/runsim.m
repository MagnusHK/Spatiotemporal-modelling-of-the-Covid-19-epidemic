function y0 = runsim(model, tspan, y0, im, beta, p, D, m, n, h, path, MakeRandom)
inter = 1;
tmax = tspan(2);
tmin = tspan(1);
data = [];
time_init = tmin;
handled = false;

totaliterations = (tmax - tmin)/inter;
curr_iteration = 0;

percentage = 0.3;       % 1% infected in a point found

for i = tmin+inter:inter:tmax
    curr_iteration = curr_iteration+1;

    tspan = [i-inter, i];

    [t, y] = ode45(model, tspan, y0, [], im, beta, p, D, m, n, h);

    % Randomized infection

    y0 = y(end,:)';

    if (randi(5) == 1) && MakeRandom
        ind = RandomIndex(im);
        
        y0(m*n+ind) = y0(m*n+ind) + percentage*y0(ind);   %Newly infected people
        y0(ind) = (1-percentage)*y0(ind);   %People removed from susceptible
    
        fprintf('Added %.3f infected at position %d\n', percentage*y0(ind), ind);
        
    end

    if i == tmax
        data = [data; t(1:end,1), y(1:end,:)];
    else
        data = [data; t(1:end-1,1), y(1:end-1,:)];
    end
    
    %Memory handling
    Gigbytes = getfield(whos('data'),'bytes')*1e-9;

    if Gigbytes >= 5
        
        name = ['Denmark_Time_', num2str(time_init), '_To_', num2str(tspan(2))];

        fprintf('Saving data under the name: %s\n', name);

        save([path, name, '.mat'], 'data', '-v7.3');
        time_init = tspan(2);

        clear data;
        data = [];
        
        if tspan(2) == tmax
            handled = true;
        end
    end

    fprintf('Iteration %d of %d \n', curr_iteration, totaliterations);
end

if ~handled
    name = ['Denmark_Time_', num2str(time_init), '_To_', num2str(tspan(2))];
    
    save([path, name, '.mat'], 'data', '-v7.3');
end


end