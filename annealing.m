function [an,path] = annealing(an,NumCity,dist_mat,path)
an.minE = inf;
an.maxE = 0;
while an.TrialN < an.MaxTrialN && an.AcceptN < an.MaxAcceptN
    % MOVE: RouletteWheelSelection:
    new_path = move_sets(NumCity,path);
    
    new_dist_energy = sum(dist_mat( ...
        (new_path-1)*NumCity+[new_path(2:NumCity) new_path(1)]));
    
    % generate a random number and accept if probability is greater
    % than the random number
    prob_accept = exp((an.dist_energy - new_dist_energy)/an.temp);
    if rand < prob_accept	% accept the dist_energy and path!
        an.dist_energy = new_dist_energy;
        path = new_path;
        an.minE = min(an.minE, an.dist_energy);
        an.maxE = max(an.maxE, an.dist_energy);
        an.AcceptN = an.AcceptN + 1;
    end
    an.TrialN = an.TrialN + 1;
end
end