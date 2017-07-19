function [Success,Problem,Population]=ResamplingStep2(Problem, Population, Success, varargin)

% lb = Problem.LB;
% ub = Problem.UB;
WA = [];
CWA_index = [];

x = Population.x;

for i=1:Population.Size     
    if Population.fy(i) < 1000
        % Create Weight Array Based on Feasible Particles
        WA(end+1) = 1.6 - Population.fy(i);
        CWA_index(end+1) = i;
    end
end

% Normalise the weight array
WA = WA / sum(WA);

% Cummulative Weight Array - for sampling
CWA = zeros(size(WA,1), 1);
CWA(1,1) = WA(1,1);
for i=2:size(WA,2)
    CWA(1,i) = CWA(1,i-1) + WA(1,i);
end

chosen_particle = 1;

for i=1:Population.Size
    rand_no = unifrnd(0,1);
    for j=1:size(CWA,2)
        if rand_no <= CWA(1,j)
            chosen_particle = CWA_index(1,j);
            break;
        end
    end
%     Metro_hastings = (1.1-Population.fy(chosen_particle)) * Population.fy(chosen_particle)/((1.1-Population.fy(i)) * Population.fy(i));
%     rand_no = unifrnd(0,1);
%     if Metro_hastings >= rand_no


%     likelihood = 1.5 - Population.fy(i);
%     rand_no = unifrnd(0,1);
%     
%     if rand_no >= likelihood
        
    Population.x(i,:) =  x(chosen_particle,:);
    Population.vx(i,:) = zeros(1,Problem.Variables);

    [Problem,ObjValue]=...
        PenaltyEval(Problem, Population.x(i,:), 0);
    % Was progress attained for the current particle?
    if Population.fy(i)>ObjValue
        % Yes. Update best particle position
        Population.fy(i)=ObjValue;
        Population.y(i,:)=Population.x(i,:);

        % Check if new leader is available
        if Population.fy(Population.Leader)>Population.fy(i)
            Population.Leader=i;
        end
    end
%     end
end
Problem.ResamplingCount = Problem.ResamplingCount + 1;

end
