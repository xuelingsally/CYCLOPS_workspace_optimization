function [Problem,Population]=Resample_Particles(Problem, Population, varargin)

% First Ensure that the inactive particles do no contribute to weights
for i=1:Population.Size
    if ~(Population.Active(i))
        Population.weight(i) = 0;
    end
end

% Normalising the weights
Population.weight = Population.weight/sum(Population.weight);

% Creating Cumulative Weight Array (distribution) for Sampling from
CWA = zeros(Population.Size, 1);
CWA(1,1) = Population.weight(1);
for i=2:Population.Size
    CWA(i,1) = Population.weight(i) + CWA(i-1, 1);
end

x = Population.x;
y = Population.y;
fy = Population.fy;

Population.Leader = 1;

% Sampling the particles
chosen_particle = 1;
for i=1:Population.Size
    rand_no = unifrnd(0,1);
    for j=1:size(CWA,1)
        if rand_no <= CWA(j,1)
            chosen_particle = j;
            break;
        end
    end
    
    Population.x(i,:) =  x(chosen_particle,:);
    Population.vx(i,:) = zeros(1,Problem.Variables);
    Population.y(i,:) = y(chosen_particle,:);
    Population.fy(i) = fy(chosen_particle);
    Population.weight(i,1) = 1/Population.Size;
    %Population.Active(i) = 1;
    
    if Population.fy(Population.Leader)>Population.fy(i) || Population.Leader==i
        Population.Leader=i;
    end
    
end



end