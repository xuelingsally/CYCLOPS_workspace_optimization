function [Success,Problem,Population]=ResamplingStep(Problem, Population, Success, varargin)

lb = Problem.LB;
ub = Problem.UB;

for i=1:Population.Size
    likelihood = 1.5 - Population.fy(i);
    rand_no = unifrnd(0,1);
    
    if rand_no >= likelihood
        Problem.ResamplingCount = Problem.ResamplingCount + 1;
        Population.x(i,:)=unifrnd(lb,ub);
        Population.vx(i,:)=zeros(1,Problem.Variables);
        
        [Problem,ObjValue]=...
            PenaltyEval(Problem, Population.x(i,:), 0);
        % Was progress attained for the current particle?
        if Population.fy(i)>ObjValue
            % Yes. Update best particle position
            Population.fy(i)=ObjValue;
            Population.y(i,:)=Population.x(i,:);

            % Check if new leader is available
            if Population.fy(Population.Leader)>Population.fy(i) || Population.Leader==i
                Population.Leader=i;
                % Particle swarm iteration declared as successful
                %Success=true;
                % Reset last success direction for pattern search
                %Problem.Poll.LastSuccess=[];
            end
        end
    end

end


end
