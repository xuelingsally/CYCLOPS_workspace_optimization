function solve_all()

% problems ir0, ir3,  removed
problems = char('ack', 'ap', 'bf1', 'bf2', 'bhs', 'bl', 'bp', 'cb3',...
    'cb6', 'cm2', 'cm4', 'da', 'em_10', 'em_5', 'ep', 'exp', 'fls', 'fr',...
    'fx_10', 'fx_5', 'gp', 'grp', 'gw', 'h3', 'h6', 'hm', 'hm1', 'hm2',...
    'hm3', 'hm4', 'hm5', 'hsk', 'hv', 'ir1', 'ir2', 'ir4',...
    'ir5', 'kl', 'ks', 'lj1_38', 'lj1_75', 'lj1_98', 'lj2_38', 'lj2_75',...
    'lj2_98', 'lj3_38', 'lj3_75', 'lj3_98', 'lm1', 'lm2_10', 'lm2_5',...
    'lms1a', 'lms1b', 'lms2', 'lms3', 'lms5', 'lv8', 'mc', 'mcp', 'mgp',...
    'mgw_10', 'mgw_2', 'mgw_20', 'ml_10', 'ml_5', 'mr', 'mrp', 'ms1',...
    'ms2', 'nf2', 'nf3_10', 'nf3_15', 'nf3_20', 'nf3_25', 'nf3_30',...
    'osp_10', 'osp_20', 'plj_38', 'plj_75', 'plj_98', 'prd',...
    'ptm', 'pwq', 'rb', 'rg_10', 'rg_2', 's10', 's5', 's7', 'sal_10',...
    'sal_5', 'sbt', 'sf1', 'sf2', 'shv1', 'shv2', 'sin_10', 'sin_20',...
    'stg', 'st_17', 'st_9', 'swf', 'sz', 'szzs', 'wf', 'xor', 'zkv_10',...
    'zkv_2', 'zkv_20', 'zkv_5', 'zlk1', 'zlk2a', 'zlk2b', 'zlk3a',...
    'zlk3b', 'zlk3c', 'zlk4', 'zlk5', 'zzs');

problems = char('ap');

[nprobs,xx]=size(problems);

nl_path='.\models\';



for j=1:1;%8 % runs for rho
    for i=1:nprobs % for all problems
%        try
            disp(problems(i,:));
            [x,fx,nfo,deg,nit,npoll,spoll,nModels,RBFSuc]=RunPSwarm(strcat(nl_path,problems(i,:)),5*10^-3);
            fid=fopen(strcat('results_pswarm',num2str(j),'.tex'),'a');
            fprintf(fid,'%s & %d & %d & %d & %d & %d & %d & %d & %d\\\\\n', problems(i,:), nfo, fx, deg, nit, npoll, spoll,nModels,RBFSuc);
            fclose(fid);
%        catch
%            fid=fopen(strcat('results_pswarm',num2str(j),'.tex'),'a');
%            fprintf(fid,'%s & --- & --- & --- & --- & --- & --- & --- & ---\\\\\n', problems(i,:));
%            fclose(fid);
%        end
    end
end
