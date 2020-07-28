
function writeClusterFile
%% simulationBiolog
% kmax = 116;
% k = 1:20:kmax;
% for i = 1:length(k)
%     subfileName = ['sub',num2str(k(i)),'.sh'];
%     fptr = fopen(subfileName,'w');
%     fprintf(fptr,'#!/bin/bash\n');
% %     fprintf(fptr,'#SBATCH -A SNIC2020-7-19\n');
%     fprintf(fptr,'#SBATCH -A C3SE2020-1-8\n');
%     fprintf(fptr,'#SBATCH -N 1\n');
%     fprintf(fptr,'#SBATCH -n 20\n');
%     fprintf(fptr,'#SBATCH -o out.txt\n');
%     fprintf(fptr,'#SBATCH -t 50:00:00\n');
%     fprintf(fptr,'#SBATCH --mail-user=cheyu@chalmers.se\n');
%     fprintf(fptr,'#SBATCH --mail-type=end\n');
%     fprintf(fptr,'module load MATLAB intel/2018b GMP\n');
%     fprintf(fptr,['i=',num2str(k(i)),'\n']);
%     if i ~= length(k)
%         for m = 1:20
%             fprintf(fptr,['a',num2str(m),'=$i+' num2str(1*(m-1)) '\n']);
%             fprintf(fptr,['b',num2str(m),'=$i+' num2str(1*m-1) '\n']);
%         end
%         for m = 1:20
%             fprintf(fptr,['matlab -nodesktop -singleCompThread -r "simulationBiolog($a',num2str(m),',$b',num2str(m),')" &\n']);
%         end
%     else
%         for m = 1:kmax-k(i)+1
%             fprintf(fptr,['a',num2str(m),'=$i+' num2str(1*(m-1)) '\n']);
%             fprintf(fptr,['b',num2str(m),'=$i+' num2str(1*m-1) '\n']);
%         end
%         for m = 1:kmax-k(i)+1
%             fprintf(fptr,['matlab -nodesktop -singleCompThread -r "simulationBiolog($a',num2str(m),',$b',num2str(m),')" &\n']);
%         end
%     end
%     fprintf(fptr,'wait;\n');
%     fclose(fptr);
% end

%% simulationLowerUptake
ion_mw_list = [39 24 56 65 40 55 64 23];
for i = 1:length(ion_mw_list)
    subfileName = ['sub',num2str(ion_mw_list(i)),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    fprintf(fptr,'#SBATCH -A SNIC2020-7-19\n');
%     fprintf(fptr,'#SBATCH -A C3SE2020-1-8\n');
    fprintf(fptr,'#SBATCH -N 1\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH -t 50:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=cheyu@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'module load MATLAB intel/2018b GMP\n');
    fprintf(fptr,['matlab -nodesktop -singleCompThread -r "simulationLowerUptake(',num2str(ion_mw_list(i)),')" \n']);
    fclose(fptr);
end

%% simulationCNPS
% [~,raw,~] = xlsread('Exp_biolog.xlsx');
% metNameID = raw(2:end,2);
% 
% kmax = sum(~ismember(metNameID,''));
% k = 1:20:kmax;
% for i = 1:length(k)
%     subfileName = ['sub',num2str(k(i)),'.sh'];
%     fptr = fopen(subfileName,'w');
%     fprintf(fptr,'#!/bin/bash\n');
% %     fprintf(fptr,'#SBATCH -A SNIC2020-7-19\n');
%     fprintf(fptr,'#SBATCH -A C3SE2020-1-8\n');
%     fprintf(fptr,'#SBATCH -N 1\n');
%     fprintf(fptr,'#SBATCH -n 20\n');
%     fprintf(fptr,'#SBATCH -o out.txt\n');
%     fprintf(fptr,'#SBATCH -t 50:00:00\n');
%     fprintf(fptr,'#SBATCH --mail-user=cheyu@chalmers.se\n');
%     fprintf(fptr,'#SBATCH --mail-type=end\n');
%     fprintf(fptr,'module load MATLAB intel/2018b GMP\n');
%     fprintf(fptr,['i=',num2str(k(i)),'\n']);
%     if i ~= length(k)
%         for m = 1:20
%             fprintf(fptr,['a',num2str(m),'=$i+' num2str(1*(m-1)) '\n']);
%             fprintf(fptr,['b',num2str(m),'=$i+' num2str(1*m-1) '\n']);
%         end
%         for m = 1:20
%             fprintf(fptr,['matlab -nodesktop -singleCompThread -r "simulationCNPS($a',num2str(m),',$b',num2str(m),')" &\n']);
%         end
%     else
%         for m = 1:kmax-k(i)+1
%             fprintf(fptr,['a',num2str(m),'=$i+' num2str(1*(m-1)) '\n']);
%             fprintf(fptr,['b',num2str(m),'=$i+' num2str(1*m-1) '\n']);
%         end
%         for m = 1:kmax-k(i)+1
%             fprintf(fptr,['matlab -nodesktop -singleCompThread -r "simulationCNPS($a',num2str(m),',$b',num2str(m),')" &\n']);
%         end
%     end
%     fprintf(fptr,'wait;\n');
%     fclose(fptr);
% end


%% simulationAnaerobicCNPS
% [~,raw,~] = xlsread('Exp_biolog.xlsx');
% metNameID = raw(2:end,2);
% 
% kmax = sum(~ismember(metNameID,''));
% k = 1:20:kmax;
% for i = 1:length(k)
%     subfileName = ['sub',num2str(k(i)),'.sh'];
%     fptr = fopen(subfileName,'w');
%     fprintf(fptr,'#!/bin/bash\n');
% %     fprintf(fptr,'#SBATCH -A SNIC2020-7-19\n');
%     fprintf(fptr,'#SBATCH -A C3SE2020-1-8\n');
%     fprintf(fptr,'#SBATCH -N 1\n');
%     fprintf(fptr,'#SBATCH -n 20\n');
%     fprintf(fptr,'#SBATCH -o out.txt\n');
%     fprintf(fptr,'#SBATCH -t 50:00:00\n');
%     fprintf(fptr,'#SBATCH --mail-user=cheyu@chalmers.se\n');
%     fprintf(fptr,'#SBATCH --mail-type=end\n');
%     fprintf(fptr,'module load MATLAB intel/2018b GMP\n');
%     fprintf(fptr,['i=',num2str(k(i)),'\n']);
%     if i ~= length(k)
%         for m = 1:20
%             fprintf(fptr,['a',num2str(m),'=$i+' num2str(1*(m-1)) '\n']);
%             fprintf(fptr,['b',num2str(m),'=$i+' num2str(1*m-1) '\n']);
%         end
%         for m = 1:20
%             fprintf(fptr,['matlab -nodesktop -singleCompThread -r "simulationAnaerobicCNPS($a',num2str(m),',$b',num2str(m),')" &\n']);
%         end
%     else
%         for m = 1:kmax-k(i)+1
%             fprintf(fptr,['a',num2str(m),'=$i+' num2str(1*(m-1)) '\n']);
%             fprintf(fptr,['b',num2str(m),'=$i+' num2str(1*m-1) '\n']);
%         end
%         for m = 1:kmax-k(i)+1
%             fprintf(fptr,['matlab -nodesktop -singleCompThread -r "simulationAnaerobicCNPS($a',num2str(m),',$b',num2str(m),')" &\n']);
%         end
%     end
%     fprintf(fptr,'wait;\n');
%     fclose(fptr);
% end





