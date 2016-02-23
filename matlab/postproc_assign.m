close all; clear all; clc;
%%	Input the following data
file_data = '/Users/yuyang/temp/outdata.txt';   % name and path of the data file
file_tab = '/Users/yuyang/temp/outtab';     % name and path of the table file
file_conn = '/Users/yuyang/temp/rhodopsin_98/rho_connection.txt';   % name and path of the connection table
file_out = '/Users/yuyang/temp/out_cstab.txt';   % name and path of the output file with 90% probability result
group_size = 100;   % total number of results
N_tab = 2;      % number of spectra
N_freq = 5;     % number of frequencies in the table file 
                % (for example, N_freq is 3 if the chemical shifts of N, C' and Ca are included in the file)
%%
%%	Begin the post-processing
str0 = 'Ng, Nb, Ne, Nu and Pareto order:';
%
% read table file and collect Ng, Nb, Ne, Nu and Pareto order information
evl = [];
for k = 1:group_size
    i = floor(k/100);
    j = floor((k-i*100)/10);
    p = k-i*100-j*10;
    file1 = [file_tab,num2str(i),num2str(j),num2str(p),'.txt'];
    fid1 = fopen(file1,'r');
    line = fgets(fid1);
    ind = findstr(line, str0);
    if ~isempty(ind)
        B = line(ind+length(str0):end);
        evl = [evl;str2num(dblnk(B))];
    end
    fclose(fid1);
end
%
ind_pick = find(evl(:,5) == 1);
Ngbeu_pick = evl(ind_pick,1:4);
%
% read data file and collect Pareto-order-1 solutions
import_temp = importdata(file_data);
data = import_temp.data;
%
out_data = [];
for k = 1:length(ind_pick)
    out_data = [out_data;data(2*ind_pick(k)-1:2*ind_pick(k),:)];
end
N_seq = size(out_data,2);
%
fprintf('There are %d Pareto-order-1 solutions \n\n',length(ind_pick));
fprintf(['File number: ',num2str(ind_pick'),'\n\n']);
%
% calculate the maximum assignment probability for each residue
tables = cell(N_seq,N_tab);
n_assign = cell(N_seq,N_tab);
show_a = zeros(N_tab,N_seq);
value_max = show_a;
%
for k1 = 1:N_seq
    for k2 = 1:N_tab
        a = out_data(k2:N_tab:end,k1);
        a1 = unique(a);
        tables{k1,k2} = a1;
        %
        temp_l = [];
        for k3 = 1:length(a1)
            ind = find(a == a1(k3));
            temp_l = [temp_l,length(ind)];
        end
        n_assign{k1,k2} = temp_l;
        show_a(k2,k1) = max(temp_l)/sum(temp_l);
        ind = find(temp_l == max(temp_l));
        value_max(k2,k1) = a1(ind(1));
    end
end
%   Then tables{p1,p2} is an array which includes all the possible peaks
%   assigning to the (p1)th residue in the (p2)th spectrum
%
%   n_assign{p1,p2} is an array which includes the number of the
%   corresponding possible peaks assigning to the (p1)th residue in the
%   (p2)th spectrum
%
%   show_a(p2,p1) is the maximum assignment probability for the (p1)th
%   residue in the (p2)th spectrum
%
%   value_max(p2,p1) is the "most possible" assigned peak number for the
%   (p1)th residue in the (p2)th spectrum
%
%%  Output the statistical evaluation of the assignments
for k = 1:N_tab
    figure,subplot(2,1,1),plot(show_a(k,:),'*-','markersize', 7);
    ylabel('Probability');
    xlabel('Residue Number');
    axis([0,N_seq+1,0,1]);
    title(['Assignment result for spectrum ',num2str(k)],'fontsize',12,'fontweight','bold');
    subplot(2,1,2), plot(value_max(k,:),'*','markersize', 7);
    axis([0,N_seq+1,0,max(value_max(k,:))*1.2]);
    ylabel('Most possible assignments');
    xlabel('Residue Number');
    %
    ind_90 = find(show_a(k,:)>0.9);
    eval(['ind_90_',num2str(k),' = ind_90;']);
    fprintf('%d out of %d residues in the sequence are with larger than 90%% assignment probability in spectrum %d\n\n',...
        length(ind_90),N_seq,k)
end
%
%%  Calculate the average chemical shifts (with 90% probability) based on all the spectra information
% read the connection table and get the chemical shifts relationship between different spectra
[Tab1,Tab2,Col1,Col2,Sft1,Sft2] = textread(file_conn,'%d%d%d%d%d%d','headerlines',1);
Tab_info = cell(N_freq,1);
Sft_info = cell(N_freq,1);
for k = 1: N_freq+1
    ind1 = find(Col1 == k);
    ind2 = find(Col2 == k);
    [Tab_temp,ind_temp] = unique([Tab1(ind1),Tab2(ind2)]);
    Sft_temp = [Sft1(ind1),Sft2(ind2)];
    Sft_temp = Sft_temp(ind_temp);
    Tab_info{k} = Tab_temp;
    Sft_info{k} = Sft_temp;
end
% 
% read one of the Pareto-order-1 solutions
ind_file = ind_pick(1);
i = floor(ind_file/100);
j = floor((ind_file-i*100)/10);
p = ind_file-i*100-j*10;
file_str = [file_tab,num2str(i),num2str(j),num2str(p),'.txt'];
%
read_format = repmat(' %f ',1,N_freq);
read_format = strcat(' %d ',read_format);
read_format = repmat(read_format,1,N_tab);
read_format = strcat('%d %s ',read_format,' %s');
%
var_str = 'peak1,C1{1:N_freq}';
for k = 2: N_tab
    var_1 = [',peak',num2str(k),',C',num2str(k),'{1:N_freq}'];
    var_str = strcat(var_str,var_1);
end
var_str = ['[num_rsd,rsd,',var_str,',S]'];
%
eval([var_str,' = textread(file_str,read_format,''headerlines'',1);']);
%[num_rsd,rsd, peak1,C{1:N_freq},peak2,D{1:N_freq},S] = textread(file_str,read_format,'headerlines', 1);
%
cs = zeros(N_seq,N_freq);
cs_times = zeros(N_seq,N_freq);
for k = 1: N_seq
    for p1 = 1: N_freq
        Tab_temp = Tab_info{p1};
        Sft_temp = Sft_info{p1};
        if isempty(Tab_temp)
            continue;
        end
        for p2 = 1: length(Tab_temp)
            ind_tab = Tab_temp(p2);
            eval(['ind_90 = ind_90_',num2str(ind_tab),';']);
            if isempty(find(ind_90 == k))
                continue;
            end          
            shift = Sft_temp(p2);
            if (k+shift)<1 || (k+shift)>N_seq
                continue;
            end
            eval(['C_temp = C',num2str(ind_tab),'{p1};']);
            cs_temp = C_temp(k);
            if cs_temp>=1000 || cs_temp == 0
                continue;
            end
            cs(k+shift,p1) = cs(k+shift,p1)+cs_temp;
            cs_times(k+shift,p1) = cs_times(k+shift,p1)+1;
        end
    end
end
for k = 1:N_seq
    for p = 1:N_freq
        if cs_times(k,p)>0
            cs(k,p) = cs(k,p)/cs_times(k,p);
        end
    end
end
%% Output the average chemical shift information into an output text file
fileID = fopen(file_out,'w');
fprintf(fileID,'Automatic resonance assignment result with 90%% probability\n\n');
for k = 1:N_seq
    fprintf(fileID, '%d\t%s\t', num_rsd(k), rsd{k});
    for k1 = 1:N_freq
        fprintf(fileID,'%8.3f\t', cs(k,k1));
    end
    fprintf(fileID, '\n');
    %
end
fclose(fileID);        

    