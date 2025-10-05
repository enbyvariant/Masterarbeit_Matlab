function [] = benchmarks_lp()
% Save benchmarks for all LPs in .json files

cd('/Users/loki/Documents/MATLAB/Masterarbeit/');
fid = fopen('LPs-new.json', 'w');
lps = '/Users/loki/Documents/MATLAB/Masterarbeit/LPs/';
myFiles = dir(lps);
cd(lps);
myFiles = myFiles(4:end, :);
len = length(myFiles);
data = [];
for k = 1:len
    baseFileName = myFiles(k).name;
    name = genvarname(baseFileName);
    fprintf(1, 'Now reading %s\n', baseFileName);
    cd(baseFileName);
    [nlp, dim] = get_input();
    [~, ~, data.(name)] = Interior_Points_LP_alt(30, nlp,dim);
    cd('..')
end
encoded = jsonencode(data);
fprintf(fid,'%s',encoded);
fclose(fid);
end
