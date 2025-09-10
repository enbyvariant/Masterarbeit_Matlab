function [] = benchmarks_qp()
% Save benchmarks for all QPs in .json file

cd('/Users/loki/Documents/MATLAB/Masterarbeit')
fid = fopen('QPs-new.json', 'w');

lps = '/Users/loki/Documents/MATLAB/Masterarbeit/QPs';
myFiles = dir(lps);
myFiles = myFiles(3:end, :);
len = length(myFiles);
cd(lps);
for k = len:-1:1
    baseFileName = myFiles(k).name;
    fprintf(1, 'Now reading %s\n', baseFileName);
    
    cd(baseFileName);
    [nlp, dim, f] = get_input();
    if f
        data.toobig = 1;
    else
        [~, ~, data] = Interior_Points_QP(30, nlp, dim);
        data.p_name = baseFileName;
        encoded = jsonencode(data);
        fprintf(fid,'%s',encoded);
    end
    cd('..')
end
fclose(fid);
end