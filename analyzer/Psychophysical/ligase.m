close all;
clear all;

filePath = '/Users/fradeshieh/Desktop/matlab/a1/test';
file = dir(fullfile(filePath,'retinalStabilization_*.mat'));
flipNameStr = 'flip';
for i = 1:length(file)
    filename = file(i).name;
end
l=length(file)
for i = 1:l
    namei = file(i).name(1:length(file(i).name)-5)
    for j = i:l
        namej = file(j).name(1:length(file(j).name)-5);
        if namei == namej
            a=load(fullfile(filePath,file(i).name));
            b=load(fullfile(filePath,file(j).name));
            a.choice=cat(1,a.choice,b.choice);
            a.choiceTime=cat(1,a.choiceTime,b.choiceTime);
            a.Conditions=cat(1,a.Conditions,b.Conditions);
            save(namei,'-struct', 'a')
            continue;
        end
    end
end
