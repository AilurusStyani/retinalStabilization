close all;
clear all;

filePath = '/Users/fradeshieh/Desktop/matlab/a3/working';
file = dir(fullfile(filePath,'retinalStabilization_*.mat'));
flipNameStr = 'flip';
for i = 1:length(file)
    filename = file(i).name;
end
l=length(file);
trap=[0 0 0 0 0];
for i = 1:l
    namei = file(i).name(1:length(file(i).name)-5);
    k=length(namei);
    stop=0;
    for m = 1:5
        if str2num(namei(k))==trap(m)
            stop=1;
            continue;
            
        end
    end
    if stop
        continue;
    end
     a=load(fullfile(filePath,file(i).name));
    for j = i+1:l
        namej = file(j).name(1:length(file(j).name)-5);
        
        if  namei(k) == namej(k)
           
            b=load(fullfile(filePath,file(j).name));
            a.choice=cat(1,a.choice,b.choice);
            a.choiceTime=cat(1,a.choiceTime,b.choiceTime);
            a.Conditions=cat(1,a.Conditions,b.Conditions);
        end
        
    end
     save(['retinalStabilization_MIX',namei(k)],'-struct', 'a');
    trap(str2num(namei(k)))=str2num(namei(k));
end
