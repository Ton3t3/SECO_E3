clear, clc, close all
%%
folder = "realOutputs";

% STRING CON TODOS LOS NOMBRES DE LOS FICHEROS %
files = dir(folder);
names_og = {files.name};
count = 0;
for i=1:length(names_og)
    if(contains(names_og{i},".txt")==1)
        count= count+1;
    end
end
names = strings(1,count);
count = 1;
for i=1:length(names_og)
    if(contains(names_og{i},".txt")==1)
        names(count) = names_og{i};
        count = count +1;
    end
end
% ·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-·-· %
for i=1:length(names)
    figure(i)
    [~,name,~] = fileparts(names(i));
    load(sprintf("%s/%s",folder, names(i)));

    % SUSTITUIR - POR _ %
    chr = convertStringsToChars(name);
    for k=1:length(chr)
        if (chr(k)=="-")
            chr(k)="_";
            break
        end
    end
    name=convertCharsToStrings(chr);
    % ·-·-·-·-·-·-·-·-· %
    eval(sprintf('t = %s(:,1);',name));
    eval(sprintf('y = %s(:,2);',name));
    plot(t,y)
    title(sprintf("Plot de fichero %s",names(i)))
    xlabel("t(s)")
    ylabel("y(t)")
end

k = length(names)+1;
figure(k)
lista_nombres=strings(1,length(names));
hold on
for i=1:length(names)
    [~,name,~] = fileparts(names(i));    
    lista_nombres(i)=name;
    % SUSTITUIR - POR _ %
    chr = convertStringsToChars(name);
    for k=1:length(chr)
        if (chr(k)=="-")
            chr(k)="_";
            break
        end
    end
    name=convertCharsToStrings(chr);
    % ·-·-·-·-·-·-·-·-· %
    eval(sprintf('t = %s(:,1);',name));
    eval(sprintf('y = %s(:,2);',name));
    plot(t,y)
end
hold off
title("Comparación ficheros")
xlabel("t(s)")
ylabel("y(t)")
legend(lista_nombres)