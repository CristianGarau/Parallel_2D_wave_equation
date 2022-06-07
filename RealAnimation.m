sizeX = 50;
sizeY = sizeX;
sizeZ = 0.2;

fileNumber = 100;
fileName = "results/data";
data = {};



for c = 1:fileNumber
    tmp = append(fileName, int2str(c));
    tmp = append(tmp, ".txt");
    data{c} = readtable(tmp);
    data{c} = data{c}(1:sizeX, 1:sizeY);
    data{c} = table2array(data{c});
end

[X,Y] = meshgrid(1:1:sizeX,1:1:sizeY);
C = X.*Y;

for j= 1:fileNumber
    j
    p1 = surf(X,Y,data{j},C);
    v = [36.883559342517408,41.399999999999999];
    [caz,cel] = view(v);
    axis ([0 sizeX 0 sizeY -sizeZ sizeZ]);
    colorbar
    
    
    pause(0.0001);
    hold on; 

    if(j<fileNumber)
    delete(p1);
    end
end