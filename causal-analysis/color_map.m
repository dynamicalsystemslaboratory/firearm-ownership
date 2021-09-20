function map = color_map(X,X_unscaled)

set(groot,'defaultLegendInterpreter','latex')
cc=floor(10*X)-min(floor(10*X))+1; % define the colors range
c=bone(max(cc)); %generates RGB colors
c=flipud(c);
ax = usamap('all');
states = shaperead('usastatelo', 'UseGeoCoords', true); states(51)=[];
colormap bone
if sum(X)~=0
    colorbar('Direction','reverse','Ticks',[(min(X)),(max(X))],'TickLabels',{string(max(X_unscaled)),string(min(X_unscaled))},'FontName','Times New Roman','FontSize',20)
else
    colorbar('FontName','Times New Roman','FontSize',20)
end
faceColors = makesymbolspec('Polygon',{'INDEX', [1 numel(states)], 'FaceColor', c(cc,:)});
geoshow(ax(1), states, 'DisplayType', 'polygon', 'SymbolSpec', faceColors)
geoshow(ax(2), states, 'DisplayType', 'polygon', 'SymbolSpec', faceColors)
geoshow(ax(3), states, 'DisplayType', 'polygon', 'SymbolSpec', faceColors)

for k = 1:3
    setm(ax(k), 'Frame', 'off', 'Grid', 'off', 'ParallelLabel', 'off', 'MeridianLabel', 'off')
end

end