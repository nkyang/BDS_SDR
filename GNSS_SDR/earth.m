load topo;
[x,y,z]= sphere(45);
s=surface(x,y,z,'FaceColor','texturemap','CData',topo);
colormap(topomap1);
brighten(.6);
campos([1.3239 -14.4250 9.4954]);
camlight;
lighting gouraud
% axis off vis3d