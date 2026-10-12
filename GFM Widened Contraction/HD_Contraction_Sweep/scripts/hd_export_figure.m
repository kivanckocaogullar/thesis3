function hd_export_figure(fig,base)
%HD_EXPORT_FIGURE Save a figure as 300-dpi PNG, vector PDF and editable FIG.
% PNG/PDF are written through an ASCII temporary folder (the project path contains
% non-ASCII characters, which the PNG/PDF writers do not always accept), then copied.
d=fileparts(base);if ~isfolder(d),mkdir(d);end
tmpd='C:\Windows\Temp';if ~isfolder(tmpd)||~canwrite(tmpd),tmpd=tempdir;end
[~,id]=fileparts(tempname);tp=fullfile(tmpd,[id '.png']);tf=fullfile(tmpd,[id '.pdf']);
c=onCleanup(@()cellfun(@(f)delete_if(f),{tp,tf}));
exportgraphics(fig,tp,'Resolution',300);copyfile(tp,[base '.png'],'f');
exportgraphics(fig,tf,'ContentType','vector');copyfile(tf,[base '.pdf'],'f');
savefig(fig,[base '.fig']);
end
function ok=canwrite(d)
f=fullfile(d,['w' char(java.util.UUID.randomUUID) '.tmp']);fid=fopen(f,'w');ok=fid>0;if ok,fclose(fid);delete(f);end
end
function delete_if(f),if isfile(f),delete(f);end,end
