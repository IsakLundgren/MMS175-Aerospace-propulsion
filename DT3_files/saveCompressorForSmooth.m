function saveCompressorForSmooth(NH, mdot, tin, pin, pout,x_igv_le, x_igv_te, x_rot_le, x_rot_te, x_stat_le, x_stat_te,...
                            r_igv_tle, r_igv_tte, r_rot_tle, r_rot_tte, r_stat_tle, r_stat_tte, r_igv_hle, r_igv_hte,...
                            r_rot_hle, r_rot_hte, r_stat_hle, r_stat_hte)
%saveCompressorForSmooth File writes output for smooth
%   Detailed explanation goes here

fid = openSmoothFile('sc90c_files/smoothInput.txt');

prat = pout/pin;
fprintf(fid, '%s %4.4f \n','prat:',prat);
fprintf(fid, '%s %4.4f \n','NH:',60.0*NH);
fprintf(fid, '%s %4.4f \n','mdot:',mdot);
fprintf(fid, '%s %4.4f \n','T0:',tin);
fprintf(fid, '%s %4.4f \n','P0:',pin);
fprintf(fid, '%s %4.4f \n','reDistFact:',1.0);
fprintf(fid, '%s %i \n','n_stages:',1);
fprintf(fid, '%s %i \n','n_span:',11); % must be odd number
fprintf(fid, '%s %i \n','fl_igv:',1);

n_bezier_pts = 4; % one at inlet, one at rotor trailing edge, one at outlet

fprintf(fid, '%s %i \n','n_bezier:',n_bezier_pts);
fprintf(fid, '%s %i \n','n_duct_inlet:',2);
fprintf(fid, '%s %i \n','n_duct_outlet:',3);
fprintf(fid, '%s %i \n','n_blade_intermediate_stations:',0);
logVec(fid,'hubStretch: ',3,ones(1,3))
logVec(fid,'tipStretch: ',3,ones(1,3))
logVec(fid,'stretch: ',3,ones(1,3))
logVec(fid,'spacing: ',3,ones(1,3))
logVec(fid,'hubSpacing: ',3,ones(1,3))
logVec(fid,'tipSpacing: ',3,ones(1,3))
fprintf(fid, '%s %i %i %i \n','blades: ',[70 100 150]);


[xbezHub,rbezHub,xbezShroud,rbezShroud] = getBezierData(n_bezier_pts, x_igv_le,x_igv_te, x_rot_le, x_stat_le,...
                                                                 x_stat_te, r_igv_tle, r_rot_tte, r_stat_tle,...
                                                                 r_stat_tte, r_igv_hle, r_rot_hte, r_stat_hle, r_stat_hte);

logVec(fid,'xbezHub: ',n_bezier_pts,xbezHub)
logVec(fid,'rbezHub: ',n_bezier_pts,rbezHub)
logVec(fid,'xbezShroud: ',n_bezier_pts,xbezShroud)
logVec(fid,'rbezShroud: ',n_bezier_pts,rbezShroud)
logVec(fid,'alpha_TE_S1: ',3,[0.0  0.0  0.0])
logVec(fid,'alpha_TE_IGV: ',3,[0.0 0.0 0.0])
fprintf(fid, '%s %4.4f %4.4f %4.4f  \n','workFact_R1: ',ones(1,3));
logVec(fid,'xHub: ',6,[x_igv_le x_igv_te x_rot_le x_rot_te x_stat_le x_stat_te])
logVec(fid,'xTip: ',6,[x_igv_le x_igv_te x_rot_le x_rot_te x_stat_le x_stat_te])
fclose(fid);

end

function logVec(fid,str,npts,vec)

fprintf(fid,'%s',str);
for i = 1:npts
  fprintf(fid, '%8.4f',vec(i));    
end
fprintf(fid,'%s \n','');



end 

function [xbezHub,rbezHub,xbezShroud,rbezShroud] = getBezierData(n_bezier_pts,...
    x_igv_le,x_igv_te,x_rot_le,x_stat_le,x_stat_te, r_igv_tle,r_rot_tte,r_stat_tle,r_stat_tte,...
    r_igv_hle,r_rot_hte,r_stat_hle,r_stat_hte)

if n_bezier_pts == 4 
% define xpts for hub and shroud bezier curve
   axChord = x_igv_te - x_igv_le;
   x1 = x_igv_le - 0.5*axChord; % one blade chord upstream of R1 leading edge
   x11 = x_rot_le;
   x2 = x_stat_le; % R1 leading edge
   x3 = x_stat_te + 0.5*(x_stat_te - x_stat_le); % one stator chord downstream of S1 trailing edge
% number of x-stations above must equal n_bezier_pts
   
   xbezHub = [x1 x11 x2 x3];
   xbezShroud = xbezHub; % no sweep is assumed 

% bezier interpolation for shroud lines 
   r1 = r_igv_tle;  % constant radius one blade chord upstream of R1 leading edge
   r11 = r_igv_tle; % R1 leading edge
   r2 = r_rot_tte; % R1 trailing edge
   r3 = r_stat_tte; % S1 trailing edge
   
% shroud interpolation to smoothen wall contour
   for i = 1:n_bezier_pts; 
      rbezShroud(i) = get_y_bc(xbezShroud(i),xbezShroud,[r1 r11 r2 r3]);
   end

% bezier interpolation for hub lines 
   r1 = r_igv_hle; % HPCgeom.Rotor.rle(1,1) % constant radius one blade chord upstream of R1 leading edge
   r11 = r_rot_hte; % HPCgeom.Rotor.rle(1,1); % R1 leading edge
   r2 = r_rot_hte; % HPCgeom.Rotor.rte(1,1); % R1 trailing edge
   r3 = r_stat_hte; % HPCgeom.Stator.rle(1,1); % S1 trailing edge
   
% hub interpolation to smoothen wall contour
   for i = 1:n_bezier_pts; 
      rbezHub(i) = get_y_bc(xbezHub(i),xbezHub,[r1 r11 r2 r3]);
   end
    
end 
end 

function fid = openSmoothFile(fn)
  fid=fopen(fn,'w');
  copyfile sc90c_files/smoothInput.txt sc90c_files/smoothInputOriginal.txt 
end 


%
%
%   2013-10-27
%   Function that calculates the y-coordinate by interpolating using Bezier curves
%
%

function y_point = get_y_bc(x_point,Px,Py)

t_in = (x_point - Px(1))/(Px(end) - Px(1)); %initial guess

t = fzero(@(t) x_point - get_B_curve(Px,t),t_in);

y_point = get_B_curve(Py,t);
end

%   2014-10-08
%   Bezier curve
%
%
function B = get_B_curve(P,t)
n = length(P);

B = 0;
for i = 0:n-1
    B = B +  nchoosek(n-1,i) * (1-t).^(n-1-i) .* t .^i * P(i+1);
end
end