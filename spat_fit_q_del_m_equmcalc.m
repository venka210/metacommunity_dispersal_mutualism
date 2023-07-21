clear all

%% parameter definition
r_x = 5; r_y = 5;

alpha_xy = 0.73; alpha_yx = 0.60; %really only affects y's local density and dispersal

K_x = 200; K_y = 200; %no point touching this

del_x = 0.01; del_y = 0.05;
%q = 0.45;
d_m = 0.3; 
a = 2.0; %capture rate of plant by frugivore
del_m = 8:0.02:10; %I'm only not starting from zero because the computational costs are absurd
%a = 0.51:0.01:0.81;%0.7-1.1 seems to work for this fig.
q = 0.91:0.01:0.94;
[Del,Q] = meshgrid(del_m,q);

Del_col = Del(:); Q_col = Q(:);

num_combinations = numel(Q_col); 

k_eff = 1; %efficiency of dispersing seeds to habitable patches

z_x = 0.7; z_y = 0.7; z_m = 0.4; %scaling factors for patch extinction rates. changing z_m relative to z_x and z_m does not change qual. change results

e_xmin = 0.05;e_ymin = 0.05; e_mmin = 0.05; e_mxmin = e_xmin;

tspan = [0,1000];

k_x = 0.08; k_y = 0.08; k_m = 0.08;

x_init = 0.1; y_init = 0.1; m_init = 0.1;
    
spp_init_no_m = [x_init y_init 0];
spp_init = [x_init y_init m_init];

%% Stable coexistence checks
%Two conditions need to be satisfied (one if mutualist is present and one
%if mutualist is not) to assume non-trivial equlibrium pop sizes at each
%patch. Else exclusion could occur.

%when mutualist is absent
no_m_upp_lim =  (r_x-del_x)/(alpha_xy*r_x); no_m_low_lim = (alpha_yx*r_y)/(r_y-del_y);
if (K_y/K_x) > no_m_low_lim && (K_y/K_x) < no_m_upp_lim
    disp('stable local coexistence w/o mutualist')
else
    error('no stable local coexistence w/o mutualist')
end
%when mutualist is present
yes_m_upp_lim = K_x*(r_x-del_x)/(alpha_xy*r_x);
yes_m_low_lim = (alpha_yx*r_y.*(Del_col+d_m))./(a.*Q_col*(r_y-del_y));

if all(K_y > yes_m_low_lim) && all(K_y < yes_m_upp_lim)
    disp('stable local coexistence w/ mutualist')
else
    error('no stable local coexistence w/ mutualist')
end

%% Local dynamics - equilibrium expressions assuming no competitive exclusion occurs

%when there is no mutualist on a patch
x_tilde = (K_x.*(1 - (del_x/r_x)) - alpha_xy.*(K_y.*(1-(del_y/r_y))))./(1-(alpha_xy)*(alpha_yx));
y_tilde = K_y.*(1 - (del_y/r_y)) - alpha_yx.*x_tilde;

%when mutualist is present on the patch

x_star = (Del_col' + d_m)./(a*Q_col');
y_star = K_y.*(1 - (del_y/r_y)) - alpha_yx.*x_star;
m_star = (1/a).*((r_x/K_x).*(K_x-x_star-(alpha_xy.*y_star))-del_x);

if (any(x_star > K_x) | any(x_star < 0) | any(imag(x_star)~=0))
    % If condition is true, throw an error and stop execution
    error('x popln is out of bounds');
end

if (any(y_star > K_y) | any(y_star < 0) | any(imag(y_star)~=0))
    error('y popln is out of bounds');
end

if (any(m_star > K_x) | any(m_star < 0) | any(imag(m_star)~=0))
    error('m popln is out of bounds');
end

%% Linking local and global dynamics

%Extinction equations
e_x = e_xmin.*((K_x./x_tilde).^z_x);e_y = e_ymin.*(K_y./(y_star)).^z_y;
e_mx = e_xmin.*(K_x./(x_star)).^z_x;
e_m = e_mmin.*(K_x./(m_star)).^z_m; %assume max population size of mutualist is similar to that of x and y

if (any(e_m < 0) | any(imag(e_m)~=0))
    % If condition is true, throw an error and stop execution
    error('e_m values are not realistic');
end


mu = e_mx - e_x;

%Colonisation equations
c_x = (k_x*del_x).*x_tilde;
c_y = k_y*del_y.*y_star;
c_m = k_m*Del_col'.*m_star; %repelem(del_m,numel(q)) can also be used in place of Del_col' if u feel funky
c_mx =(k_x*del_x+k_m*k_eff*a.*(1-Q_col').*Del_col'.*m_star).*x_star;

lambda = c_mx-c_x;

%% Global dynamics - equilibrium expressions

gamma = (0.5.*((lambda.*(1+(e_m./c_m))+c_x)-(e_x+mu)))./(c_x+lambda);%0.5*(1-((mu-e_x-e_m)/(c_x+c_m)));

eta2 = (e_m./c_m).*((mu-lambda)./(c_x+lambda));

px_star_pos = gamma + sqrt(gamma.^2+eta2);
px_star_neg = gamma - sqrt(gamma.^2+eta2);
py_star = 1 - (e_y./c_y);
pm_star = px_star_pos - (e_m./c_m);
 %Note: 'any' and 'all' commands work on something that is already in logical
 %form i.e. arrays if zeros and ones which are logical
 
if all(px_star_pos <= 1) && all(px_star_pos >= 0) && all(imag(px_star_pos) == 0)
    disp('roots of x make sense');
else
    disp('roots of x dont make sense');
end

if  all(pm_star <= 1) && all(pm_star >= 0) && all(imag(pm_star)== 0)
    disp('roots of m make sense');
else
    disp('roots of m dont make sense');
end

if (all(py_star <= 1) && all(py_star >= 0) && all(imag(py_star)==0))
    disp('roots of y make sense');
else
    disp('roots of y dont make sense');
end

px_star_old = px_star_pos; py_star_old = py_star; pm_star_old = pm_star;

px_star_pos(px_star_old < 0 | px_star_old > 1 | imag(px_star_old) ~= 0) = 0;%technically this should be whatever the p* of x is in a non-mutualist patch (which is 0 in this case of local exclusion)
pm_star(pm_star_old < 0 | pm_star_old > 1 | imag(pm_star_old) ~= 0) = 0;
%% creating matrix for plotting

px_star_pos_plot(px_star_pos > py_star) = 1; px_star_pos_plot(px_star_pos <= py_star) = -1;%this is only done to make the plot look pretty
plotdata = reshape(px_star_pos_plot,[numel(q),numel(del_m)]);% - reshape(py_star,[numel(q),numel(del_m)]);
plotdata_diff = reshape((px_star_pos-py_star),[numel(q),numel(del_m)]);

%% plotting

figure()
h = heatmap(del_m,fliplr(q), flipud(plotdata));
colormap parula
customlabels_x = string(del_m);
customlabels_x(mod(del_m,0.5) ~= 0) = " ";
h.XDisplayLabels = customlabels_x;
% surf(Del, Q, reshape(frac_occup_3d(1,:),numel(q),numel(del_m))) colormap
% parula; freezeColors; hold on surf(Del, Q,
% reshape(frac_occup_3d(2,:),numel(q),numel(del_m))) colormap winter;
% freezeColors;
xlabel('mutualist dispersal rate (\delta_m)')
ylabel ('consumption fraction (q)')
%zlabel('fraction of patches occupied')
title('Fraction of patches occupied vs mutualist dispersal rate and predation rate')
%legend('Species with mutualist (x)')%, 'Species without mutualist (y)',
%'mutualist (m)', 'location', 'best' )
print('spat_fit_q_del_m_equmcalc','-djpeg','-r600')
   
    
