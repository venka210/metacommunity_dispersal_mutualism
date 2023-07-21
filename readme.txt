25/4/23

I'm pausing this project to work on submitting the microbes manuscript. This doc is a summary of the functions/files that are in the folder

BetweenPatchDynamics_* --  Self explanatory and has multiple versions. All of these give single parameter set value outputs. The vectorized version is the one currently being
used 

vectorized_BetweenPatchDynamics -- Use this for large parameter ranges

disp_mutualism_metacommunity* -- The og version has only a del_m sweep. The one w/o 'og' has a nested loop and 
sweeps across predation rate and del_m

spat_fit_q_del_m_equmcalc --  In this code, I assume steady state has been reached at the local and global level
Then, I calculate the the regions of q-del_m space where x is fitter than y. Need to validate this by getting
outcomes similar to the vectorized_q_del_m_disp_mutualism. Very Important.

vectorized_q_del_m_disp_mutualism -- similar to above but actually solves the ODEs. Very Important!

vectorized_a_del_m_varied -- similar to above but in a-del_m space. Both this and the above need to be parametrized 
to obtain parameter space where this makes sense (and is stable)

new_params_coexist, new_params_coexist_nonneg --  only del_m sweep but shows coexistence when x is a poorer competitor
and disperser relative to y. The second one has a nonneg odeset option that makes it slightly better I guess

Lot of random pieces of code are present which where helpful in making a surface plot with two colours (jicolorbar,
unfreezeColors etc.)

NEXT STEPS:

-Get the analytical expression for spatial overlap and Beta-diversity discussed with Allison written down. Substitute
the expressions for px and py in there and write code to plot it wrt relevant parameters

- Find parameter range that stabilises the q-del_m and a-del_m plots and identify how to prevent e_m taking complex
values. Develop a schema that checks for this. Ultimately the vectorized_q_del_m_disp_mutualism and
spat_fit_q_del_m_equmcalc results need to kinda match










