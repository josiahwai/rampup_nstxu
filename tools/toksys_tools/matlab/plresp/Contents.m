%DIII-D Plasma Response Functions and Data
%
% (/users/humphrys/matlab/plresp)
%
% rrigd3d - Make OLD format rz plasma response model and calc eigenvalues
% calc_rzresp - Make NEW format rz plasma response model and save set
%		(drzdi.mat)
% compare_isobgd_rig - Compare EFIT isoflux/Bgrid responses to model-derived
%      => best general comparison code: does radial/vertical comparison too!
% compare_efit_rig - Compare EFIT radial/vertical responses to model-derived
% calc_rzresp - Calculates zr response model (replaces rrigd3d)  8/98
% plasma_dynamics - Calculates dynamic plasma response model objects
%
% Doc file:
%   /users/humphrys/matlab/plresp/plresp.info
%
% Other plasma response model codes connecting to the scripts in plresp can
%  be found in directory:
%   /users/walker/matlab/plasma_models/*.m
%
% Unix scripts usedwith these Matlab scripts and functions:
%	make_modelDAH - DAH version of Mike's make_model. Does entire build
%		of model objects and files needed for doing MIMO design. Starts
%		up Matlab as needed and generates everything.
%
