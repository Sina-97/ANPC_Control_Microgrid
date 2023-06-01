function [JJ,dJJ]=calcjjdjj1(u_vec,Ni1,Nu1,Nj1,N21,Ai,Ts1,ref,tiu,rho1,dUtilde_dU1,Normalize1,minp1,maxp1)
%CALCJJDJJ Calculate the cost function JJ and the change in the cost function dJJ
%with respect to U for the NN Predictive Controller.
%
%  Synopsis
%
%    [JJ,dJJ]=calcJJdJJ(u_vec,Ni,Nu,Nj,N2,Ai,Ts,ref,tiu,rho,dUtilde_dU,Normalize)
%
%  Description
%
%    This function predict the output of the system based on a previously trained
%    NN. With that information, calculates the cost function and the change in the
%    cost function based on the inputs to the plant. The cost function is:
%
%    JJ = evec'*evec + rho*(duvec'*duvec);
%
%    where: evec = Error during the cost horizon time.
%           duvec = Change in the control action during the control horizon time.
%           rho = Control weighting factor.
%
%    [JJ,dJJ]=CALCJJDJJ(u_vec,Ni,Nu,Nj,N2,Ai,Ts,ref,tiu,rho,dUtilde_dU,Normalize) takes,
%      u_vec - Vector with sequence of control actions.
%      Ni    - Number of delayed plant inputs.
%      Nu    - Control Horizon.
%      Nj    - Number of delayed plant outputs.
%      N2    - Cost Horizon.
%      Ai    - Initial layer delay conditions.
%      Ts    - Time steps.
%      ref   - Reference input.
%      tiu   - Initial time for U.
%      rho   - Control weighting factor.
%      dUtlde_dU - Derivate of the difference of U(t)-U(t-1) respect U.
%      Normalize - Indicate if the NN has input-output normalized.
%    and returns,
%      JJ    - Cost function value.
%      dJJ   - Derivate of the cost function respect U.
%
%    This function is being called by the line search functions (CSRCHBAC, CSRCHGOL, 
%    CSRCHHYB, CSRCHCHA, CSRCHBRE) inside the function PREDOPT, that is located in
%    the Simulink block for the NN Predictive Controller.

% Orlando De Jesus, Martin Hagan, 1-30-00
% Copyright 1992-2002 The MathWorks, Inc.

assignin('base','cont_u',evalin('base','cont_u')+1);

set_param('ptest3sim21/Subsystem','u_init1',mat2str(u_vec(Ni1),20),'ud_init1',mat2str(u_vec(Ni1-1:-1:1),20), ...
                                  'y_init1',mat2str(Ai{2,Nj1},20),'yd_init1',mat2str(cat(1,Ai{2,Nj1-1:-1:1}),20));
[time,xx0,Ac1,Ac2,E,gU,gUd,dY_dU] = sim('ptest3sim21',[0 N21*Ts1],[],[(0:Ts1:(N21-2)*Ts1)' u_vec(1:N21-1) ref(1:N21-1)]);

yhat_vec=Ac1(1:N21+1,1)';

E=E(2:N21+1,:);
gU=gU(1:N21,:)';
gUd=gUd(1:N21,:)';

evec=E;

if tiu==1
   duvec = [0; u_vec(tiu+1:tiu+Nu1-1)-u_vec(tiu:tiu+Nu1-2)];
else   
   duvec = u_vec(tiu:tiu+Nu1-1)-u_vec(tiu-1:tiu+Nu1-2);
end
 
JJ = evec'*evec + rho1*(duvec'*duvec);

% Forward Perturbation
dY_dU=dY_dU(2:N21+1,:)';
dJJ   = 2*(-dY_dU*evec + rho1*(dUtilde_dU1*duvec));
if Normalize1
   dJJ=dJJ/(maxp1-minp1);
end