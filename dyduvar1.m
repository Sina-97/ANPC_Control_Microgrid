function [sys1,x01,str1,ts1] = dyduvar1(t1,x1,u1,flag1,Nu1,Ni1,Nj1,Ts1,minp1,maxp1,mint1,maxt1,Normalize1)
%DYDUVAR This function calculates the partial derivative of the 
%   output Y respect to the input U. It's intended to be used with
%   the function ptest2sim2 of the Neural Network Predictive Controller.
%   
%   See sfuntmpl.m for a general S-function template.
%
%   See also SFUNTMPL.
    
% Orlando De Jesus, Martin Hagan, 1-30-00
% Copyright 1992-2010 The MathWorks, Inc.

switch flag1,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
   [sys1,x01,str1,ts1]=mdlInitializeSizes(Nu1,Ni1,Nj1,Ts1);
  
  %%%%%%%%%%  
  % Update %
  %%%%%%%%%%
  case 2,                                               
    sys1 = mdlUpdate(t1,x1,u1,Nu1,Ni1,Nj1);
    
  %%%%%%%%%%
  % Output %
  %%%%%%%%%%
  case 3,  
    sys1 = mdlOutputs(t1,x1,u1,Nu1,Ni1,Nj1,minp1,maxp1,mint1,maxt1,Normalize1);    

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9,                                               
    sys1 = [];

  otherwise
    nnerr.throw('Control',['unhandled flag = ',num2str(flag1)]);
end

%end sfundsc1

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys1,x01,str1,ts1]=mdlInitializeSizes(Nu1,Ni1,Nj1,Ts1)

sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = Nu1*(Ni1+Nj1);
sizes.NumOutputs     = Nu1;
sizes.NumInputs      = -1;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;

sys1 = simsizes(sizes);

x01  = zeros(Nu1*(Ni1+Nj1),1);
x01(1)=1;
str1 = [];
ts1  = [Ts1 0]; % Inherited sample time

% end mdlInitializeSizes

%
%=======================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=======================================================================
%
function sys1 = mdlOutputs(t1,x1,u1,Nu1,Ni1,Nj1,minp1,maxp1,mint1,maxt1,Normalize1)
sys1 = x1(Ni1+1);
for k=1:Nu1-1
   sys1 = [sys1;x1(k*(Ni1+Nj1)+Ni1+1)];%u;    
   if Normalize1
      sys1=sys1*(maxt1-mint1)/(maxp1-minp1);
   end
end

%end mdlUpdate

%
%=======================================================================
% mdlOutputs
% Return the output vector for the S-function
%=======================================================================
%
function sys1 = mdlUpdate(t1,x1,u1,Nu1,Ni1,Nj1)

kk=(Nu1-1)*(Ni1+Nj1);
out=u1(1:Ni1)'*x1(kk+1:kk+Ni1)+u1(Ni1+1:Ni1+Nj1)'*x1(kk+Ni1+1:kk+Ni1+Nj1);
x1(kk+2:kk+Ni1)=x1(kk+1:kk+Ni1-1);           
if x1((Nu1-2)*(Ni1+Nj1)+1)==1 | x1(kk+1)==1
  x1(kk+1)=1;            
else
  x1(kk+1)=0;
end
x1(kk+Ni1+2:kk+(Ni1+Nj1))=x1(kk+Ni1+1:kk+(Ni1+Nj1)-1);           
x1(kk+Ni1+1)=out;    

for k=Nu1-1:-1:1
  kk=(k-1)*(Ni1+Nj1);
  out=u1(1:Ni1)'*x1(kk+1:kk+Ni1)+u1(Ni1+1:Ni1+Nj1)'*x1(kk+Ni1+1:kk+Ni1+Nj1);
  if k~=1
    x1(kk+1:kk+Ni1)=x1((k-2)*(Ni1+Nj1)+1:(k-2)*(Ni1+Nj1)+Ni1);
  else
    x1(kk+2:kk+Ni1)=x1(kk+1:kk+Ni1-1);           
    x1(kk+1)=0;            
  end
  x1(kk+Ni1+2:kk+Ni1+Nj1)=x1(kk+Ni1+1:kk+Ni1+Nj1-1);           
  x1(kk+Ni1+1)=out;    
end


sys1=x1;

%end mdlUpdate


