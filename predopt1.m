function [sys1,x01,str1,ts1] = predopt1(t1,x1,u1,flag1,N21,Ts1,Nu1,maxiter1,csrchfun1,rho1,alpha1,S11,IW1,LW1_21,LW2_11,B11,B21,Ni1,Nj1,min_i1,max_i1,minp1,maxp1,mint1,maxt1,Normalize1)
%PREDOPT Executes the Predictive Controller Approximation based on Gauss Newton.
%   
    
% Copyright 1992-2010 The MathWorks, Inc.
% Orlando De Jesus, Martin Hagan, 1-25-00

switch flag1,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    load_system('ptest3sim21');
    if Normalize1
       IW_gU1=((maxt1-mint1)/(maxp1-minp1))*IW1;
    else
       IW_gU1=IW1;
    end
    set_param('ptest3sim21/Subsystem','B21',num2str(B21,20),'B11',mat2str(B11,20),'LW2_11',mat2str(LW2_11,20), ...
                                      'LW1_21',mat2str(LW1_21,20),'IW1',mat2str(IW1,20),'IW_gU1',mat2str(IW_gU1,20), ...
                                      'Ts1',num2str(Ts1),'S11',num2str(S11),'Ni1',num2str(Ni1), ...
                                      'Nj1',num2str(Nj1),'minp1',num2str(minp1,20),'maxp1',num2str(maxp1,20), ...
                                      'minp1',num2str(minp1,20),'mint1',num2str(mint1,20),'maxt1',num2str(maxt1,20), ...
                                      'Normalize1',num2str(Normalize1),'Nu1',num2str(Nu1));
    assignin('base','t_init',cputime);
    assignin('base','cont_u',0);
    [sys1,x01,str1,ts1]=mdlInitializeSizes(N21,Ts1,Nu1,alpha1,S11,Ni1,Nj1,min_i1,max_i1);
  
  %%%%%%%%%%  
  % Update %
  %%%%%%%%%%
  case 2,      
      IW_gU1=IW1;
  
    set_param('ptest3sim21/Subsystem','B21',num2str(B21,20),'B11',mat2str(B11,20),'LW2_11',mat2str(LW2_11,20), ...
                                      'LW1_21',mat2str(LW1_21,20),'IW1',mat2str(IW1,20),'IW_gU1',mat2str(IW_gU1,20), ...
                                      'Ts1',num2str(Ts1),'S11',num2str(S11),'Ni1',num2str(Ni1), ...
                                      'Nj1',num2str(Nj1),'minp1',num2str(minp1,20),'maxp1',num2str(maxp1,20), ...
                                      'minp1',num2str(minp1,20),'mint1',num2str(mint1,20),'maxt1',num2str(maxt1,20), ...
                                      'Normalize1',num2str(Normalize1),'Nu1',num2str(Nu1));
    assignin('base','t_init',cputime);
    assignin('base','cont_u',0);
    sys1 = mdlUpdate(t1,x1,u1,N21,Ts1,Nu1,maxiter1,csrchfun1,rho1,alpha1,S11,Ni1,Nj1,min_i1,max_i1,minp1,maxp1,mint1,maxt1,Normalize1);
    
  %%%%%%%%%%
  % Output %
  %%%%%%%%%%
  case 3,  
    sys1 = mdlOutputs(t1,x1,u1,Nu1,Ni1);    

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9,                                               
    close_system('ptest3sim21',0);
    assignin('base','t_end',cputime);
    sys1 = [];

  otherwise
    nnerr.throw('Control',['unhandled flag1 = ',num2str(flag1)])
end

%end sfundscl

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys1,x01,str1,ts1]=mdlInitializeSizes(N21,Ts1,Nu1,alpha1,S11,Ni1,Nj1,min_i1,max_i1)

global tiu1 dUtilde_dU1
global N11 d1 alpha21 upi1 uvi1

sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = Ni1+Nu1-1+(S11+1)*(Nj1-1);
sizes.NumOutputs     = 1;
sizes.NumInputs      = -1;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;

sys1 = simsizes(sizes);

% State Index:
%
%             x(1:Ni-1) = Previous Plant input u - Controller output (Size Ni-1).
%                 x(Ni) = Actual Plant input u - Controller output (Size 1).
%       x(Ni+1:Nu+Ni-1) = Next Plant input u - Controller output (Size Nu-1).
%              x(Nu+Ni) = Previous NN 2nd layer output (Size 1). 
%   x(Nu+Ni+1:Nu+Ni+S1) = Previous NN 1st layer output (Size S1).
%
%   Last two variables will repeat in case of multiple outputs. Not tested yet.
%
x01  = zeros(Ni1+Nu1-1+(S11+1)*(Nj1-1),1);
% ODJ 1-31-00 We place initial Plant input u - Controller output at mid range.
x01(Ni1:Nu1+Ni1-1) = (max_i1-min_i1)/2;
str1 = [];
ts1  = [Ts1 0]; % Inherited sample time

tiu1=Ni1;
dUtilde_dU1 = eye(Nu1);
dUtilde_dU1(1:Nu1-1,2:Nu1)=dUtilde_dU1(1:Nu1-1,2:Nu1)-eye(Nu1-1);
N11=1;
d1=1;
alpha21     = alpha1*alpha1;
upi1   = [1:Nu1-1 Nu1(ones(1,N21-d1-Nu1+2))];
uvi1   = [tiu1:N21-N11+Ni1];

% end mdlInitializeSizes

%
%=======================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=======================================================================
%
function sys1 = mdlOutputs(t1,x1,u1,Nu1,Ni1)
sys1 = x1(Ni1);

%end mdlUpdate

%
%=======================================================================
% mdlOutputs
% Return the output vector for the S-function
%=======================================================================
%
function sys1 = mdlUpdate(t1,x1,u1,N21,Ts1,Nu1,maxiter1,csrchfun1,rho1,alpha1,S11,Ni1,Nj1,min_i1,max_i1,minp1,maxp1,mint1,maxt1,Normalize1)

global tiu1 dUtilde_dU1
global N11 d1 alpha21 upi1 uvi1

Ai=num2cell(zeros(2,Nj1));
for k=1:Nj1-1
  Ai{1,k}=x1(Nu1+Ni1+1+(k-1)*(S11+1):Nu1+Ni1+S11+(k-1)*(S11+1));
  Ai{2,k}=x1(Nu1+Ni1+(k-1)*(S11+1));                             % delayed plant output
end
Ai{1,Nj1}=u1(4:3+S11);

ref(1:N21,1)=u1(1);
initval = '[upmin(Nu1)]';

upmin=[x1(Ni1+1:Nu1+Ni1-1);x1(Nu1+Ni1-1)];
u_vec(1:Ni1-1,1)=x1(2:Ni1);
if Normalize1
   ref=((ref-mint1)*2/(maxt1-mint1)-1);
   Ai{2,Nj1}=((u1(3)-mint1)*2/(maxt1-mint1)-1);           % Actual NN output
   upmin=((upmin-minp1)*2/(maxp1-minp1)-1); 
   u_vec=((u_vec-minp1)*2/(maxp1-minp1)-1); 
else
   Ai{2,Nj1}=u1(3);
end

upmin0   = upmin;             
einitval = eval(initval);     % Evaluate inival string

for tr=1:length(einitval),
  up=upmin0;                  % Initial value for numerical search for a new u  
  up(Nu1)=einitval(tr);
  u_vec(uvi1,1) = up(upi1);  
  dw = 1;                     % Flag specifying that up is new
  lambda = 0.1;               % Initialize Levenberg-Marquardt parameter
  
  
  %>>>>>>>>>>>>>>> COMPUTE PREDICTIONS FROM TIME t+N1 TO t+N2 <<<<<<<<<<<<<<<<
  assignin('base','cont_u',evalin('base','cont_u')+1);

  set_param('ptest3sim21/Subsystem','u_init1',mat2str(u_vec(Ni1),20),'ud_init1',mat2str(u_vec(Ni1-1:-1:1),20), ...
                                  'y_init1',mat2str(Ai{2,Nj1},20),'yd_init1',mat2str(cat(1,Ai{2,Nj1-1:-1:1}),20));
  [time,xx0,Ac1,Ac2,E,gU,gUd,dY_dU] = sim('ptest3sim21',[0 N21*Ts1],[],[(0:Ts1:(N21-2)*Ts1)' u_vec(1:N21-1) ref(1:N21-1)]);

  yhat_vec=Ac1(1:N21+1,1)';

  E=E(2:N21+1,:);

  gU=gU(1:N21,:)';
  gUd=gUd(1:N21,:)';

  evec=E;

  if tiu1==1
     duvec = [0; u_vec(tiu1+1:tiu1+Nu1-1)-u_vec(tiu1:tiu1+Nu1-2)];
  else   
     duvec = u_vec(tiu1:tiu1+Nu1-1)-u_vec(tiu1-1:tiu1+Nu1-2);
  end
 
  JJ = evec'*evec + rho1*(duvec'*duvec);

  % Forward Perturbation
  dY_dU=dY_dU(2:N21+1,:)';
  dJJ   = 2*(-dY_dU*evec + rho1*(dUtilde_dU1*duvec));
  if Normalize1
    dJJ=dJJ/(maxp1-minp1);
  end
  
  %>>>>>>>>>>>>>>>>>>>>>>    EVALUATE CRITERION    <<<<<<<<<<<<<<<<<<<<<<
  J = JJ;
    
    
  %>>>>>>>>>>>>>>>>>>>>>>>>      DETERMINE dyhat/du       <<<<<<<<<<<<<<<<<<<<<<<<<

  %>>>>>>>>>>>>>>>>>>>>>>>>>>>>    DETERMINE dJ/du     <<<<<<<<<<<<<<<<<<<<<<<<<<<<
  dJdu   = dJJ;


  %>>>>>>>>>>>>>>>>>>>>>>    DETERMINE INVERSE HESSIAN    <<<<<<<<<<<<<<<<<<<<<<<<<
  B = eye(Nu1);                  % Initialize Hessian to I


  delta=1;
  tol=1/delta;
  ch_perf = J;      % for first iteration.
  %>>>>>>>>>>>>>>>>>>>>>>>     BEGIN SEARCH FOR MINIMUM      <<<<<<<<<<<<<<<<<<<<<<    
  for m = 1:maxiter1,
  
  
    %>>>>>>>>>>>>>>>>>>>>>>>   DETERMINE SEARCH DIRECTION   <<<<<<<<<<<<<<<<<<<<<<<
    dX = -B*dJdu;
    
    if dX'*dJdu>0    % We reset the gradient if positive.
        %>>>>>>>>>>>>>>>>>>>>>>    DETERMINE INVERSE HESSIAN    <<<<<<<<<<<<<<<<<<<<<<<<<
      B = eye(Nu1);                  % Initialize Hessian to I
      delta=1;
      tol=1/delta;
      ch_perf = J;      % for first iteration.
        %>>>>>>>>>>>>>>>>>>>>>>>   DETERMINE SEARCH DIRECTION   <<<<<<<<<<<<<<<<<<<<<<<
      dX = -B*dJdu;
    end

    if Normalize1
     switch csrchfun1,
      case 1, %'csrchgol',
        [up_delta,J,dJdu_old,dJdu,retcode,delta,tol]=csrchgol(up,u_vec,ref,Ai,Nu1,N11,N21,d1,Ni1,Nj1,dX,dJdu,J,dX'*dJdu,delta,rho1,dUtilde_dU1,alpha1,tol,Ts1,-1,1,Normalize1,minp1,maxp1);
      case 2  %'csrchbac1',
        [up_delta,J,dJdu_old,dJdu,retcode,delta,tol]=csrchbac1(up,u_vec,ref,Ai,Nu1,N11,N21,d1,Ni1,Nj1,dX,dJdu,J,dX'*dJdu,delta,rho1,dUtilde_dU1,alpha1,tol,Ts1,-1,1,Normalize1,minp1,maxp1);
     case 3  %'csrchhyb'
        [up_delta,J,dJdu_old,dJdu,retcode,delta,tol]=csrchhyb(up,u_vec,ref,Ai,Nu1,N11,N21,d1,Ni1,Nj1,dX,dJdu,J,dX'*dJdu,delta,rho1,dUtilde_dU1,alpha1,tol,Ts1,-1,1,Normalize1,minp1,maxp1);
      case 4  %'csrchbre'
        [up_delta,J,dJdu_old,dJdu,retcode,delta,tol]=csrchbre(up,u_vec,ref,Ai,Nu1,N11,N21,d1,Ni1,Nj1,dX,dJdu,J,dX'*dJdu,delta,rho1,dUtilde_dU1,alpha1,tol,Ts1,-1,1,Normalize1,minp1,maxp1);
      case 5  %'csrchcha'
        J_old=J;
        [up_delta,J,dJdu_old,dJdu,retcode,delta,tol]=csrchcha(up,u_vec,ref,Ai,Nu1,N11,N21,d1,Ni1,Nj1,dX,dJdu,J,dX'*dJdu,delta,rho1,dUtilde_dU1,alpha1,tol,Ts1,-1,1,Normalize1,minp1,maxp1,ch_perf);
        ch_perf = J - J_old;
      otherwise
        J_old=J;
        [up_delta,J,dJdu_old,dJdu,retcode,delta,tol]=feval(csrchfun1,up,u_vec,ref,Ai,Nu1,N11,N21,d1,Ni1,Nj1,dX,dJdu,J,dX'*dJdu,delta,rho1,dUtilde_dU1,alpha1,tol,Ts1,-1,1,Normalize1,minp1,maxp1,ch_perf);
        ch_perf = J - J_old;
     end
    else
     switch csrchfun1,
      case 1, %'csrchgol',
        [up_delta,J,dJdu_old,dJdu,retcode,delta,tol]=csrchgol(up,u_vec,ref,Ai,Nu1,N11,N21,d1,Ni1,Nj1,dX,dJdu,J,dX'*dJdu,delta,rho1,dUtilde_dU1,alpha1,tol,Ts1,min_i1,max_i1,Normalize1,minp1,maxp1);
      case 2  %'csrchbac1',
        [up_delta,J,dJdu_old,dJdu,retcode,delta,tol]=csrchbac1(up,u_vec,ref,Ai,Nu1,N11,N21,d1,Ni1,Nj1,dX,dJdu,J,dX'*dJdu,delta,rho1,dUtilde_dU1,alpha1,tol,Ts1,min_i1,max_i1,Normalize1,minp1,maxp1);
     case 3  %'csrchhyb'
        [up_delta,J,dJdu_old,dJdu,retcode,delta,tol]=csrchhyb(up,u_vec,ref,Ai,Nu1,N11,N21,d1,Ni1,Nj1,dX,dJdu,J,dX'*dJdu,delta,rho1,dUtilde_dU1,alpha1,tol,Ts1,min_i1,max_i1,Normalize1,minp1,maxp1);
      case 4  %'csrchbre'
        [up_delta,J,dJdu_old,dJdu,retcode,delta,tol]=csrchbre(up,u_vec,ref,Ai,Nu1,N11,N21,d1,Ni1,Nj1,dX,dJdu,J,dX'*dJdu,delta,rho1,dUtilde_dU1,alpha1,tol,Ts1,min_i1,max_i1,Normalize1,minp1,maxp1);
      case 5  %'csrchcha'
        J_old=J;
        [up_delta,J,dJdu_old,dJdu,retcode,delta,tol]=csrchcha(up,u_vec,ref,Ai,Nu1,N11,N21,d1,Ni1,Nj1,dX,dJdu,J,dX'*dJdu,delta,rho1,dUtilde_dU1,alpha1,tol,Ts1,min_i1,max_i1,Normalize1,minp1,maxp1,ch_perf);
        ch_perf = J - J_old;
      otherwise
        J_old=J;
        [up_delta,J,dJdu_old,dJdu,retcode,delta,tol]=feval(csrchfun1,up,u_vec,ref,Ai,Nu1,N11,N21,d1,Ni1,Nj1,dX,dJdu,J,dX'*dJdu,delta,rho1,dUtilde_dU1,alpha1,tol,Ts1,min_i1,max_i1,Normalize1,minp1,maxp1,ch_perf);
        ch_perf = J - J_old;
     end
    end

    
    %>>>>>>>>>>>>>>>>>>>>>>>>   UPDATE FUTURE CONTROLS   <<<<<<<<<<<<<<<<<<<<<<<<<
    up_old = up;
    up = up_delta; 
     
     
     %>>>>>>>>>>>>>>>>>>>>>>>>     CHECK STOP CONDITION     <<<<<<<<<<<<<<<<<<<<<<<
    dup = up-up_old;
    if (dup'*dup < alpha21) | (ch_perf==0),
      break;
    end 
       
       
     %>>>>>>>>>>>>>>>>>>>     BFGS UPDATE OF INVERSE HESSIAN    <<<<<<<<<<<<<<<<<<
    dG  = dJdu - dJdu_old;
    BdG = B*dG;
    dupdG = dup'*dG;
    fac = 1/dupdG;
    diff = dup - BdG;
    dupfac=dup*fac;
    diffdup = diff*(dupfac'); 
    B = B + diffdup + diffdup' - (diff'*dG)*(dupfac*dupfac');
  end


    %>>>>>>>>>>>>>>>>>>>>>>>     SELECT BEST MINIMUM     <<<<<<<<<<<<<<<<<<<<<<<<<
  if tr==1,
    Jmin_old = J;
    upmin = up;
  else
    if J<Jmin_old,
      upmin = up;
    end
  end
end

x1(1:Ni1-1)=x1(2:Ni1);           % State 1 to Nu = actual controls
if upmin(1)>1 | upmin(1)<-1
   upmin(1)=upmin(1);
end
if Normalize1
   upmin=(upmin+1)*(maxp1-minp1)/2+minp1;
end
x1(Ni1:Nu1+Ni1-1)=upmin;           % State 1 to Nu = actual controls
for k=1:Nj1-2
  x1(Nu1+Ni1+1+(k-1)*(S11+1):Nu1+Ni1+S11+(k-1)*(S11+1))=x1(Nu1+Ni1+1+(k)*(S11+1):Nu1+Ni1+S11+(k)*(S11+1));
  x1(Nu1+Ni1+(k-1)*(S11+1))=x1(Nu1+Ni1+(k)*(S11+1));                             % delayed plant output
end
if Nj1>=2
   if Normalize1
      x1(Nu1+Ni1+(Nj1-2)*(S11+1))=((u1(3)-mint1)*2/(maxt1-mint1)-1);            % state Nu+1 = NN output
   else
      x1(Nu1+Ni1+(Nj1-2)*(S11+1))=u1(2);
   end
   x1(Nu1+Ni1+1+(Nj1-2)*(S11+1):Nu1+Ni1+S11+(Nj1-2)*(S11+1))=Ai{1,Nj1};    % State Nu+2... = delayed layer 1 output.
end
 

sys1=x1;


%end mdlUpdate


x1(Nu1+Ni1+1+(Nj1-2)*(S11+1):Nu1+Ni1+S11+(Nj1-2)*(S11+1))=Ai{1,Nj1};    % State Nu+2... = delayed layer 1 output.
sys1=x1;


%end mdlUpdate


