 %
%  USAGE:  >> run_divertor_model
%
%  PURPOSE: Script to run divertor_model.m to evolve divertor/core dynamics
%		with specified inputs, intended to execute open loop OR 
%		closed loop  evolution. For open loop, inputs waveforms (Gamec(t), 
%		Gamrc(t),Gamrd(t)) are defined explicitly below. For open loop, 
%		input waveforms are determined dynamic by feedback laws defined
%		below. 
%
%  INPUTS:
%
%  OUTPUTS:
%
%  RESTRICTIONS:
%
%  METHOD:  

%  WRITTEN BY:  Dave Humphreys  ON	1/1/00
%
%  MODIFICATION HISTORY:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Explicit Inputs:
  idoopenloop = 1;
  idoclosedloop = 0; 
  tmin = 0;
  tmax = 10.0;    %max time for evolution [s]
  nt = 100;     %# of time points
  dt = (tmax-tmin)/(nt-1);
  Gamec0 = 2e19;  %constant OL e flow rate [electrons/m^3/s]
  Gamrc0 = 4e17;  %constant OL e flow rate [electrons/m^3/s]
  Gamrd0 = 4e18;  %constant OL e flow rate [electrons/m^3/s]
  nec0 = 1e19;    %initial nec
  nrc0 = 1e17;     %initial nrc
  nrd0 = 1e18;    %initial nrd


% Prelims and Constants:
   mu0 = 0.4*pi;

% Derived Values:
  t = linspace(tmin,tmax,nt)';
  Gamect = Gamec0*ones(nt,1);
  Gamrct = Gamrc0*ones(nt,1);
  Gamrdt = Gamrd0*ones(nt,1);
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evolve divertor model Open Loop with input waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if idoopenloop
  nect=zeros(nt,1); nect(1)=nec0;
  nrct=zeros(nt,1); nrct(1)=nrc0;
  nrdt=zeros(nt,1); nrdt(1)=nrd0;
  Wtht = zeros(nt,1);
  Qtardt = zeros(nt,1);
  Pradct = zeros(nt,1);
  Praddt = zeros(nt,1);
  for ii=2:nt
     disp(['Step ',int2str(ii)]);
     [nec_out,Wth,Qtard,nec2,nrc2,nrd2,output_data] = ...
        divertor_model(Gamect(ii),Gamrct(ii),Gamrdt(ii), ...
                       nect(ii-1),nrct(ii-1),nrdt(ii-1),dt);
     nect(ii) = nec2;
     nrct(ii) = nrc2;
     nrdt(ii) = nrd2;
     Wtht(ii) = Wth;
     Qtardt(ii) = Qtard;
     Pradct(ii) = output_data.Pradc;
     Praddt(ii) = output_data.Pradd;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evolve divertor model Closed Loop with feedback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if idoclosedloop
 



end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1),clf,hold off
nrows=3; ncols=1;
subplot(nrows,ncols,1)
  plot(t,nect,'b')
  grid on
  ylabel('Core density','FontSize',14)
title('Plasma state density evolution','FontSize',15)
subplot(nrows,ncols,2)
  plot(t,nrct,'b')
  grid on
  ylabel('Core RI density','FontSize',14)
subplot(nrows,ncols,3)
  plot(t,nrdt,'b')
  grid on
  ylabel('Div RI density','FontSize',14)
xlabel('t [s]','FontSize',15)

figure(2),clf,hold off
nrows=3; ncols=1;
subplot(nrows,ncols,1)
  plot(t,Gamect,'b')
  grid on
  ylabel('Core fueling rate','FontSize',14)
title('Inputs: Fueling & Impurity Flow Rates','FontSize',15)
subplot(nrows,ncols,2)
  plot(t,Gamrct,'b')
  grid on
  ylabel('Core RI flow rate','FontSize',14)
subplot(nrows,ncols,3)
  plot(t,Gamrdt,'b')
  grid on
  ylabel('Div RI flow rate','FontSize',14)
xlabel('t [s]','FontSize',15)

figure(3),clf,hold off
nrows=3; ncols=1;
subplot(nrows,ncols,1)
  plot(t,nect,'b')
  grid on
  ylabel('Core e density','FontSize',14)
title('Outputs: nec,Wth,Qtard','FontSize',15)
subplot(nrows,ncols,2)
  plot(t,Wtht,'b')
  grid on
  ylabel('Core stored energy','FontSize',14)
subplot(nrows,ncols,3)
  plot(t,Qtardt,'b')
  grid on
  ylabel('Div tar ht flux','FontSize',14)
xlabel('t [s]','FontSize',15)

figure(4),clf,hold off
nrows=2; ncols=1;
subplot(nrows,ncols,1)
  plot(t,Pradct,'b')
  grid on
  ylabel('Core radiation','FontSize',14)
title('Radiation Summary: core & div','FontSize',15)
subplot(nrows,ncols,2)
  plot(t,Praddt,'b')
  grid on
  ylabel('Div radiation','FontSize',14)
xlabel('t [s]','FontSize',15)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   disp('All done.')
               
