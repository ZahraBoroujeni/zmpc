ctrl = {'lqru', 'pid', 'mpc', 'mpc2', 'fastmpc'};

%controller_type = 'lqru';
%controller_type = 'pid'; 
%controller_type = 'mpc'; 


for ctrl_idx = [2, 5]
    
controller_type = ctrl{ctrl_idx};


clear vehicle;
qzVehicle;
cvx_solver sedumi;
cvx_quiet(true); 

%target Fx Fz
%Ftarget = [15; -1.2*vehicle.weight];
Ftarget = [3; -1*vehicle.weight];
uLoop = true;
uTarget = 2;




Tmax = 3;
t = 0:vehicle.dt:Tmax;
n = length(t);
U = zeros(4,n);
FM = zeros(3,n);
Y = zeros(size(vehicle.sysSim.c,1),n);
Mycmd = zeros(1,n);
thetaCmd = Mycmd;
%estimator
Mydist = Mycmd;
ThetaQEst = zeros(2,n);
Fhist = zeros(2,n);


%initial attitude
vehicle.xsim(1) = 5*pi/180;
Y(:,1) = vehicle.sysSim.C * vehicle.xsim;

%intialize disturbance estimator
vehicle.estimator_dist.xd = [vehicle.xsim(1); 0 ; 0; 0];
ThetaQEst(:,1) = vehicle.estimator_dist.xd(1:2);
%vehicle.estimator_dist.Myd = 1;
vehicle = dynamics(vehicle,zeros(4,1));

tstart = tic;
for i = 1:n-1
    if(uLoop)
        vehicle = control_u_pid(vehicle,uTarget);
        Ftarget(1) = vehicle.control_u_pid.FxTarget;
    end
    Fhist(:,i) = Ftarget;
    switch controller_type
        case 'fastmpc'
            vehicle = control_fastmpc(vehicle,Ftarget);
            Mycmd(:,i) = 0;
            thetaCmd(:,i) = 0;
            
        case 'mpc2'
            vehicle = control_cvx_noineq(vehicle,Ftarget);
            Mycmd(:,i) = vehicle.lcm(3,:)  * vehicle.U;
            thetaCmd(:,i) = 0;
            
        case 'mpc' 
            vehicle = control_cvx(vehicle,Ftarget);
            Mycmd(:,i) = 0;
            thetaCmd(:,i) = 0;
            
        case 'pid'
            vehicle = control_pid(vehicle,Ftarget);
            Mycmd(:,i) = vehicle.control_pid.Mycmd;
            thetaCmd(:,i) = vehicle.control_pid.cmd;
            
        case 'lqru'
            vehicle = control_lqru(vehicle,Ftarget);
            Mycmd(:,i) = vehicle.control_lqru.Mycmd;
            thetaCmd(:,i) = vehicle.control_lqru.theta_cmd;
    end
    U(:,i) = vehicle.U;
    vehicle = dynamics(vehicle,U(:,i));
    vehicle = estimator(vehicle); 

    Y(:,i+1) = vehicle.ysim; 
    FM(:,i+1) = vehicle.FM;
    
    ThetaQEst(:,i+1) = vehicle.estimator_dist.xd(1:2);
    Mydist(:,i) = vehicle.estimator_dist.Myd;
   
    
    if(mod(i,50) == 0)
        fprintf(1,'%2.1f%%\n',100*i/n);
    end
end

deltaT = toc(tstart);

if(isfield(vehicle.control_cvx,'numSol'))
    numSol = vehicle.control_cvx.numSol;
    deltaT_persolve =   vehicle.control_cvx.deltaT/vehicle.control_cvx.numSol;
else
    numSol = 0;
    deltaT_persolve = 0;
end
fprintf(1,'runttime was %f (total solves = %i , per solve = %3.0fms) \n', deltaT, numSol, deltaT_persolve*1000);

%%
U = U/(vehicle.weight/size(U,1));
thetaCmd_d = thetaCmd * 180/pi;
theta_d = Y(vehicle.thetaIndex,:)*180/pi;
q_d = Y(vehicle.qIndex,:)*180/pi;
u = Y(vehicle.uIndex,:);
thetaEst_d = ThetaQEst(1,:)*180/pi;
Fx = FM(1,:);
Fz = FM(2,:);
My = FM(3,:);

figuren(controller_type); clf;
subplot(2,2,1);
plot(t,theta_d,'b',t,thetaCmd_d,'b--',t,thetaEst_d,'k--');
hold on;
plot(t,q_d,'g--');
plot(t,u,'r');
ylim([-20 20]); grid on;
xlabel('time'); ylabel('theta');
legend('Theta','ThetaCmd','ThetaEst','q','u','Location','EastOutside');


subplot(2,2,2);
plot(t,U(1,:), t, U(2,:) , t, U(3,:), t, U(4,:));
ylim([0 1.3]); grid on;
xlabel('time'); ylabel('Normalized control');
%legend({'t1','t2','t3','t4'});

subplot(2,2,3);
plot(t,Fx./vehicle.weight , t,Fhist(1,:)./vehicle.weight, t, Fz./Ftarget(2));
ylim([-1 1.5]); grid on;
xlabel('time'); ylabel('Normalized Forces');
legend({'fx', 'fxcmd','fz'});

subplot(2,2,4);
plot(t,My,'b', t, Mycmd,'b--', t, -Mydist, 'k--'); 
ylim([-3 3]); grid on;
xlabel('time'); ylabel('My');
legend('My','Mycmd','-Mydist');


end