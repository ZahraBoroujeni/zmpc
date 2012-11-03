ctrl = {'lqru', 'pid', 'mpc'};

%controller_type = 'lqru';
%controller_type = 'pid'; 
%controller_type = 'mpc'; 


for i = 3
    
controller_type = ctrl{i};

tic;
clear vehicle;
qzVehicle;
toc;
cvx_solver sedumi;
cvx_quiet(true); 



Tmax = 3;
t = 0:vehicle.dt:Tmax;
n = length(t);
U = zeros(4,n);
X = zeros(2,n);
FM = zeros(3,n);
Mycmd = zeros(1,n);
thetaCmd = Mycmd;
%estimator
Mydist = Mycmd;
Xest = X;

%target Fx Fz
%Ftarget = [15; -1.2*vehicle.weight];
Ftarget = [3; -1*vehicle.weight];

%initial attitude
vehicle.x = [5*pi/180 ; -1*pi/180];
%intialize disturbance estimator
vehicle.estimator_dist.xd = [vehicle.x ; 0];

X(:,1) = vehicle.x;
FM(:,1) = vehicle.FM;
Xest(:,1) = X(:,1);


for i = 1:n-1
    switch controller_type
        case 'mpc' 
            vehicle = control_cvx(vehicle,Ftarget);
            Mycmd(:,i) = 0;
            thetaCmd(:,i) = 0;
            
        case 'pid'
            vehicle = control_pid(vehicle,Ftarget);
            Mycmd(:,i) = vehicle.control_pid.Mycmd;
            thetaCmd(:,i) = vehicle.control_pid.theta_cmd;
            
        case 'lqru'
            vehicle = control_lqru(vehicle,Ftarget);
            Mycmd(:,i) = vehicle.control_lqru.Mycmd;
            thetaCmd(:,i) = vehicle.control_lqru.theta_cmd;
    end
    U(:,i) = vehicle.U;
    vehicle = dynamics(vehicle,U(:,i));
    X(:,i+1) = vehicle.x; 
    FM(:,i+1) = vehicle.FM;
    
    vehicle = estimator_dist(vehicle); 
    Xest(:,i+1) = vehicle.estimator_dist.xd(1:2);
    Mydist(:,i) = vehicle.estimator_dist.Myd;
   
    
    if(mod(i,50) == 0)
        fprintf(1,'%2.1f%%\n',100*i/n);
    end
end
%%
toc;
U(:,end) = U(:,end-1);
U = U/(vehicle.weight/size(U,1));
thetaCmd_d = thetaCmd * 180/pi;
theta_d = X(1,:)*180/pi;
thetaEst_d = Xest(1,:)*180/pi;
Fx = FM(1,:);
Fz = FM(2,:);
My = FM(3,:);

figuren(controller_type); clf;
subplot(2,2,1);
plot(t,theta_d,'b',t,thetaCmd_d,'b--',t,thetaEst_d,'k--');
ylim([-20 20]); grid on;
xlabel('time'); ylabel('theta');
legend('Theta','ThetaCmd','ThetaEst');


subplot(2,2,2);
plot(t,U(1,:), t, U(2,:) , t, U(3,:), t, U(4,:));
ylim([0 1.3]); grid on;
xlabel('time'); ylabel('Normalized control');
%legend({'t1','t2','t3','t4'});

subplot(2,2,3);
plot(t,Fx/Ftarget(1) , t, Fz/Ftarget(2));
ylim([-1 1.5]); grid on;
xlabel('time'); ylabel('Normalized Forces');
legend({'fx', 'fz'});

subplot(2,2,4);
plot(t,My,'b', t, Mycmd,'b--', t, -Mydist, 'k--'); 
ylim([-3 3]); grid on;
xlabel('time'); ylabel('My');
legend('My','Mycmd','-Mydist');
toc;


end