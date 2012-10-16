

qzVehicle;

Tmax = 5;
t = 0:vehicle.dt:Tmax;
n = length(t);
U = zeros(4,n);
X = zeros(2,n);
FM = zeros(3,n);
Mycmd = zeros(1,n);

%target Fx Fz
Ftarget = [1; -1.3*vehicle.weight];
%initial attitude
vehicle.x = [0 ; 0];

X(:,1) = vehicle.x;
FM(:,1) = vehicle.FM;
for i = 1:n-1    
    vehicle = control_pid(vehicle,Ftarget);
    Mycmd(:,i) = vehicle.control_pid.Mycmd;
    
    U(:,i) = vehicle.U;
    vehicle = dynamics(vehicle,U(:,i));
    X(:,i+1) = vehicle.x; 
    FM(:,i+1) = vehicle.FM;
end
U(:,end) = U(:,end-1);
U = U/vehicle.tmax;
%%
theta_d = X(1,:)*180/pi;
Fx = FM(1,:);
Fz = FM(2,:);
My = FM(3,:);

figure(1); clf;
subplot(2,2,1);
plot(t,theta_d);
grid on;
xlabel('time'); ylabel('theta');


subplot(2,2,2);
plot(t,U(1,:), t, U(2,:) , t, U(3,:), t, U(4,:));
grid on;
xlabel('time'); ylabel('Normalized control');
legend({'t1','t2','t3','t4'});

subplot(2,2,3);
plot(t,Fx/Ftarget(1) , t, Fz/Ftarget(2));
grid on;
xlabel('time'); ylabel('Normalized Forces');
legend({'fx', 'fz'});

subplot(2,2,4);
plot(t,My, t, Mycmd); 
grid on;
xlabel('time'); ylabel('My');
