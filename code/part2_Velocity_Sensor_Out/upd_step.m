function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the  state
    %uEst - estimated mean of the state
    
    % Initialising values for n, vt, R, alpha, k, beta
    R =  0.00001* eye(3); % System Noise
    vt = normrnd(0,(0.0005),[3,1]); 
    n = 15; % n = size of state vector
    alpha = 0.001;
    k = 1;
    beta = 2;
    sqrtCovarEst = chol(covarEst,"lower");

    % Computing spread of sigma points:
    lambda = ((alpha^2)*(n + k)) - n;

    % Next we propogate the sigma points
    X_t = zeros(15, ((2*n)+1));  % We create a 15X31 zero matrix for plugging in vals

    % Setting First Coloumn of the X_t (27X31).
    X_t(:,1) = uEst;
    
    
    % Find all 2n+1 = 31 sigma points
    for i = 1:n
        X_t(:, i+1) = uEst + (sqrt(n+lambda) * sqrtCovarEst(:,i));
        X_t(:, i+n+1) = uEst - (sqrt(n+lambda) * sqrtCovarEst(:, i));
    end


    % % Rotation (IMU to camera frame)
    Rb2c = transpose([0.707 -0.707 0; -0.707 -0.707 0; 0 0 -1]); % eul2rotm([-pi/4,pi,0])
    % % Translation (IMU to camera frame)
    % Tc2b = [0.0283; -0.0283; -0.0300]; % eul2rotm([-pi/4,pi,0]) * [-0.04, 0.0, -0.03]';
    % Skewing the obtained Translation matrix
    SkwTc2b = [      0    0.0300   -0.0283;
                -0.0300         0   -0.0283;
                0.0283    0.0283         0];
    % Rotation matrix to go from camera to IMU frame
    Rc2b = transpose(Rb2c);
    % Rotation matrix to go from IMU to world frame
    Rb2w=[cos(uEst(6,1))*cos(uEst(5,1)), cos(uEst(6,1))*sin(uEst(4,1))*sin(uEst(5,1)) - cos(uEst(4,1))*sin(uEst(6,1)), sin(uEst(6,1))*sin(uEst(4,1)) + cos(uEst(6,1))*cos(uEst(4,1))*sin(uEst(5,1));
             cos(uEst(5,1))*sin(uEst(6,1)), cos(uEst(6,1))*cos(uEst(4,1)) + sin(uEst(6,1))*sin(uEst(4,1))*sin(uEst(5,1)), cos(uEst(4,1))*sin(uEst(6,1))*sin(uEst(5,1)) - cos(uEst(6,1))*sin(uEst(4,1));
                           -sin(uEst(5,1)),                                                cos(uEst(5,1))*sin(uEst(4,1)),                                               cos(uEst(4,1))*cos(uEst(5,1))];
      
    
    
    % Initialise update vector matrix to accomodate all sigma points after
    % being passed through non-linear update model
    Zt = zeros(3, ((2*n)+1));
    % Find all update vectors through iteration
    for i = 1:((2*n)+1)
        x1 = X_t(1:3, i);
        x2 = X_t(4:6, i); 
        x3 = X_t(7:9, i);
        x4 = X_t(10:12, i);
        x5 = X_t(13:15, i);
        l_x2_x3 = Rb2w \ x3;
        % Non-linear update model
        Zt(:,i) = (Rb2c * l_x2_x3) - (Rb2c * SkwTc2b * (Rc2b * z_t(4:6,1))) + vt;
        % Z_t(:,i)=R_c2b*(R_b2w')*x_aug(7:9,i)-R_c2b*skew_b2c*(R_c2b')*z_t(4:6);
    end
   
    
    % Find weights for mean and covariances 
    Wm = zeros(1, ((2*n)+1));
    Wm(1) = lambda / (n + lambda);
    Wc = zeros(1, ((2*n)+1));
    Wc(1) = ((lambda) / (n + lambda)) + (1 - (alpha^2) + beta);
    

   for i = 2:((2*n)+1)
        Wm(i) = 1 / (2 * (n + lambda));
        Wc(i) = 1 / (2 * (n + lambda));
    end
    
    % Initialise estimated measurement update vector
    Zut = zeros(3,1);
    % Find weighted sum of all update vectors obtained by passing all sigma 
    % points through non-linear update model
    
    for i = 1:((2*n)+1)
        Zut = Zut + (Wm(i) * Zt(:,i));
    end
    
    % Initialise covariance and cross-covariance matrices
    Ct = zeros(15,3);
    St = zeros(3,3);
    % Computed weighted sum of covariance and cross-covariance matrices
    for i = 1:((2*n)+1)
        Ct = Ct + (Wc(i) * (X_t(:,i) - uEst) * (Zt(:,i) - Zut)');
        St = St + (Wc(i) * (Zt(:,i) - Zut) * (Zt(:,i) - Zut)');
    end
    St = St + R;

   K = Ct / St;
   uCurr = uEst + K * (z_t(1:3) - Zut);
   covar_curr = covarEst - (K * St * K');
end
