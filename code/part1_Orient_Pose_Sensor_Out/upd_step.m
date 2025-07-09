function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the  state
    %uEst - estimated mean of the state
    
    C_t=horzcat(eye(6),zeros(6,9));
    C_t_trnsp=C_t';
    R=eye(6)*0.00000001;
    K_t=covarEst*C_t_trnsp/((C_t*covarEst*C_t_trnsp+R));
    uCurr=uEst+K_t*(z_t-C_t*uEst);
    covar_curr=covarEst-K_t*C_t*covarEst;
end

