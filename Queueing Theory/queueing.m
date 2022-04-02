clc;
clear all;
lamda = 10:10:90; 
mu = 100;  
entity = 1000000;
for kk =1:length(lamda)
    arr_time = zeros(entity,1);
    fin_time = zeros(entity,1);
    r = rand(entity,1);
    iat = -1/lamda(kk) * log(1-r);  
    arr_time(1) = iat(1);    
    for i=2:entity
       arr_time(i) = arr_time(i-1) + iat(i);  
    end
    r = rand(entity,1);
    ser_time = -1/mu * log(1-r);
    fin_time(1) = arr_time(1)+ser_time(1); 
    for i=2:entity
       fin_time(i) = max(arr_time(i)+ser_time(i),  fin_time(i-1)+ser_time(i));
    end
    total_time = fin_time - arr_time;       
    wait_time  = total_time - ser_time;   
    ave_service_time = sum(ser_time)/entity;
    ave_wait_time = sum(wait_time)/entity;
    ave_total_time = sum(total_time)/entity;
    q_length = zeros(ceil(fin_time(entity))+1,1);
    for i=1:entity
      for t=ceil(arr_time(i)):floor(arr_time(i)+wait_time(i)) 
        q_length(t+1) = q_length(t+1) + 1;
      end;
    end;
    if (kk==4)
        snap = q_length;
    end
    int_q_len =zeros(entity,1); 
    for gg = 2:entity
        if arr_time(gg)<fin_time(gg-1)
            int_q_len(gg) = int_q_len(gg)+1;
        else
            if int_q_len(gg-1)>1
                int_q_len(gg) = int_q_len(gg) -1;
            end
        end
    end
    int_frac_time_q_emp(kk) = length(find(int_q_len==0))/(length(int_q_len)) ;       
    avg_ser_time(kk) = ave_service_time;
    avg_wait_time(kk) = ave_wait_time;
    avg_total_time(kk) = ave_total_time;
    avg_q_length(kk) = mean(q_length);
    act_total_time(kk) = 1/(mu - lamda(kk));
    rho(kk) = lamda(kk)/mu; 
    frac_q_empty_actual(kk) = 1-rho(kk);
    act_q_length(kk) = rho(kk)^2/(1-rho(kk));
    n_baar(kk) =  avg_wait_time(kk)*lamda(kk); %by little formula
    up_conf_level(kk) = avg_q_length(kk) + 1.96*var(q_length)/sqrt(length(q_length));
    lw_conf_level(kk) = avg_q_length(kk) - 1.96*var(q_length)/sqrt(length(q_length));
end
figure()

plot(rho,avg_q_length)
hold on
plot(rho, act_q_length)
hold on
plot(rho,n_baar)
% plot(rho, up_conf_level)
%hold on
% plot(rho, lw_conf_level)
legend('Simulated','Theoretical','by little formula')
title('Average Queue Length')
xlabel('Traffic Density')
ylabel('Average Queue Length')

hold off 

figure()

plot(rho,int_frac_time_q_emp)
hold on 
plot(rho,frac_q_empty_actual)
xlabel('Traffic Density')
ylabel('Fraction of time queue empty')
legend('Simulated','Theoretical')
title('Fraction of time Queue is Empty')
hold off
