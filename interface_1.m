function interface(animationSkip,q,h,x,L_flat,deltaX,t,het,P_het,wave_dom_lsa,e,N,realization,h_adjusted,t_rupt,c)

trace=1;
first_int = 0;
interface=0;
deltaT = deltaX^c;
rupt_index = t_rupt./ deltaT;
N_nodes_het = round(N/(L_flat/P_het));  % No of nodes in one heterogeneous stripe
node_het = 1:N_nodes_het:N;
iter=1;  % iteration index for dydx
pin = zeros(1,2); % pinning index
% int_l = 0;
% int_r = 0;
% vel_inst = zeros(1,q-rupt_index);
for iter_t = rupt_index:q   %16804:16806 %
    
    Y = h(:,iter_t);
    dYdx(:,iter) = ((Y(4:h_adjusted-1)-Y(2:h_adjusted-3))./(2.*deltaX));
    
        if first_int ==1
            
            int_l =0;  % left interface index
            int_r=0;   % right interface index
            % for next time step initializing back to zero
            
            for iter_x = het_index(1):1:het_index(2)
%                [min_dydt,min_index] = min((dYdx(het_index(1):1:het_index(2),iter)))
%                x(het_index(1)+min_index) 
                if  (dYdx(iter_x+1,iter))>= -0.1 && int_l == 0
    
                    fprintf('Position of the %dst interface is %f and the derivative is %f at time %f\n',1, x(iter_x),dYdx(iter_x,iter),t(iter_t))
                    %dYdx(iter_x-10:iter_x+10,iter)
                    
%                     vel_inst(1,iter)= (x(iter_x)-int_left)/deltaT;
                    int_left = x(iter_x);
                    int_l_index = iter;
                    int_l =1;
%                     vel_inst(1,iter)
                    if iter_x == het_index(1) && pin(1) ==0
                        vel(1) = (x(iter_x)-interface)/trace
                        pin(1) = 1;
                        tl_pin = trace;
%                         break
                    
                    end
                                  
                end
               
                
                if (dYdx(iter_x+1,iter))>= 0.1 && (dYdx(iter_x,iter))>= 0 && int_r ==0 && int_l == 1
    
                    fprintf('Position of the %dnd interface is %f and the derivative is %f at time %f\n',2, x(iter_x),dYdx(iter_x,iter),t(iter_t))
                    %dYdx(iter_x-10:iter_x+10,iter)
                    
%                     vel_inst(2,iter)= (x(iter_x)-int_right)/deltaT;
%                     vel_inst(2,iter)
                    int_right =x(iter_x);   
                    int_r = 1;
                    
                    if iter_x == het_index(2) && pin(2) ==0
                        vel(2) = (x(iter_x)-interface)/trace;
                        pin(2) = 1;
                        tr_pin = trace;
                                             
                    end
                    break              
                end
                
                %%direct piece of if statement to just determine velocity
%                 if dYdt(iter_x,iter)>=max_dydt && (iter_x == het_index(1) || iter_x == het_index(2))
%                     dYdt(iter_x,iter)
%                     x(iter_x)
%                     vel = (interface-x(iter_x))/trace
%                     pin = 1;
%                     break
%                                   
%                 end
            end
            
             if int_l == 0 
                    fprintf('left interface not found with current criterion\n')
                    pause(3)
             else if int_r == 0
                  fprintf('right interface not found with current criterion\n')
                    pause(3)   
                 end
             end 
            if pin(2) == 0
                trace = trace + deltaT;
                
            else
                break 
            end
        else
%             [max_dydt,max_index] = max(abs(dYdx(floor(N_nodes_het/4):end,iter)))
%             x(floor(N_nodes_het/4)+max_index)
            for iter_x =floor(N_nodes_het/4):1:(N+1-floor(N_nodes_het/4)) % looking only in the middle of the film
              
                if (dYdx(iter_x,iter))<=0 && (dYdx(iter_x+1,iter))>= 0 %Condition: interface is when the derivative sign changes
                    dYdx(iter_x-10:iter_x+10,iter);
                    dYdx(iter_x,iter);
                    interface = x(iter_x)
                    first_int =1;    % indicating that interface has been locked
                    trace = deltaT;       %To trace the interations till pinning
                    int_right =interface;
                    int_left = interface;
                    if interface >= P_het
                    
                        rem_het = mod(interface,P_het); % remainder to determine how far it is from a wettable patch 
                        rem_index = mod(iter_x-1,N_nodes_het); % -1 to account for the extra 0 element in x
                        het_x = [(interface-rem_het-(0.25*P_het)),(interface-rem_het+(0.25.*P_het))]; % previous and next patch
                        het_index = [(iter_x-rem_index-(N_nodes_het*0.25)),(iter_x-rem_index+(0.25*N_nodes_het))]; % previous and next patch 
                    %(looks complicated as it is a cosine curve
                    else
                        rem_het = P_het - interface; % remainder to determine how far it is from a wettable patch 
                        rem_index = N_nodes_het - iter_x +1; % -1 to account for the extra 0 element in x
                        het_x = [(interface+rem_het-(0.25*P_het)),(interface+rem_het+(0.25.*P_het))]; % previous and next patch
                        het_index = [(iter_x+rem_index-(N_nodes_het*0.25)),(iter_x+rem_index+(0.25*N_nodes_het))]; % previous and next patch 
                    %(looks complicated as it is a cosine curve
                    end
                    
                    break
                end
            end
        end
%         iter = iter +1;
        
        
%         plot(x,dYdt(:,iter))
%         hold on
        
        %[maximum,max_index]=max(abs(dYdt(:,iter)))
        %face(1,all(face==0))
%         dYdt(556,iter)
%         x(556)
        
    
    iter = iter+1;
end
% vel_inst(1,int_l_index-10:int_l_index)
% vel_inst(2,iter-10:iter)
% plot(t(rupt_index+1:iter_t),vel_inst(1,2:iter)) % first element is neglected as velocity is not calculated at t_rupt
% hold on
% plot(t(rupt_index+1:iter_t),vel_inst(2,2:iter))

fprintf('Pinning velocity is %f (%f s) for left side \nPinning velocity is %f (%f s) for the right side \n',vel(1), tl_pin, vel(2), tr_pin)
% vel
% tr_pin
% tl_pin
end