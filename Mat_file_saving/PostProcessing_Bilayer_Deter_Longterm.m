%% Post Processing
clear;clc;

deOr = 10;
deOrT = 2.1;

A11 = 1; A12 = 0.1;
A21 = 1; A22 = 0.1;
A1 = 1; A2 = 0.2;
% A11 = 1; A12 = 0.1;
% A21 = -11.5625; A22 = 0;
% A1 = 1.125; A2 = 0.16875 ;
L = 15;
T = 0.0025;
mu = 1;
gamma1 = 1;
gamma2 = 1;
N = 2^deOr;
DeltaX = L/2^deOr;
DeltaT = DeltaX^deOrT;
x = linspace(0,L,N+1);
q = 2*pi./(N*DeltaX) .* [1:N/2];
%%
figure(1)
timeframe = [0, 30, 10000, 11400, 22800, 24700];
for i = 1:6
    h_post = eval(['H_storage_want',num2str(i)]);
    h1 = h_post(:,3:end/2-3)'; % detect lower layer
    h2 = h_post(:,end/2+3:end-3)'; % detect upper layer
    h = h2 - h1;
    
    ax(i) = subaxis(2,3,i,'SpacingVert',0,'SpacingHoriz',0.02,'MR',0.02,'ML',0.05); 
    %ax(i) = subplot(3,3,i);
    if i == 1
        area(x(1:end-1),[h1(:,1), h(:,1)]);
    else
        area(x(1:end-1),[h1(:,end), h(:,end)]);
    end
    if any(i == 1:3)
        set(ax(i),'XTickLabel',[])
    else
        xlabel('x(-)')
    end
    if all(i ~= [1,4])
        set(ax(i),'YTickLabel',[])
    else
        ylabel('h(-)')
    end
    axis equal
    axis([0, L, 0, 7])
    %xlabel('x (-)'); ylabel('h (-)')
    title(['t = ', num2str(round(timeframe(i),1))])
end
linkaxes(ax,'y')

%%
yyaxis left
h_post = eval(['H_storage_want',num2str(6)]);
    h1 = h_post(:,3:end/2-3)'; % detect lower layer
    h2 = h_post(:,end/2+3:end-3)'; % detect upper layer
    h = h2 - h1;
    area(x(1:end-1),[h1(:,end), h(:,end)]);
   
for i = 1:711
    h_post = eval(['H_storage_want',num2str(i)]);
    h1 = h_post(:,3:end/2-3)'; % detect lower layer
    h2 = h_post(:,end/2+3:end-3)'; % detect upper layer
    h = h2 - h1;
    area(x(1:end-1),[h1(:,end), h(:,end)]);
    drawnow
end

% modify the file
filename = ['Bi.Deter.Longterm','_N_',num2str(deOr),'_OrderT_',num2str(deOrT),'_L_',num2str(L),...^,'_T_',num2str(T),...
    '_A11_',num2str(A11),'_A12_',num2str(A12),'_A21_',num2str(A21),...
    '_A22_',num2str(A22),'_A1_',num2str(A1),'_A2_',num2str(A2),...
    '_mu_',num2str(mu),'_gamma1_',num2str(gamma1),'_gamma2_',num2str(gamma2),'.mat']; %,'_test',num2str(1),'.mat'];
%S = load(filename);
%fields = fieldnames(S);

count = 1;
filename_new = ['Bi.Deter.Longterm.Modified','_N_',num2str(deOr),'_OrderT_',num2str(deOrT),'_L_',num2str(L),'_',num2str(count),'.mat'];
%% modify the missing data at the first row
for i = 1:100
    S = load(filename,['H_storage',num2str(i)]);
    S_next = load(filename,['H_storage',num2str(i+1)]);
    fields = fieldnames(S);
    fields_next = fieldnames(S_next);
    if i == 1
        temp = S.(fields{1});
        assignin('base', ['H_storage_new',num2str(i)], temp );
    else
        h_post = S.(fields{1});
        h_post_next = S_next.(fields_next{1});
        temp = [h_post(i:1000,:); h_post_next(2:i,:)];
        assignin('base', ['H_storage_new',num2str(i)], temp );
    end
    
    if mod(i,100) == 1
        save(filename_new,['H_storage_new',num2str(i)])
    else
        save(filename_new,['H_storage_new',num2str(i)],'-append')
    end
    
    if mod(i,100) == 0
        count = count + 1;
        filename_new = ['Bi.Deter.Longterm.Modified','_N_',num2str(deOr),'_OrderT_',num2str(deOrT),'_L_',num2str(L),'_',num2str(count),'.mat'];
    end
    clear(['H_storage_new',num2str(i)])
end

%
% for i = 1:109
%     if i == 109
%         assignin('base', ['H_storage',num2str(i)], H_storage((i-1)*1000+1:108891,:));
%     else
%         assignin('base', ['H_storage',num2str(i)], H_storage((i-1)*1000+1:i*1000,:));
%     end
%     if i == 1
%         save(filename,['H_storage',num2str(i)])
%     else
%         save(filename,['H_storage',num2str(i)],'-append')
%     end
% end

%% evolution inspection

count = 1;
filename_new = ['Bi.Deter.Longterm.Modified','_N_',num2str(deOr),'_OrderT_',num2str(deOrT),'_L_',num2str(L),'_',num2str(count),'.mat'];

for i = 1:712
    if i < 
    if mod(i,100) == 1
        S = load(filename_new);
        fields = fieldnames(S);
    end
    if mod(i,100) == 0
        h_post = S.(fields{100});
    else
        h_post = S.(fields{mod(i,100)});
    end
    
    %     for j = 1:floor(size(h_post,1)/50)
    %         area(linspace(0,L,N),[h_post(50*j,3:N+2)',(h_post(50*j,N+7:2*N+6)-h_post(50*j,3:N+2))']);
    %         text(L/2,5,['t = ',num2str((i-1)*100*1000*DeltaT +5000*j*DeltaT)])
    %         xlabel('x (-)')
    %         ylabel('h (-)')
    %         drawnow
    %     end
    h1 = h_post(end,1:end/2);
    h2 = h_post(end,end/2+1:end);
    h = h2 - h1;
    
    dE1 = 1/2*( ( h1(5:end-1)- h1(3:end-3) )/2/DeltaX ).^2;
    dE2 = gamma2/2*( ( h2(5:end-1)- h2(3:end-3) )/2/DeltaX ).^2;
    dE3 = potential(h1(4:end-2),A21,A22);
    dE4 = potential(h2(4:end-2),A11,A12);
    dE5 = potential(h(4:end-2),A1,A2);
    E1(i) = trapz(x(2:end),dE1');
    E2(i) = trapz(x(2:end),dE2');
    E3(i) = trapz(x(2:end),dE3');
    E4(i) = trapz(x(2:end),dE4');
    E5(i) = trapz(x(2:end),dE5');
    %E = trapz(x(2:end),(dE1+dE2+dE3+dE4+dE5)');
    
    
%     P1 = -(h1(:,5:end-1)-2*h1(:,4:end-2)+h1(:,3:end-3))/DeltaX^2 - gamma2*(h2(:,5:end-1)-2*h2(:,4:end-2)+h2(:,3:end-3))/DeltaX^2 +...
%         disjoiningP(h1(:,4:end-2),A11,A12,0) + disjoiningP(h2(:,4:end-2),A21,A22,0);
%     P2 =  - gamma2*(h2(:,5:end-1)-2*h2(:,4:end-2)+h2(:,3:end-3))/DeltaX^2 +...
%         disjoiningP(h(:,4:end-2),A1,A2,0) + disjoiningP(h2(:,4:end-2),A21,A22,0);
%     [pks, locs] = findpeaks(h1(end,4:end-2),'MinPeakHeight',2);
%     diffP.P1(i) = P1(end,locs(1))-P1(end,locs(2));
%     diffP.time(i) = 100*DeltaT*1000*i
    if mod(i,100) == 0
        clear('S')
        count = count + 1;
        filename_new = ['Bi.Deter.Longterm.Modified','_N_',num2str(deOr),'_OrderT_',num2str(deOrT),'_L_',num2str(L),'_',num2str(count),'.mat'];
    end
    
end

%% mean pressure
subplot(2,1,1)
area(x(1:end-1),[h_post(end,3:N+2)',(h_post(end,N+7:2*N+6)-h_post(end,3:N+2))']);
subplot(2,1,2)
boundedline(x(1:end-1),mean(P1,1),std(P1,1),'b','alpha',x(1:end-1),mean(P2,1),std(P2,1),'r','alpha')

yyaxis left
area(x(1:end-1),[h_post(end,3:N+2)',(h_post(end,N+7:2*N+6)-h_post(end,3:N+2))']);
yyaxis right
boundedline(x(1:end-1),mean(P1,1),std(P1,1),'b','alpha',x(1:end-1),mean(P2,1),std(P2,1),'r','alpha')


function phi = potential(h,A1,A2)
    phi = -A1./2./h.^2 + A2./3./h.^3;
end