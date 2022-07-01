function [f,reacstruc] = reacoptim(x,reacstruc)
global var counter relcost 

counter=counter+1;
Price = reacstruc.optim.Price;

%Parametre til optim
for i = 1:length(reacstruc.optim.var)
    % Scale back from normalization
    z(i) = x(i).*(reacstruc.optim.UB(i)-reacstruc.optim.LB(i)) + ...
                  reacstruc.optim.LB(i);
              
    % Update model with new z 
    eval([reacstruc.optim.var{i},'=',num2str(z(i),10),';']);
end
z=z(:)';

% Run simulation
reacstruc = reacsim(reacstruc);

% Unwrap solution from structure array
t         = reacstruc.out.t;
y         = reacstruc.out.y;
ysum      = 0*y(:,1:4);             % Allocation
ysum(:,1) = y(:,3);                 % Component A
ysum(:,2) = sum(y(:,4:6),2);        % Monoacylated
ysum(:,3) = sum(y(:,7:9),2);        % Diacylated
ysum(:,4) = y(:,10);                % Triacylated
xsim      = 1-ysum(:,1)-ysum(:,2);

if reacstruc.optim.flag==0
        %clf
        
        f=100;
        %    colvecall={'k','y-','b-','g-','g--','g.-','r-','r--','r.-','m-','k--'};
        colsum={'b-','g-','r-','m-'};
        yfac=[1 1 1 1];
        
        figure(1)
%        subplot(subP(1),subP(2),k)
 %       title(num2str(reacstruc.optim.runvec(k)))
        hold on

        for i=1:4
            plot(t,f*yfac(i)*ysum(:,i),colsum{i})
            ylabel('Content (%, scaled)')
            hold on
        end
%         legend(reacstruc.process.comps{3:6},'Location','NorthEast' )
%         legend boxoff
        hold on
        
       
        %     xlim([0 max(t)])
                 ylim([0 f])
        
        figure(2)
 %       subplot(subP(1),subP(2),k)
%        title(num2str(reacstruc.optim.runvec(k)))
        hold on
        
        
            plot(f*yfac(4)*ysum(:,4),f*yfac(3)*ysum(:,3),colsum{3});
        
            ylabel('Diacylated (%)')
            xlabel('Triacylated (%)')
            hold on            
        
             ylim([0 f])
             xlim([0 .02*yfac(4)*f])
      
    
      ysum(end,1:4)
    
%        'SC reacted:'       
%    reacstruc.out.yreac    
end %if no flag

% f=(1-yield+Price.SCrel*reacstruc.process.lambda0)/yield;    %cost
% Compute objective function
yield     = ysum(end,3);
f = (1 + Price.SCrel*reacstruc.process.lambda0) / yield;    %cost/cost precursor

var(counter,:)   = z;
relcost(counter) = f;
% 
% figure(5)
% hold on
% xlabel('Iterations')
% ylabel('Residual')
% plot(1:counter,relcost,'o')
% ylim([0 2*min(relcost)])
% 
% drawnow
