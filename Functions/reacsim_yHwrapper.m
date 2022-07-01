function yH = reacsim_yHwrapper(x, reacstruc)

%     reacstruc.process.lambda0 = x(1);
%     reacstruc.process.pH      = x(2);
%     reacstruc.process.T       = x(3);
%     reacstruc.process.Co      = x(4);
% 
%     reacstruc = reacsim(reacstruc);
%     
%     yH = reacstruc.out.y(end,10);

    for i = 1:size(x,1)
%         reacstruc.process.lambda0 = x(i,1);
%         reacstruc.process.pH      = x(i,2);
%         reacstruc.process.T       = x(i,3);
%         reacstruc.process.Co      = x(i,4);
%         reacstruc.process.tdose   = x(i,5);
        reacstruc.process.T       = x(i,1);
        reacstruc.process.pH      = x(i,2);
        reacstruc.process.Co      = x(i,3);
        reacstruc.process.lambda0 = x(i,4);
        reacstruc.process.tdose   = x(i,5);

        reacstruc = reacsim(reacstruc);
        
        yH(i) = reacstruc.out.y(end,10);
    end
end