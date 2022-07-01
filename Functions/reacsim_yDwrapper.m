function yD = reacsim_yDwrapper(x, reacstruc)
    yD = zeros(size(x,1),1);

    for i = 1:size(x,1)
        reacstruc.process.T       = x(i,1); % C
        reacstruc.process.pH      = x(i,2); %
        reacstruc.process.Co      = x(i,3); % g/L
        reacstruc.process.lambda0 = x(i,4); % mol HEP/mol N9
        reacstruc.process.tdose   = x(i,5);

%         reacstruc.process.lambda0 = x(i,1);
%         reacstruc.process.pH      = x(i,2);
%         reacstruc.process.T       = x(i,3);
%         reacstruc.process.Co      = x(i,4);
%         reacstruc.process.tdose   = x(i,5);

        reacstruc = reacsim(reacstruc);
    
        yD(i) = reacstruc.out.y(end,7);
    end
end