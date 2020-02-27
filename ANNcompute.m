function [PSOResults, R2, MSENET, Net ] = ANNcompute(X,Y)
Y4=Y;
CONt=1;
Sanitizer=0;
for NNEU1=7:-1:1
    if Sanitizer==1
        break
    end
    for NNEU2=7:-1:0
        if Sanitizer==1
            break
        end
        if NNEU2>0
        Net=newff(X,Y4,{NNEU1,NNEU2});
        disp('2 layers');
        else
        Net=newff(X,Y4,{NNEU1});
        disp('1 layer');
        end
        Net.trainParam.epochs=1000;
        Net.trainParam.goal=1e-8;
        Net.trainParam.min_grad=1e-8;
        %Net.performParam.regularization=0.1;
        Net.trainFcn = 'trainbr';
        %Net.trainParam.max_fail=1000;
        Net.TrainParam.ShowWindow=false;
for i = 1:10
  fprintf('Training ANN %d/%d\n', i, 10);
  Net= train(Net,X,Y4);
  NN{i}=Net;
  YY=sim(Net,X);
  perfs(i) = mse(Net, YY, Y4);
  SSres=sum( (Y4-YY).^2 );
  SStot=sum( (Y4-mean(Y4)).^2 );
  R2(i)=1-SSres/SStot;
end
[MSENET, Ind_Max]=min(perfs);
Net=NN(Ind_Max);
R2=R2(Ind_Max);
CONt=CONt+1;
R2G(1)=R2;
R2G(CONt)=R2;
MSENETA(1)=MSENET;
MSENETA(CONt)=MSENET;
NNN{1}=Net;
NNN(CONt)=Net;
if R2G(CONt) <0.98*R2G(CONt-1)
    Net=NNN(CONt-1);
    R2=R2G(CONt-1);
    MSENET=MSENETA(CONt-1);
    Sanitizer=1;
end
    end
end
    
%ANNsim
disp('Finished');




PSOResults=psoCyclic(Net{1});

end

