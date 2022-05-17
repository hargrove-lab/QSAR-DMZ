function [R2]=Rsquare(y,yhat)

R2=1-sum((y-yhat).*(y-yhat))/sum((y-mean(y)).*(y-mean(y)));
end
