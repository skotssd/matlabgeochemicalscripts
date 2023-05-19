function [G,H,I]=resample(data,em,ex,samplefreq)

% try to resample so things look better

[N,M]=size(data); % M is em points N is ex points

for i=1:N
    Fem=data(i,:); c=0;
    for j=1:round(M/samplefreq):M
       c=c+1; Femred(c)=Fem(j); emred(c)=em(j);
    end
    datare(i,:)=Femred;
end

[N,M]=size(datare); % M is em points N is ex points.  sample into 10 steps.

for i=1:M
    Fex=datare(:,i); c=0;
    for j=1:round(N/samplefreq):N
       c=c+1;
       Fexred(c)=Fex(j); exred(c)=ex(j);
    end
    datarere(:,i)=Fexred;
end

G=datarere;
H=emred;
I=exred;

end