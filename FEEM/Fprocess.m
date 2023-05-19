function [F,em,ex]=Fprocess(data,fraction,exmin,exmax,flag)

[N,M]=size(data); c=0; if flag==1; data(data==0)=NaN; end

for i=2:1:M
    c=c+1; MIN=min(data(2:N,i)); F(:,c)=data(2:N,i)-min(data(2:N,i));
end

F=F'; ex=data(1,2:M); em=data(2:N,1); [N,M]=size(F);



for i=1:N
    for j=1:M
        EM=em(j);
        EX=ex(i);
        if EM>=EX*(1-fraction)
            if EM<=EX*(1+fraction)
                F(i,j)=NaN;
            end
        end
    end
end
  
for i=1:N
    for j=1:M
        EM=em(j);
        EX=ex(i);
        if EM>=(2*EX)*(1-fraction)
            if EM<=(2*EX)*(1+fraction)
                F(i,j)=NaN;
            end
        end
    end
end

% now modify for ex max and min

c=0;
for k=1:length(ex);
    extst=ex(k);
    if extst>exmin
        if extst<exmax
        c=c+1;
        EXX(c)=extst;
        FX(c,:)=F(k,:);
        end
    end
end

F=FX; ex=EXX;

end