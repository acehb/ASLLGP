function [X] = NestedLHD(n,m,d)
%Inputs n:the number of runs for the low-fidelity experiment;m:the number of runs for the high-fidelity experiment;d:dimension of the inputs.
%Outputs X:the maximin nested Latin hypercube design
X = nestedsample(n,m,d);

   maxdist = mindist(X);
   %for i=2:1000
   for i=2:500
      x = nestedsample(n,m,d);
      
      newdist = mindist(x);
      if newdist > maxdist
         X = x;
         maxdist = newdist;
      end
   end
end

function x = nestedsample(n,m,d)
pin=zeros(d,n);
tao=zeros(d,m);

for j=1:d
pim=randperm(m);
t=n/m;

for i=1:m  
    tao(j,i)=randi([(pim(i)-1)*t+1,pim(i)*t],1);
end
pin(j,:)=[tao(j,:), setdiff(randperm(n),tao(j,:),'stable') ];
end
x=(pin-rand(d,n))'./n;
end

function s = mindist(x)

   [~,dist] = knnsearch(x,x,'k',2);
   s = min(dist(:,2));
 
end
