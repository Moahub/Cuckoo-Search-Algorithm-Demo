% Cuckoo Search Algorithm (CSA) for demo function:
% f(x)=x1^2-x1*x2+x2^2+2*X1+4*x2+3 minmiztion into interval of [-5,5]
% this code consists of 5 phases according to the cuckoo algorithm concept

% Phase 1 : Constants
D = 2;      %Dimension of the problem
lb= [-10,-10]; %Lower boundary of variables
ub= [10,10];  %Upper boundary of variables
N = 20;     %Population size
n = N;
pa= 0.25;   %Discovery rate of alien eggs/solutions
max_iter=100;
beta=1.5;
%--------------------------------------------------------------------------
% Note Phase 2 is the objective function generation
%--------------------------------------------------------------------------

% Phase 3 : Generate initial population randomly
nest=zeros(N,D);
for i=1:N
    for j=1:D
    nest(i,j)=lb(:,j)+rand.*(ub(:,j)-lb(:,j));
    end
end

fx = Obj(nest);

sigma=(gamma(1+beta)*sin(pi*beta/2)/...
    (gamma((1+beta)/2)*beta*2.^((beta-1)/2))).^(1/beta);
% Phase 4 :Cuckoo Search Main Loop

for iter= 1:max_iter
    [fxmin,ind]=min(fx);
    best= nest(ind,:);
    for j=1:N
        s=nest(j,:);
        X=s;
        %%Levy flights by  Mantegna's algorithm
        u = randn(size(s))*sigma;
        v=randn(size(s));
        step=u./abs(v).^(1/beta);
        Xnew = X+randn(size(s)).*0.01.*step.*(X-best);
        %%check bounds
        for i2=1:D
            if Xnew(i2)>ub(i2)
                Xnew(i2)=ub(i2);
            elseif Xnew(i2)<lb(i2)
                Xnew(i2)=lb(i2);
            end
        end
        %%Perform Greedy Selection
        fxnew=Obj(Xnew);
        if fxnew<fx(j,:)
            nest(j,:)=Xnew;
            fx(j,:)=fxnew;
        end
    end
    [fxmin,k1]=min(fx);
    best= nest(k1,:);
    %%Phase 5: Replace some nests by constructing new solutions/nests
    k=rand(size(nest))<pa;
    stepsizek=rand*(nest(randperm(n),:)-nest(randperm(n),:));
    new_nest=nest+stepsizek.*k;
    for i3=1:size(nest,1)
        s=new_nest(i3,:);
        %%check bounds
        for i2=1:D
            if s(i2)>ub(i2)
                s(i2)=ub(i2);
            elseif s(i2)<lb(i2)
                s(i2)=lb(i2);
            end
        end
        %%Perform Greedy Selection
        fxnew=Obj(s);
        if fxnew<fx(i3,:)
            nest(i3,:)=s;
            fx(i3,:)=fxnew;
        end  
        
    end
    %%Phase 6: Memorize the Best
    [optval, optind]= min(fx);
    BestFx(iter)= optval;
    BestX(iter,:)= nest(optind,:);
    
    %----------------------------------------------------------------------
    [m,h]=meshgrid(-10:0.5:10);
    ne=[m h];
    fxd = m.^2-m.*h+h.^2+2.*m+4.*h+3;
    zmin=min(min(fxd));
    [r,c]=find(fxd==zmin);

    figure(2)
    surf (m,h,fxd)
    hold on
    plot3(nest(optind,1),nest(optind,2),optval,'r*','Linewidth',5)
    hold off
    %----------------------------------------------------------------------

end
%%Phase 7: Plotting the Results
figure
plot(BestFx,'Linewidth',2);
xlabel('Iteration Number');
ylabel('fitness values');
title('Convergence Vs Iterataions')
grid on



        