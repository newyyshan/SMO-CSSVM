classdef CSSVM
    %CSSVM Summary of this class goes here
    %Detailed explanation goes here
    
    properties
        alpha=[];   %dual variables
        b=0;        %offset
        pms=[0;0;0];       %two parameters C,C+ and C-
        objective=0;
    end
    methods(Static = true)
        function obj=CSSVM(C,C_plus,C_minus)
            global x;
            global y;
            H=CSSVM.Kernel(x,x);
            H=H.*(y*y');  %Q
            size_training=length(y);
            omega=1*(y==1)+(1/(2*C_minus-1))*(y==-1);
            f=-ones(size_training,1);
            f=f.*omega;
            ub=C*C_plus*(y==1)+C*(2*C_minus-1)*(y==-1);  %set the upper bound
            [alpha,b,objective,~]=CSSVM.quadsmo(H,f,ub,C,C_plus,C_minus,omega);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.pms(1)=C;
            obj.pms(2)=C_plus;
            obj.pms(3)=C_minus;
            obj.alpha=alpha;
            obj.b=b;
            obj.objective=objective;
        end
        function k = Kernel(x, y)%the kernel function
            % function k = kernel(x, y);
            %
            %	x: (Lx,N) with Lx: number of points; N: dimension
            %	y: (Ly,N) with Ly: number of points
            %	k: (Lx,Ly)
            %
            %	KTYPE = 1:      linear kernel:      x*y'
            %	KTYPE = 2,3,4:  polynomial kernel:  (x*y'*KSCALE+1)^KTYPE
            %	KTYPE = 5:      sigmoidal kernel:   tanh(x*y'*KSCALE)
            %	KTYPE = 6:	gaussian kernel with variance 1/(2*KSCALE)
            %
            %       assumes that x and y are in the range [-1:+1]/KSCALE (for KTYPE<6)
            
            global KTYPE
            global KSCALE
            
            k = x*y';
            if KTYPE == 1				% linear
                % take as is
            elseif KTYPE <= 4			% polynomial
                k = (k*KSCALE+1).^KTYPE;
            elseif KTYPE == 5			% sigmoidal
                k = tanh(k*KSCALE);
            elseif KTYPE == 6			% gaussian
                [Lx,~] = size(x);       % the number of x rows
                [Ly,~] = size(y);
                k = 2*k;
                k = k-sum(x.^2,2)*ones(1,Ly);   %sum(A,2) means compute the sum of the elements in each row
                k = k-ones(Lx,1)*sum(y.^2,2)';
                k = exp(k*KSCALE);
            end
        end
        function [local_alpha,b,object,exitflag] = quadsmo(initQIn,initfIn,initUbIn,C,C_plus,C_minus,omega)
            global fake_zero  y  max_iteration_smo epsilon
            max_iteration_smo=100000;
            exitflag=1;
            local_length=length(initfIn);%n
            local_alpha=CSSVM.InitialValue(C,C_plus,C_minus);%n*1
            gi=CSSVM.CaculateF2(local_alpha,initQIn,omega);%n*1
            index_iteration=0;
            while index_iteration < max_iteration_smo
                [local_index,min_value1,~,max_value1,~] = CSSVM.RandomSelSmo2(local_alpha,gi,initUbIn);
                if max_value1-min_value1<=epsilon   %stopping condition
                    b=(max_value1+min_value1)/2;
                    break;
                end
                if local_index(1)==local_index(2)   %the max and the min at the same point
                    break;
                end
                other_local_index=(1:local_length);%1*n
                other_local_index(local_index)=[];  %1*(n-2)
                initQ=initQIn(local_index,local_index); %[Qii Qij;Qji Qjj]
                initf=initfIn(local_index);%2*1
                otherQ=initQIn(local_index,other_local_index);  %2*£¨n-2£©
                local_tmp3=otherQ*local_alpha(other_local_index);   %[2*£¨n-2£©]*[(n-2)*1]
                initf=initf+local_tmp3; %(-e_B+Q_BN*alpha_n_k)
                sum_two_alpha=sum(local_alpha(local_index).*y(local_index)); %y1alpha1+y2alpha=zeta(constant)
                fval0=0.5*local_alpha(local_index)'*initQ*local_alpha(local_index)+initf'*local_alpha(local_index);%the old value of two original alpha
                [initAlpha] = CSSVM.OneSmo(initQ,initf,sum_two_alpha,local_index,initUbIn);
                fval1=0.5*initAlpha'*initQ*initAlpha+initf'*initAlpha;  %the new value of the two updated alpha
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if fval0-fval1<-fake_zero   %error
                    error('myquadprog do not convergence!! \n')
                end
                fprintf('%e\n ',fval0-fval1);
                local_alpha(local_index)=initAlpha; %update the two new alpha
                gi=CSSVM.CaculateF2(local_alpha,initQIn,omega);
                index_iteration = index_iteration+1;
            end
            fprintf('%d\n',index_iteration);
            if index_iteration >= max_iteration_smo  
                exitflag=0;
            end
            object=0.5*local_alpha'*initQIn*local_alpha+initfIn'*local_alpha;  %objective value
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%check the KKT conditions
            local_obj.alpha=local_alpha;
            local_obj.b=b;
            outflag=CSSVM.TestKKT(local_obj,C,C_plus,C_minus,omega);
            if outflag==1
                fprintf('meet the KKT conditions\n');
            else
                fprintf('don not meet the KKT conditions\n');
            end
        end
        function intial_alpha=InitialValue(C,C_plus,C_minus)%initialize alpha
            global y;
            length_plus=sum(y==1);
            length_minus=sum(y==-1);
            intial_alpha=zeros(length(y),1);
            if C*C_plus*length_plus > length_minus*C*(2*C_minus-1)
                intial_alpha(y==-1)= length_minus*C*(2*C_minus-1)/length_minus;
                intial_alpha(y==1) = length_minus*C*(2*C_minus-1)/length_plus;
            else
                intial_alpha(y==-1)=C*C_plus*length_plus/length_minus;
                intial_alpha(y==1)=C*C_plus*length_plus/length_plus;
            end
        end
        function f=CaculateF2(local_alpha,initQIn,omega)%calculate gi
            global y;
            f=initQIn*local_alpha-omega;    % f=initQIn*local_alpha-omega omega is length(initQIn) by 1
            f=-f.*y;
        end
        function [res,min_value1,min_index,max_value1,max_index] = RandomSelSmo2(local_alpha,gi,initUbIn)   %select the two alpha
            % min max find two samples
            global  fake_zero y
            local_length=length(local_alpha);
            index=(1:local_length);
            flag_up=((local_alpha>fake_zero) & y==-1) | ((local_alpha<initUbIn -fake_zero) & y==1);%I_up
            index_up=index(flag_up);    
            flag_low=((local_alpha>fake_zero) & y==1) | ((local_alpha<initUbIn -fake_zero) & y==-1);%I_low
            index_low=index(flag_low);
            [min_value1,min_index]=min(gi(flag_low));
            [max_value1,max_index]=max(gi(flag_up));
            res=[index_low(min_index);index_up(max_index)]; %min & max index
        end
        function [initAlpha] = OneSmo(initQ,initf,sum_two_alpha,two_index,ub) %update the two alpha
            global y
            first_index=two_index(1);
            second_index=two_index(2);
            y1=y(first_index);
            y2=y(second_index);
            yy=y1*y2;
            %convert to quadratic equation of one variable
            a=(initQ(1,1)+initQ(2,2)-2*yy*initQ(1,2))/2;
            b=-sum_two_alpha*y1*initQ(2,2)+initQ(1,2)*y2*sum_two_alpha+initf'*[1;-yy];
            % c=0.5*initQ(2,2)*sum_two_alpha*sum_two_alpha + initf(2)*sum_two_alpha*y2;
            if yy==1    %y1=y2  y1*alpha1+y2*alpha2=zeta --y1(y1*alpha1+y2*alpha)=alpha1+alpha2
                left_bound=max(0,y1*sum_two_alpha-ub(second_index));     %L
                right_bound=min(ub(first_index),y1*sum_two_alpha);      %H
            else    %y1<>y2  y1(y1*alpha1+y2*alpha)=alpha1-alpha2
                left_bound=max(0,y1*sum_two_alpha);
                right_bound=min(ub(first_index),ub(second_index)+y1*sum_two_alpha);
            end
            opitimal_alpha1=-b/(2*a);
            if opitimal_alpha1 > right_bound
                opitimal_alpha1=right_bound;
            end
            if opitimal_alpha1 < left_bound
                opitimal_alpha1=left_bound;
            end
            initAlpha=[opitimal_alpha1; y2*sum_two_alpha-yy*opitimal_alpha1];
        end
        function out_KKT=TestKKT(os,C,C_plus,C_minus,omega)
            global fake_zero
            global y
            out_KKT=0;
            local_g=CSSVM.CaculateF(os,omega);
            alpha=os.alpha;
            zero=alpha'*y;
            if zero<fake_zero && zero>-fake_zero
                local_tmp=CSSVM.SubKKT(alpha,local_g,C,C_plus,C_minus);
                if local_tmp==1
                    out_KKT=1;
                end
            end
        end
        function f=CaculateF(os,omega)
            global x;
            global y;
            Q=CSSVM.Kernel(x,x).*(y*y');
            local_Q=[y Q];
            local_alpha=[os.b;os.alpha];
            local_f=local_Q*local_alpha;
            f=local_f-omega;
        end
        function out_KKT=SubKKT(alpha,local_g,C,C_plus,C_minus)
            global fake_zero
            global y;
            out_KKT=1;
            sum_length=length(alpha);
            C=C*C_plus*(y==1)+C*(2*C_minus-1)*(y==-1);
            for i=1:sum_length
                if alpha(i)<fake_zero
                    if local_g(i)<-fake_zero
                        out_KKT=0;
                        break;
                    end
                else
                    if alpha(i)>C(i)-fake_zero
                        if local_g(i)>fake_zero
                            out_KKT=0;
                            break;
                        end
                    else
                        if local_g(i)>fake_zero || local_g(i)<-fake_zero
                            out_KKT=0;
                            break;
                        end
                    end
                end
            end
        end
    end
end

