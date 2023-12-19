   N=6;M=2;                                                                %number of elements, number of sources
   phi_tr=[50,100]/180*pi;                                                 %true DOA value
   A_tr=zeros(N,M);                                                        %true array matrix
   for m1=1:M                                                              %
       for n1=1:N
           A_tr(n1,m1)=exp(-1i*pi*cos(phi_tr(m1))*(n1-1));
       end
   end
   ES=[1,1,sqrt(2)];                                                       %M+1个相互独立随机变量的标准差
   P_tr=[2,2;2,2];                                                         %true source covariance matrix
   no_tr=[1,2,3,4,2,10];                                                   %true noise variances
   Cy_tr=A_tr*P_tr*(A_tr')+diag(no_tr);                                    %true array covariance matrix
   invCy_tr=inv(Cy_tr);                                                    %inversion of the true array covariance matrix
   
   inv_S=diag(1./sqrt(no_tr));                                             %CRLB
   A_C=inv_S*A_tr;
   De_C=zeros(N,M);                                                        %
   for m6=1:M                                                              %
       for n6=1:N
       De_C(n6,m6)=exp(-1i*pi*cos(phi_tr(m6))*(n6-1))*1i*(n6-1)*pi*sin(phi_tr(m6));
       end
   end
   De_C=inv_S*De_C;
   Pj_D=A_C*inv((A_C')*A_C)*(A_C');
   Pj_C=eye(N)-Pj_D;                                                       %
   invCy_C=diag(sqrt(no_tr))*invCy_tr*diag(sqrt(no_tr));
   
   UU1=P_tr*(A_C')*invCy_C*A_C*P_tr;
   UU2=((De_C')*Pj_C*invCy_C*De_C).';
   UU_M=2*real(((invCy_C*A_C*P_tr).').*((De_C')*Pj_C));
   UW=Pj_D*invCy_C;
   UU_T=inv(real(conj(invCy_C).*invCy_C-conj(UW).*UW));
   
   CRLB_M=inv(2*real(UU1.*UU2)-UU_M*UU_T*(UU_M.'));
   
   KK=10;                                                                  %number of 
   X1=zeros(1,KK);X2=X1;Y1=X1;Y2=X1;X3=X1;Y3=X1;X4=X1;Y4=X1;
   CRLB1=zeros(1,KK);CRLB2=zeros(1,KK);
   number=2000;                                                            %number of independent trials
   
   for ii=1:KK
       T=10*ii;                                                            %快拍数
       C_i=CRLB_M/T;
       CRLB1(ii)=sqrt(C_i(1,1))/pi*180;
       CRLB2(ii)=sqrt(C_i(2,2))/pi*180;                                    %CRLB obtained
       S_tr=zeros(M,T);                                                    %signal matrix
       W_tr=zeros(N,T);                                                    %noise matrix
   for num=1:number                                                        %num-th trial
   
   S_B=normrnd(0,1/sqrt(2),M+1,T)+1i*normrnd(0,1/sqrt(2),M+1,T);           %定义M+1行相互独立的标准循环复序列
   for m=1:M
       S_tr(m,:)=ES(M+1)*S_B(M+1,:);
   end                                                                     %true signal matrix
   for n2=1:N                                                              
       W_tr(n2,:)=sqrt(no_tr(n2)/2)*(normrnd(0,1,1,T)+1i*normrnd(0,1,1,T));%true noise matrix
   end
   Y=A_tr*S_tr+W_tr;                                                       %array signal matrix
   
   S_y=zeros(N,N);                                                         %array covariance matrix estimate
   for t1=1:T
       S_y=S_y+Y(:,t1)*(Y(:,t1)');
   end
   R_y=S_y/T;
   
   %SAGE algorithm
   phi_ite=[45,95]/180*pi;                                                 %initial DOA estimate 
   A_ite=zeros(N,M);                                                       %
   for m4=1:M                                                              %
       for n4=1:N
           A_ite(n4,m4)=exp(-1i*pi*cos(phi_ite(m4))*(n4-1));
       end
   end
   P_ite=ones(1,M);No_ite=ones(1,N);                                       %initial source covariance matrix estimate, initial noise variance estimates
   C_y=diag(No_ite);                                                       %定义用于迭代的样本协方差矩阵
   for m_4=1:M
       C_y=C_y+P_ite(m_4)*A_ite(:,m_4)*(A_ite(:,m_4)');
   end                                                                     %初始化样本协方差矩阵
   invC_y=inv(C_y);                                                        %使用逆矩阵
   k=0;                                                                    %初始化迭代数
   k=k+1;                                                            %第k次迭代开始
   phi_ite1=phi_ite;                                                 %存储上一次
   for m=1:M                                                         %第m个源
       for im=1:M                                                    %更新其他源的功率
           if im~=m
              b_i=invC_y*A_ite(:,im)*P_ite(im);
              gr=real((b_i')*R_y*b_i);
              P_ite(im)=gr+P_ite(im)*max(1-real((A_ite(:,im)')*b_i),0);
           end
       end
       C_m=P_ite(m)*A_ite(:,m)*(A_ite(:,m)')+diag(No_ite);           %当前源的协方差矩阵
       R_R=C_m*invC_y*R_y*invC_y*C_m+C_m-C_m*invC_y*C_m;             %得到当前源的完整数据协方差矩阵 复数矩阵
       inv_N=diag(1./No_ite);                                        %对噪声协方差矩阵求逆
       R_ite=inv_N*R_R*inv_N;                                        %得到当前源的变形完整数据协方差矩阵 复数矩阵

       Sum_g=0;                                                      %计算当前源的比较函数值  梯度下降法开始
       for g1=1:N                                                    %第一部分与角无关，可以重复使用，均匀线阵的函数
           Sum_g=Sum_g+real(R_ite(g1,g1));                           %取实部保证实数
       end
       Sum_G=0;                                                      %第二部分与角有关，每次更新方位角后都要重新计算
       for p_1=1:N-1
           s=0;
           for q_1=1:N-p_1
               s=s+R_ite(q_1+p_1,q_1);
           end
           Sum_G=Sum_G+exp(1i*pi*cos(phi_ite(m))*p_1)*s;             %这里有当前源的方位角
       end
       g_ite=Sum_g+2*real(Sum_G);                                    %得到当前源的比较函数值 实数

       sum_de=0;                                                     %计算当前源比较函数的导数值
       for n3=1:N-1
           s=0;
           for qq=1:N-n3
               s=s+R_ite(qq+n3,qq);
           end
           sum_de=sum_de+n3*exp(1i*pi*cos(phi_ite(m))*n3)*s;
       end
       g_de_ite=2*pi*sin(phi_ite(m))*imag(sum_de);                   %得到当前源比较函数的导数值

       while abs(g_de_ite)>0.001                                     %判断导数值大小
             if g_de_ite>0                                           %初始化步长
                t=(pi-phi_ite(m))/g_de_ite/10;
             else
                t=-phi_ite(m)/g_de_ite/10;
             end                                                     %得到初始步长
             phi_ite0=phi_ite(m)+t*g_de_ite;                         %将初始步长代入后的当前源的试探性自变量
             Sum_G=0;                                                %计算这个试探性自变量的比较函数值
             for p3=1:N-1
                 s=0;
                 for q3=1:N-p3
                     s=s+R_ite(q3+p3,q3);
                 end
                 Sum_G=Sum_G+exp(1i*pi*cos(phi_ite0)*p3)*s;          %这里有自变量
             end
             g_ite0=Sum_g+2*real(Sum_G);                             %得到第一个源的试探性比较函数值

             while g_ite0<g_ite+0.3*t*(g_de_ite)^(2)                 %线性回溯法
                 t=t*0.5;
                 phi_ite0=phi_ite(m)+t*g_de_ite;                     %更新试探性自变量
                 Sum_G=0;                                            %更新当前试探性自变量的比较函数值
                 for p4=1:N-1
                     s=0;
                     for q4=1:N-p4
                         s=s+R_ite(q4+p4,q4);
                     end
                     Sum_G=Sum_G+exp(1i*pi*cos(phi_ite0)*p4)*s;      %这里有自变量
                 end
                 g_ite0=Sum_g+2*real(Sum_G);                         %得到新的试探性比较函数值
             end                                                     %得到想要的步长

             phi_ite(m)=phi_ite0;                                    %得到当前源下一次迭代的自变量
             g_ite=g_ite0;                                           %得到下一次迭代的目标函数值

             sum_de=0;                                               %更新比较函数的导数值
             for n5=1:N-1
                 s=0;
                 for q5=1:N-n5
                     s=s+R_ite(q5+n5,q5);
                 end
                 sum_de=sum_de+n5*exp(1i*pi*cos(phi_ite(m))*n5)*s;
             end
             g_de_ite=2*pi*sin(phi_ite(m))*imag(sum_de);             %得到当前源比较函数的导数值
       end                                                           %得到当前源的方位角 梯度下降法结束
       P_ite(m)=max((g_ite/sum(1./No_ite)-1)/sum(1./No_ite),0);      %得到当前源的功率

       for n_4=1:N                                                   %更新当前源的阵列导向矢量
           A_ite(n_4,m)=exp(-1i*pi*cos(phi_ite(m))*(n_4-1));
       end
       C_y=diag(No_ite);                                             %定义用于迭代的样本协方差矩阵
       for m_4=1:M
           C_y=C_y+P_ite(m_4)*A_ite(:,m_4)*(A_ite(:,m_4)');
       end                                                           %初始化样本协方差矩阵
       invC_y=inv(C_y);                                              %使用逆矩阵
   end                                                               %得到本次迭代中所有源的方位角和功率和方向矩阵
   
   while norm((phi_ite-phi_ite1)/pi*180)>0.001                         %算法停止条件
         Bb_ite=invC_y*A_ite*diag(P_ite);                                %更新源和噪声的协方差矩阵
         Cc=(Bb_ite')*(R_y-C_y)*Bb_ite;
         for mm=1:M
             P_ite(mm)=real(Cc(mm,mm))+P_ite(mm);
         end
         D_z=diag(No_ite);                                                 %更新噪声参数估计值
         R_z=D_z*invC_y*R_y*invC_y*D_z+D_z-D_z*invC_y*D_z;
         for nn=1:N
             No_ite(nn)=real(R_z(nn,nn));
             %No_ite(nn)=0.9*No_ite(nn)+0.1*real(R_z(nn,nn));
         end
         C_y=diag(No_ite);                                                 %定义用于迭代的样本协方差矩阵
         for mm3=1:M
             C_y=C_y+P_ite(mm3)*A_ite(:,mm3)*(A_ite(:,mm3)');
         end                                                               %初始化样本协方差矩阵
         invC_y=inv(C_y);                                                  %协方差逆矩阵
                                                                           %本次迭代结束
         k=k+1;                                                            %下一次迭代开始
         phi_ite1=phi_ite;                                                 %存储上一次值
         for m=1:M                                                         %第m个源
             for im=1:M                                                    %更新其他源的功率
                 if im~=m
                    b_i=invC_y*A_ite(:,im)*P_ite(im);
                    gr=real((b_i')*R_y*b_i);
                    P_ite(im)=gr+P_ite(im)*max(1-real((A_ite(:,im)')*b_i),0);
                 end
             end
             C_m=P_ite(m)*A_ite(:,m)*(A_ite(:,m)')+diag(No_ite);           %当前源的协方差矩阵
             R_R=C_m*invC_y*R_y*invC_y*C_m+C_m-C_m*invC_y*C_m;             %得到当前源的完整数据协方差矩阵 复数矩阵
             inv_N=diag(1./No_ite);                                        %对噪声协方差矩阵求逆
             R_ite=inv_N*R_R*inv_N;                                        %得到当前源的变形完整数据协方差矩阵 复数矩阵

             Sum_g=0;                                                      %计算当前源的比较函数值  梯度下降法开始
             for g1=1:N                                                    %第一部分与角无关，可以重复使用，均匀线阵的函数
                 Sum_g=Sum_g+real(R_ite(g1,g1));                           %取实部保证实数
             end
             Sum_G=0;                                                      %第二部分与角有关，每次更新方位角后都要重新计算
             for p_1=1:N-1
                 s=0;
                 for q_1=1:N-p_1
                     s=s+R_ite(q_1+p_1,q_1);
                 end
                 Sum_G=Sum_G+exp(1i*pi*cos(phi_ite(m))*p_1)*s;             %这里有当前源的方位角
             end
             g_ite=Sum_g+2*real(Sum_G);                                    %得到当前源的比较函数值 实数

             sum_de=0;                                                     %计算当前源比较函数的导数值
             for n3=1:N-1
                 s=0;
                 for qq=1:N-n3
                     s=s+R_ite(qq+n3,qq);
                 end
                 sum_de=sum_de+n3*exp(1i*pi*cos(phi_ite(m))*n3)*s;
             end
             g_de_ite=2*pi*sin(phi_ite(m))*imag(sum_de);                   %得到当前源比较函数的导数值

             while abs(g_de_ite)>0.001                                     %判断导数值大小
                   if g_de_ite>0                                           %初始化步长
                      t=(pi-phi_ite(m))/g_de_ite/10;
                   else
                      t=-phi_ite(m)/g_de_ite/10;
                   end                                                     %得到初始步长
                   phi_ite0=phi_ite(m)+t*g_de_ite;                         %将初始步长代入后的当前源的试探性自变量
                   Sum_G=0;                                                %计算这个试探性自变量的比较函数值
                   for p3=1:N-1
                       s=0;
                       for q3=1:N-p3
                           s=s+R_ite(q3+p3,q3);
                       end
                       Sum_G=Sum_G+exp(1i*pi*cos(phi_ite0)*p3)*s;          %这里有自变量
                   end
                   g_ite0=Sum_g+2*real(Sum_G);                             %得到第一个源的试探性比较函数值

                   while g_ite0<g_ite+0.3*t*(g_de_ite)^(2)                 %线性回溯法
                       t=t*0.5;
                       phi_ite0=phi_ite(m)+t*g_de_ite;                     %更新试探性自变量
                       Sum_G=0;                                            %更新当前试探性自变量的比较函数值
                       for p4=1:N-1
                           s=0;
                           for q4=1:N-p4
                               s=s+R_ite(q4+p4,q4);
                           end
                           Sum_G=Sum_G+exp(1i*pi*cos(phi_ite0)*p4)*s;      %这里有自变量
                       end
                       g_ite0=Sum_g+2*real(Sum_G);                         %得到新的试探性比较函数值
                   end                                                     %得到想要的步长

                   phi_ite(m)=phi_ite0;                                    %得到当前源下一次迭代的自变量
                   g_ite=g_ite0;                                           %得到下一次迭代的目标函数值

                   sum_de=0;                                               %更新比较函数的导数值
                   for n5=1:N-1
                       s=0;
                       for q5=1:N-n5
                           s=s+R_ite(q5+n5,q5);
                       end
                       sum_de=sum_de+n5*exp(1i*pi*cos(phi_ite(m))*n5)*s;
                   end
                   g_de_ite=2*pi*sin(phi_ite(m))*imag(sum_de);             %得到当前源比较函数的导数值
             end                                                           %得到当前源的方位角 梯度下降法结束
             P_ite(m)=max((g_ite/sum(1./No_ite)-1)/sum(1./No_ite),0);      %得到当前源的功率

             for n_4=1:N                                                   %更新当前源的阵列导向矢量
                 A_ite(n_4,m)=exp(-1i*pi*cos(phi_ite(m))*(n_4-1));
             end
             C_y=diag(No_ite);                                             %定义用于迭代的样本协方差矩阵
             for m_4=1:M
                 C_y=C_y+P_ite(m_4)*A_ite(:,m_4)*(A_ite(:,m_4)');
             end                                                           %初始化样本协方差矩阵
             invC_y=inv(C_y);                                              %使用逆矩阵
         end                                                               %得到本次迭代中所有源的方位角和功率和方向矩阵
   end                                                                     %实验结束
   X1(ii)=X1(ii)+((phi_ite(1)-phi_tr(1))/pi*180)^(2);                      %第num次实验的第一个源的方位角 用于画图
   Y1(ii)=Y1(ii)+((phi_ite(2)-phi_tr(2))/pi*180)^(2);                      %第num次实验的第二个源的方位角 用于画图
   
   %ECME算法
   phi_ite=[45,95]/180*pi;                                                 %两个源的方位角
   A_ite=zeros(N,M);                                                       %定义用于迭代的阵列流型矩阵
   for m2=1:M                                                              %初始化阵列流型矩阵
       for n2=1:N
           A_ite(n2,m2)=exp(-1i*pi*cos(phi_ite(m2))*(n2-1));
       end
   end
   P_ite=eye(M);Q_ite=eye(N);                                              %初始化信号和噪声的功率协方差矩阵并用于迭代
   C_ite=A_ite*P_ite*(A_ite')+Q_ite;invC_ite=inv(C_ite);
   k=0;                                                                    %初始化迭代数
   k=k+1;                                                              %第k次迭代开始
   phi_ite3=phi_ite;                                                   %存储上一次
                                                                       %EM算法
   Bb_ite=invC_ite*A_ite*P_ite;                                        %更新信号和噪声的协方差矩阵
   P_ite=(Bb_ite')*(R_y-C_ite)*Bb_ite+P_ite;
   QC_ite=Q_ite*invC_ite;
   %Q_ite=QC_ite*(R_y-C_ite)*(QC_ite')+Q_ite;                           %有色噪声
   Gh=QC_ite*(R_y-C_ite)*(QC_ite')+Q_ite;                              %非均匀噪声
   for n=1:N
       Q_ite(n,n)=real(Gh(n,n));
   end                                                                 %EM算法结束 开始更新到达方向
   C_ite=A_ite*P_ite*(A_ite')+Q_ite;invC_ite=inv(C_ite);               %更新协方差矩阵
   f_ite=real(log(det(C_ite))+trace(invC_ite*R_y));                    %更新目标值
                                                                           %现在使用梯度下降法更新角度
       %A_com=zeros(N,M);                                                   %定义用于比较的阵列流型矩阵
       %for i1=1:Nu-1                                                       %比较网格的角度 二维遍历
       %    for j1=1:Nu-1
       %        phi_com=da*[i1,j1];                                         %处于某个网格上两个源的方位角
       %        for m3=1:M                                                  %更新比较阵列流型矩阵
       %            for n3=1:N
       %                A_com(n3,m3)=exp(-1i*pi*cos(phi_com(m3))*(n3-1));
       %            end
       %        end
       %        C_com=A_com*P_ite*(A_com')+Q_ite;invC_com=inv(C_com);       %更新比较阵列协方差矩阵
       %        f_com=real(log(det(C_com))+trace(invC_com*R_y));            %更新比较目标值
       %        if f_com<f_ite                                              %
       %            phi_ite=phi_com;A_ite=A_com;
       %            C_ite=C_com;invC_ite=invC_com;
       %            f_ite=f_com;
       %        end
       %    end
       %end                                                                 %得到接近最优的角度

   D=zeros(N,M);                                                       %定义阵列导数矩阵 复数矩阵
   for n4=1:N
       for m4=1:M
       D(n4,m4)=1i*(n4-1)*pi*sin(phi_ite(m4))*exp(-1i*(n4-1)*pi*cos(phi_ite(m4)));
       end
   end                                                                 %得到阵列导数阵列
   de_ite=zeros(1,M);                                                  %定义角度梯度
   for m=1:M
   de_ite(m)=2*real(P_ite(m,:)*(A_ite')*invC_ite*(C_ite-R_y)*invC_ite*D(:,m));
   end

   while norm(de_ite)>0.001                                            %判断导数范数的大小
         %t=da/max(abs(de_ite));                                        %初始步长
         t_i=zeros(1,M);
         for ti=1:M
             if de_ite(ti)<0                                           %初始化步长
                t_i(ti)=-(pi-phi_ite(m))/de_ite(ti);
             else
                t_i(ti)=phi_ite(m)/de_ite(ti);
             end
         end
         t=min(t_i)/10;                                                %得到初始步长
         phi_ite0=phi_ite-t*de_ite;                                    %试探性自变量 负梯度方向
         A_ite0=zeros(N,M);                                            %定义试探性阵列矩阵
         for m5=1:M                                                    %更新试探性阵列流型矩阵
             for n5=1:N
                 A_ite0(n5,m5)=exp(-1i*pi*cos(phi_ite0(m5))*(n5-1));
             end
         end
         C_ite0=A_ite0*P_ite*(A_ite0')+Q_ite;invC_ite0=inv(C_ite0);
         f_ite0=real(log(det(C_ite0))+trace(invC_ite0*R_y));

         while f_ite0>f_ite-0.3*t*(norm(de_ite)^(2))                   %线性回溯法
               t=t*0.5;
               phi_ite0=phi_ite-t*de_ite;                              %更新试探性自变量
               for m6=1:M                                              %更新阵列流型矩阵
                   for n6=1:N
                       A_ite0(n6,m6)=exp(-1i*pi*cos(phi_ite0(m6))*(n6-1));
                   end
               end
               C_ite0=A_ite0*P_ite*(A_ite0')+Q_ite;invC_ite0=inv(C_ite0);
               f_ite0=real(log(det(C_ite0))+trace(invC_ite0*R_y));
         end                                                           %得到满足要求的步长

         phi_ite=phi_ite0;A_ite=A_ite0;                                %得到当前源下一次的迭代自变量
         C_ite=C_ite0;invC_ite=invC_ite0;
         f_ite=f_ite0;                                                 %得到下一次迭代的目标函数值
                          
         for n7=1:N                                                    %更新阵列导数矩阵
             for m7=1:M
             D(n7,m7)=1i*(n7-1)*pi*sin(phi_ite(m7))*exp(-1i*(n7-1)*pi*cos(phi_ite(m7)));
             end
         end
         for m=1:M                                                     %更新梯度
         de_ite(m)=2*real(P_ite(m,:)*(A_ite')*invC_ite*(C_ite-R_y)*invC_ite*D(:,m));
         end
   end                                                                 %得到方位角 梯度下降法结束
                                                                       %本次迭代结束
   while norm((phi_ite3-phi_ite)/pi*180)>0.001                             %算法停止条件
       k=k+1;                                                              %第k次迭代开始
       phi_ite3=phi_ite;                                                   %存储上一次值
                                                                           %EM算法
       Bb_ite=invC_ite*A_ite*P_ite;                                        %更新信号和噪声的协方差矩阵
       P_ite=(Bb_ite')*(R_y-C_ite)*Bb_ite+P_ite;
       QC_ite=Q_ite*invC_ite;
       %Q_ite=QC_ite*(R_y-C_ite)*(QC_ite')+Q_ite;                           %有色噪声
       Gh=QC_ite*(R_y-C_ite)*(QC_ite')+Q_ite;                              %非均匀噪声
       for n=1:N
           Q_ite(n,n)=real(Gh(n,n));
       end                                                                 %EM算法结束
       C_ite=A_ite*P_ite*(A_ite')+Q_ite;invC_ite=inv(C_ite);               %更新协方差矩阵
       f_ite=real(log(det(C_ite))+trace(invC_ite*R_y));                    %更新目标值
                                                                           %现在使用梯度下降法更新角度
       %A_com=zeros(N,M);                                                   %定义用于比较的阵列流型矩阵
       %for i1=1:Nu-1                                                       %比较网格的角度 二维遍历
       %    for j1=1:Nu-1
       %        phi_com=da*[i1,j1];                                         %处于某个网格上两个源的方位角
       %        for m3=1:M                                                  %更新比较阵列流型矩阵
       %            for n3=1:N
       %                A_com(n3,m3)=exp(-1i*pi*cos(phi_com(m3))*(n3-1));
       %            end
       %        end
       %        C_com=A_com*P_ite*(A_com')+Q_ite;invC_com=inv(C_com);       %更新比较阵列协方差矩阵
       %        f_com=real(log(det(C_com))+trace(invC_com*R_y));            %更新比较目标值
       %        if f_com<f_ite                                              %
       %            phi_ite=phi_com;A_ite=A_com;
       %            C_ite=C_com;invC_ite=invC_com;
       %            f_ite=f_com;
       %        end
       %    end
       %end                                                                 %得到接近最优的角度

       D=zeros(N,M);                                                       %定义阵列导数矩阵 复数矩阵
       for n4=1:N
           for m4=1:M
           D(n4,m4)=1i*(n4-1)*pi*sin(phi_ite(m4))*exp(-1i*(n4-1)*pi*cos(phi_ite(m4)));
           end
       end                                                                 %得到阵列导数阵列
       de_ite=zeros(1,M);                                                  %定义角度梯度
       for m=1:M
       de_ite(m)=2*real(P_ite(m,:)*(A_ite')*invC_ite*(C_ite-R_y)*invC_ite*D(:,m));
       end

       while norm(de_ite)>0.001                                            %判断导数范数的大小
             %t=da/max(abs(de_ite));                                        %初始步长
             t_i=zeros(1,M);
             for ti=1:M
                 if de_ite(ti)<0                                           %初始化步长
                    t_i(ti)=-(pi-phi_ite(m))/de_ite(ti);
                 else
                    t_i(ti)=phi_ite(m)/de_ite(ti);
                 end
             end
             t=min(t_i)/10;                                                %得到初始步长
             phi_ite0=phi_ite-t*de_ite;                                    %试探性自变量 负梯度方向
             A_ite0=zeros(N,M);                                            %定义试探性阵列矩阵
             for m5=1:M                                                    %更新试探性阵列流型矩阵
                 for n5=1:N
                     A_ite0(n5,m5)=exp(-1i*pi*cos(phi_ite0(m5))*(n5-1));
                 end
             end
             C_ite0=A_ite0*P_ite*(A_ite0')+Q_ite;invC_ite0=inv(C_ite0);
             f_ite0=real(log(det(C_ite0))+trace(invC_ite0*R_y));

             while f_ite0>f_ite-0.3*t*(norm(de_ite)^(2))                   %线性回溯法
                   t=t*0.5;
                   phi_ite0=phi_ite-t*de_ite;                              %更新试探性自变量
                   for m6=1:M                                              %更新阵列流型矩阵
                       for n6=1:N
                           A_ite0(n6,m6)=exp(-1i*pi*cos(phi_ite0(m6))*(n6-1));
                       end
                   end
                   C_ite0=A_ite0*P_ite*(A_ite0')+Q_ite;invC_ite0=inv(C_ite0);
                   f_ite0=real(log(det(C_ite0))+trace(invC_ite0*R_y));
             end                                                           %得到满足要求的步长

             phi_ite=phi_ite0;A_ite=A_ite0;                                %得到当前源下一次的迭代自变量
             C_ite=C_ite0;invC_ite=invC_ite0;
             f_ite=f_ite0;                                                 %得到下一次迭代的目标函数值
                          
             for n7=1:N                                                    %更新阵列导数矩阵
                 for m7=1:M
                 D(n7,m7)=1i*(n7-1)*pi*sin(phi_ite(m7))*exp(-1i*(n7-1)*pi*cos(phi_ite(m7)));
                 end
             end
             for m=1:M                                                     %更新梯度
             de_ite(m)=2*real(P_ite(m,:)*(A_ite')*invC_ite*(C_ite-R_y)*invC_ite*D(:,m));
             end
       end                                                                 %得到方位角 梯度下降法结束
   end                                                                     %实验结束
   X2(ii)=X2(ii)+((phi_ite(1)-phi_tr(1))/pi*180)^(2);                      %第num次实验的第一个源的方位角 用于画图
   Y2(ii)=Y2(ii)+((phi_ite(2)-phi_tr(2))/pi*180)^(2);                      %第num次实验的第二个源的方位角 用于画图
   
   %新的较差子空间法
   D_3=diag(R_y);                                                          %取对角元素 矢量
   [F_3,sD]=sort(D_3);                                                     %将对角元素按从小到大排序
   R3=R_y-diag(D_3);                                                       %构造对角元素为0的类协方差矩阵
   [V3,D3]=eig(R3);                                                        %特征值分解 V3存储特征向量 D3是特征值对角矩阵
   D3=abs(diag(D3));                                                       %特征值转为实数 存储在一个向量中
   [F3,index3]=sort(D3);                                                   %对特征值从小到大排序
   U_n3=zeros(N,N-M);                                                      %构造噪声特征向量
   for ii3=1:N-M                                                           %排序后前N-M个特征值对应的向量
       U_n3(:,ii3)=V3(:,index3(ii3));
   end
   
   k3=sD(1);ev=zeros(N,1);ev(k3,1)=1;                                      %单位矢量
   UU3=U_n3*(U_n3');Q_0=abs((ev')*R_y*UU3*ev/(UU3(k3,k3)));
   Q_g=diag(D_3)-D3(1)*eye(N);
   Q=Q_g+Q_0*eye(N);
   [V4,D4]=eig(R_y,Q);                                                     %广义特征值分解 V4存储广义特征向量 D4是广义特征值对角矩阵
   D4=abs(diag(D4));                                                       %特征值转为实数 存储在一个向量中
   [F4,index4]=sort(D4);                                                   %对特征值从小到大排序
   U_n4=zeros(N,N-M);                                                      %构造噪声特征向量
   for ii4=1:N-M                                                           %排序后前N-M个特征值对应的向量
       U_n4(:,ii4)=V4(:,index4(ii4));
   end
   
   c=U_n4(1,:);                                                             %第一行向量
   U_nn=U_n4(2:N,:);                                                        %剩下的行
   v=zeros(N,1);v(1,1)=1;v(2:N,1)=U_nn*(c')/(norm(c)^2);                      %得到最优向量
   r=roots(fliplr(v'));                                                    %求根
   [FF,ind]=sort(abs(abs(r)-ones(N-1,1)));                                 %靠近单位圆排序
   angel1=acos(angle(r(ind(1)))/(-pi));
   angel2=acos(angle(r(ind(2)))/(-pi));
   ass1=sort([angel1,angel2]);
   X3(ii)=X3(ii)+((ass1(1)-phi_tr(1))/pi*180)^(2);                      %第num次实验的第一个源的方位角 用于画图
   Y3(ii)=Y3(ii)+((ass1(2)-phi_tr(2))/pi*180)^(2);                      %第num次实验的第二个源的方位角 用于画图
   
   %最新的较差子空间法
   D_ite=diag(diag(R_y));                                                  %用快拍协方差矩阵的对角元素初始化噪声协方差矩阵
   for iq=1:20
       D_1=sqrt(D_ite);D_2=inv(D_1);
       [V5,D5]=eig(R_y,D_ite);                                             %广义特征值分解 V5存储广义特征向量 D5是广义特征值对角矩阵
       D5=abs(diag(D5));                                                   %特征值转为实数 存储在一个向量中
       [F5,index5]=sort(D5);                                               %对特征值从小到大排序
       U_n5=zeros(N,N-M);                                                  %构造噪声特征向量
       for ii5=1:N-M                                                       %排序后前N-M个特征值对应的向量
           U_n5(:,ii5)=V5(:,index5(ii5));
       end
       UU5=U_n5*(U_n5');                                                   %得到噪声子空间
       D_ite=1/2*abs(diag(diag(R_y*UU5+UU5*R_y))*inv(diag(diag(UU5))));    %更新噪声协方差矩阵
   end
   
   c=U_n5(1,:);                                                             %第一行向量
   U_nn=U_n5(2:N,:);                                                        %剩下的行
   v=zeros(N,1);v(1,1)=1;v(2:N,1)=U_nn*(c')/(norm(c)^2);                      %得到最优向量
   r=roots(fliplr(v'));                                                    %求根
   [FF,ind]=sort(abs(abs(r)-ones(N-1,1)));                                 %靠近单位圆排序
   angel1=acos(angle(r(ind(1)))/(-pi));
   angel2=acos(angle(r(ind(2)))/(-pi));
   ass1=sort([angel1,angel2]);
   X4(ii)=X4(ii)+((ass1(1)-phi_tr(1))/pi*180)^(2);                      %第num次实验的第一个源的方位角 用于画图
   Y4(ii)=Y4(ii)+((ass1(2)-phi_tr(2))/pi*180)^(2);                      %第num次实验的第二个源的方位角 用于画图
   end                                                                     %全部实验结束
   end
   
   U=10*(1:10);
   subplot(2,1,1);semilogy(U,sqrt(X1/number));hold on;subplot(2,1,2);semilogy(U,sqrt(Y1/number));
   hold on;subplot(2,1,1);semilogy(U,sqrt(X2/number));hold on;subplot(2,1,2);semilogy(U,sqrt(Y2/number));
   hold on;subplot(2,1,1);semilogy(U,sqrt(X3/number));hold on;subplot(2,1,2);semilogy(U,sqrt(Y3/number));
   hold on;subplot(2,1,1);semilogy(U,sqrt(X4/number));hold on;subplot(2,1,2);semilogy(U,sqrt(Y4/number));
   hold on;subplot(2,1,1);semilogy(U,CRLB1);hold on;subplot(2,1,2);semilogy(U,CRLB2);