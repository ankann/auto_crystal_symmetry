clear variables

%symmetry operators
a_3=[-1/2 -sqrt(3)/2 0; sqrt(3)/2 -1/2 0; 0 0 1];
a_32=[-1/2 sqrt(3)/2 0; -sqrt(3)/2 -1/2 0; 0 0 1];
a_2=[-1 0 0; 0 -1 0;0 0 1];
a_s=[-1/2 -sqrt(3)/2 0; sqrt(3)/2 -1/2 0; 0 0 1]*[1 0 0; 0 1 0;0 0 -1];
a_h=[1 0 0; 0 1 0;0 0 -1];
a_v=[-1 0 0; 0 1 0;0 0 1];

%defining moduli as 'dd'
 dd=sym('dd',[3 6]);

for l=1:3
    for m=1:3
        for n=1:3
            d(l,m,n)=dd(l,ff(m,n));
        end
    end
end

sum_3=0;
sum_2=0;
sum_h=0;
sum_s=0;
sum_v=0;
sum_32=0;

for i=1:3
    for j=1:3
        for k=1:3
            sum_3=0;
            sum_2=0;
            sum_h=0;
            sum_s=0;
            sum_v=0;
            sum_32=0;
            for l=1:3
                for m=1:3
                    for n=1:3
                        sum_3=sum_3+a_3(i,l)*a_3(j,m)*a_3(k,n)*d(l,m,n);
                        sum_2=sum_2+a_2(i,l)*a_2(j,m)*a_2(k,n)*d(l,m,n);
                        sum_h=sum_h+a_h(i,l)*a_h(j,m)*a_h(k,n)*d(l,m,n);
                        sum_s=sum_s+a_s(i,l)*a_s(j,m)*a_s(k,n)*d(l,m,n);
                        sum_v=sum_v+a_v(i,l)*a_v(j,m)*a_v(k,n)*d(l,m,n);
                        sum_32=sum_32+a_32(i,l)*a_32(j,m)*a_32(k,n)*d(l,m,n);
                    end
                end
            end
            d_dash_3(i,j,k)=sum_3;
            d_dash_2(i,j,k)=sum_2;
            d_dash_h(i,j,k)=sum_h;
            d_dash_s(i,j,k)=sum_s;
            d_dash_32(i,j,k)=sum_32;
            d_dash_v(i,j,k)=sum_v;
        end
    end
end


%equating different combination of the symmetry operator
for i=1:3
    for j=1:3
        for k=1:3
%             eq(i,ff(j,k))= d_dash_2(i,j,k) == d_dash_3(i,j,k);
%             eq(3+i,3+ff(j,k))= d_dash_2(i,j,k) == d_dash_h(i,j,k);
             eq(6+i,6+ff(j,k))= d_dash_h(i,j,k) == d_dash_3(i,j,k);
              eq(9+i,9+ff(j,k))= d(i,j,k) == d_dash_3(i,j,k);
%             eq(12+i,12+ff(j,k))= d(i,j,k) == d_dash_2(i,j,k);
             eq(15+i,15+ff(j,k))= d(i,j,k) == d_dash_h(i,j,k);
            
%             eq(18+i,18+j)= d_dash_s(i,j,k) == d_dash_h(i,j,k);
%             eq(21+i,21+j)= d_dash_s(i,j,k) == d(i,j,k);
%             eq(24+i,24+j)= d_dash_s(i,j,k) == d_dash_3(i,j,k);
%             eq(27+i,27+j)= d_dash_s(i,j,k) == d_dash_2(i,j,k);
% 
%             eq(30+i,30+j)= d_dash_v(i,j,k) == d_dash_32(i,j,k);
%             eq(33+i,33+j)= d_dash_v(i,j,k) == d_dash_h(i,j,k);
%             eq(36+i,36+j)= d_dash_v(i,j,k) == d_dash_3(i,j,k);
%             eq(39+i,39+j)= d_dash_v(i,j,k) == d_dash_2(i,j,k);
%             eq(42+i,42+j)= d_dash_v(i,j,k) == d_dash_s(i,j,k);
%             eq(45+i,45+j)= d_dash_v(i,j,k) == d(i,j,k);
            
            %eq(39+i,39+j)= d_dash_32(i,j,k) == d_dash_2(i,j,k);
            %eq(42+i,42+ff(j,k))= d_dash_32(i,j,k) == d_dash_h(i,j,k);
            %eq(45+i,45+j)= d_dash_32(i,j,k) == d_dash_v(i,j,k);
            %eq(48+i,48+j)= d_dash_32(i,j,k) == d_dash_s(i,j,k);
            %eq(51+i,51+ff(j,k))= d_dash_32(i,j,k) == d_dash_3(i,j,k);
            %eq(54+i,54+ff(j,k))= d_dash_32(i,j,k) == d(i,j,k);
            
            %eq(54+i,54+j)= d_dash_v(i,j,k) == d_dash_2(i,j,k);
            eq(57+i,57+j)= d_dash_v(i,j,k) == d_dash_3(i,j,k);
            %eq(60+i,60+j)= d_dash_v(i,j,k) == d_dash_s(i,j,k);
            eq(63+i,63+j)= d_dash_v(i,j,k) == d_dash_h(i,j,k);
            eq(63+i,63+j)= d_dash_v(i,j,k) == d(i,j,k);
        end
    end
end

syms x
eq2= x == 0;
for i=1:length(eq)
    eq2=[eq2; eq(:,i)];
end

eq2=eq2(2:length(eq2));


s=0;
for i=1:length(eq2)
    if(eq2(i)~=0)
        s=s+1;
        eqq(s)=eq2(i);
    end
end

eqq=simplify(eqq);

s1=0;
for i=1:length(eqq)
    if(~eqq(i))
        s1=s1+1;
        eqqf(s1)=eqq(i);
    end
end

 
eqqf=simplify(eqqf);

eqqf_new=eqqf;

dd_new=dd;
clear dd

dd=sym('dd',[3 6]);

%equation list variable 'eqqf'
%equating different combination of the symmetry operator
for k=1:length(eqqf)
    s(k)=0;
    for i=1:3
        for j=1:6
            if(~eqqf(k))
                if(~isempty(solve(eqqf(k),dd_new(i,j))))
                    ssq=solve(eqqf(k),dd_new(i,j));
                    if(isstruct(ssq))
                    for l=1:3
                        for m=6
                               if(l~=i && m~=j)
                                    if(isfield(ssq,"dd"+int2str(l)+"_"+int2str(m)))
                                        aa=extractfield(ssq,"dd"+int2str(l)+"_"+int2str(m));
                                        if(~isempty(aa{1}))
                                            %if(aa{1}~=0)
                                                eqqf=subs(eqqf,dd_new(l,m),aa{1});
                                                eqqf=simplify(eqqf);
                                                dd=subs(dd,dd_new(l,m),aa{1});
                                                dd=simplify(dd);
                                        end
                                    end
                               end
                        end
                    end
                    else
                            eqqf=subs(eqqf,dd_new(i,j),ssq);
                            eqqf=simplify(eqqf);
                            %eqqf'
                            dd=subs(dd,dd_new(i,j),ssq);
                            dd=simplify(dd);
                            dd
                    end
                end
                eqqf=simplify(eqqf);
           end
       end
    end
end


%mapping in 3D
function g=ff(m,n)
    if (m==n)
        g=m;
    else
        if(m==1 && n==2 || m==2 && n==1)
            g=6;
        elseif(m==1 && n==3||m==3 && n==1)
            g=5;
        elseif(m==2 && n==3 ||m==3 && n==2)
            g=4;
        end
    end    
end