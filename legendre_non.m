function p=legendre_non(x,n)
%这里的x是地心余纬
%输入参数:
%输出参数：P，
p=zeros(n+1,n+1);
t=cos(x);
u=sin(x);
p(1,1)=1;
p(2,1)=sqrt(3)*t;
p(2,2)=sqrt(3)*u;

for i=3:n+1
    for j=1:2
        p(i,j)=sqrt((2*i-3)*(2*i-1)/((i-j)*(i+j-2)))*t*p(i-1,j)-sqrt((2*i-1)*(i+j-3)*(i-j-1)/((2*i-5)*(i+j-2)*(i-j)))*p(i-2,j);
    end
end

for i=3:n+1
    for j=3:i
        if((i-j)>=1)
            if(j==3)
                p(i,j)=sqrt((2*i-1)*(i-j)*(i-j-1)/((2*i-5)*(i+j-2)*(i+j-3)))*p(i-2,j)+sqrt(2)*sqrt((2*i-1)*(i+j-4)*(i+j-5)/((2*i-5)*(i+j-2)*(i+j-3)))*p(i-2,j-2)-sqrt(2)*sqrt((i-j+1)*(i-j+2)/((i+j-2)*(i+j-3)))*p(i,j-2);
            else
                p(i,j)=sqrt((2*i-1)*(i-j)*(i-j-1)/((2*i-5)*(i+j-2)*(i+j-3)))*p(i-2,j)+sqrt((2*i-1)*(i+j-4)*(i+j-5)/((2*i-5)*(i+j-2)*(i+j-3)))*p(i-2,j-2)-sqrt((i-j+1)*(i-j+2)/((i+j-2)*(i+j-3)))*p(i,j-2);
            end
        else
            if(j==3)
                p(i,j)=sqrt(2)*sqrt((2*i-1)*(i+j-4)*(i+j-5)/((2*i-5)*(i+j-2)*(i+j-3)))*p(i-2,j-2)-sqrt(2)*sqrt((i-j+1)*(i-j+2)/((i+j-2)*(i+j-3)))*p(i,j-2);
            else
                p(i,j)=sqrt((2*i-1)*(i+j-4)*(i+j-5)/((2*i-5)*(i+j-2)*(i+j-3)))*p(i-2,j-2)-sqrt((i-j+1)*(i-j+2)/((i+j-2)*(i+j-3)))*p(i,j-2);
            end
        end
    end
end
end
