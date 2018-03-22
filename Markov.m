clc
clear
close all

fieldx=0:1:100;
fieldy=0:1:100;
diag=[10 0;0 10];
mu=[0;50];

n=20;

v1=-15:35/n:20;
for i=1:length(v1)
    fv1(i)=1/(sqrt(2*pi)*5)*exp(-1/2*((v1(i)-2)/5)^2);
end

figure
plot(v1,fv1)

v2=-5:11/n:6;
for i=1:length(v2)
    fv2(i)=1/(sqrt(2*pi)*1)*exp(-1/2*((v2(i)-1)/1)^2);
end
figure
plot(v2,fv2)

theta2=-5:10/n:5;
for i=1:length(theta2)
    ftheta2(i)=1/(sqrt(2*pi)*(pi/8))*exp(-1/2*((theta2(i)+pi/6)/(pi/8))^2);
end
figure
plot(theta2,ftheta2)

theta1=-pi/4:(-pi/8+pi/4)/n:-pi/8;
ftheta1=ones(1,n+1)*1/(n+1);
figure
plot(-pi/4:(pi/8+pi/4)/n:pi/8,ftheta1)


for i=1:length(fieldx)
    for j=1:length(fieldy)
    	f(j,i)=1/(2*pi*det(diag^(1/2)))*exp(-1/2*([fieldx(i);fieldy(j)]-mu)'*inv(diag)*([fieldx(i);fieldy(j)]-mu));
    end
end




f(:,:,1)=f;

figure
mesh(fieldx,fieldy,f(:,:,1))

f(:,:,2:10)=zeros(101,101,9);



tempth=0;
m=2;

for m=2:10
    for i=1:length(fieldx)
        for j=1:length(fieldy)
            tempv=0;
            for k=1:n+1
                for kk=1:n+1        
                    temppos=[fieldx(i);fieldy(j)]-v1(kk).*[cos(theta1(k));sin(theta1(k))]*1;
                    idx=round(temppos/1)+1;
                    if idx(1)>101||idx(1)<1||idx(2)>101||idx(2)<1
                        temp1=0;
                    else
                        temp1=f(idx(2),idx(1),m-1)*ftheta1(k)*fv1(kk)*35/n;
                    end
                    tempv=tempv+temp1;
                end
            end
            f(j,i,m)=tempv;
        end
    end
    figure
    mesh(fieldx,fieldy,f(:,:,m))
    xlabel('x position (m)')
    ylabel('y position (m)')
    zlabel('Probability')
    xlim([0 100])
    ylim([0 100])
    zlim([0 6*10^-3])
%     view(50,50);
    colorbar;
end

for m=11:25
    for i=1:length(fieldx)
        for j=1:length(fieldy)
            tempv=0;
            for k=1:n+1
                for kk=1:n+1        
                    temppos=[fieldx(i);fieldy(j)]-v2(kk).*[cos(theta2(k));sin(theta2(k))]*1;
                    idx=round(temppos/1)+1;
                    if idx(1)>101||idx(1)<1||idx(2)>101||idx(2)<1
                        temp1=0;
                    else
                        temp1=f(idx(2),idx(1),m-1)*ftheta2(k)*fv2(kk)*11/n*10/n;
                    end
                    tempv=tempv+temp1;
                end
            end
            f(j,i,m)=tempv;
        end
    end
    figure
    mesh(fieldx,fieldy,f(:,:,m))
    xlabel('x position (m)')
    ylabel('y position (m)')
    zlabel('Probability')
    xlim([0 100])
    ylim([0 100])
    zlim([0 6*10^-3])
%     view(50,50);
    colorbar
end
%     figure
%     mesh(fieldx,fieldy,f1(:,:,m))
% end



