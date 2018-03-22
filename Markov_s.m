% Markov motion model
mu = [0 50];
sigma = [10 0; 0 10];
L = 100;
delt = 1;
[X, Y] = meshgrid(0:0.5:L, 0:0.5:L);
P = mvnpdf([X(:) Y(:)],mu,sigma);
P = reshape(P,size(X));

theta_min = -pi/4;
theta_max = -pi/8;
n_theta = 20;
theta_1 = linspace(theta_min, theta_max, n_theta);
f_theta = 1/n_theta;
mu_v = 2;
sigma_v = 5;
n_v=30;
v_1 = linspace((mu_v-3*sigma_v),(mu_v+3*sigma_v),n_v);
delv = v_1(2) - v_1(1);
 
for i=1:delt:10
    P_new = zeros(size(P));
    for index_x = 1:length(X(1,:))
        for index_y = 1:length(Y(:,1))
            temp = 0;
            for index_theta=1:length(theta_1)
                for index_v=1:length(v_1)
                    p_theta = f_theta;
                    p_v = normpdf(v_1(index_v), mu_v, sigma_v)*delv;
                    
                    x_prev = X(1,index_x) - v_1(index_v)*cos(theta_1(index_theta))*delt;
                    y_prev = Y(index_y,1) - v_1(index_v)*sin(theta_1(index_theta))*delt;
                    
                    index_x_prev = round(x_prev/0.5)+1;
                    index_y_prev = round(y_prev/0.5)+1;
                    
                    if (index_x_prev*index_y_prev < 1 || index_x_prev > length(X(1,:)) || index_y_prev > length(Y(:,1)))
                        index_x_prev = 1;
                        index_y_prev = 1;
                        p_v = 0;
                        p_theta = 0;
                    end
                    
                    temp = temp + P(index_y_prev, index_x_prev)*p_v*p_theta;
                end
            end
            P_new(index_y,index_x) = temp;
        end
    end
    P = P_new;
    figure();
    mesh(X, Y, P);
    xlabel('x_T [m]');
    ylabel('y_T [m]');
    zlabel('Probability');
    view(50,50);
    title(['time = ' num2str(i) '(sec)']);
    colorbar;
end
