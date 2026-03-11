function main
    plot_flag = true; %choose whether to show plots or not

    [t, y_L, v_L, y_R, v_R] = load_motor_data(plot_flag);
    range = [100:200];
    single_v_L = v_L(range);
    single_t = t(range) - t(range(1));
    [alpha, c] = fitData(single_t, single_v_L, plot_flag)
    [t_gyro, theta] = load_pendulum_data(plot_flag);
    

    % U_step * beta = c
    U_step = 300; %just 1 bc -300,300 is -1,1 in motor space?
    beta = c/U_step;
    tau = 1/alpha;

    start_post = 1;
    end_post = 1;
    while t_gyro(start_post) < 5
        start_post = start_post + 1;
    end
    while t_gyro(end_post) < 15
        end_post = end_post + 1;
    end

    t_fit_gyro = t_gyro(start_post:end_post) - t_gyro(start_post);
    theta_fit_gyro = theta(start_post:end_post);

    [wd, sigma, wn, zeta] = fit_underdamped(t_fit_gyro, theta_fit_gyro, plot_flag)
    g = 9.81; % gravitational acceleration in m/s^2
    leff = g/(wn^2)

end


function [a, c] = fitData(t, v_L, plot_flag)
    [xData, yData] = prepareCurveData( t, v_L );
    
    % Set up fittype and options.
    ft = fittype( 'c*(1-exp(-a*x))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [1 0.0461713906311539];
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    
    % Plot fit with data.
    if plot_flag
        figure(1);
        h = plot( fitresult, xData, yData );
        legend( h, 'v_L vs. t', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
        % Label axes
        xlabel( 't', 'Interpreter', 'none' );
        ylabel( 'v_L', 'Interpreter', 'none' );
        grid on
    end
    a = fitresult.a;
    c = fitresult.c;
end

function [t, y_L, v_L, y_R, v_R] = load_motor_data(plot_flag)
    %path and file name of data
    fpath = './'; %path (change this!)
    fname_in = 'motor_calibration.txt'; %file name (change this!)
    %load the motor calibration data
    motor_data = importdata([fpath,fname_in]);
    %unpack the motor calibration data
    t = motor_data(:,1);
    y_L = motor_data(:,2); v_L = motor_data(:,3);
    y_R = motor_data(:,4); v_R = motor_data(:,5);
    %plot the motor calibration data
    if plot_flag
        figure(2);
        subplot(2,1,1);
        hold on
        plot(t,v_L,'c','linewidth',1);
        plot(t,v_R,'r','linewidth',1);
        xlabel('time (sec)'); ylabel('wheel speed (m/sec)');
        title('Motor Calibration Data');
        h1 = legend('Left Wheel','Right Wheel');
        set(h1,'location','southeast');
        subplot(2,1,2);
        hold on
        plot(t,y_L,'c','linewidth',1);
        plot(t,y_R,'r--','linewidth',1);
        xlabel('time (sec)'); ylabel('wheel command (-)');
        title('Motor Calibration Data');
        h2 = legend('Left Wheel','Right Wheel');
        set(h2,'location','southeast');
    end
end

%loads and plots the pendulum calibration data
function [t, theta] = load_pendulum_data(plot_flag)
    %path and file name of data
    fpath = './'; %path (change this!)
    fname_in = 'pendulum_calibration_data.txt'; %file name (change this!)
    %load the pendulum calibration data
    pendulum_data = importdata([fpath,fname_in]);
    %unpack the pendulum calibration data
    t = pendulum_data(:,1) - pendulum_data(1, 1); theta = pendulum_data(:,2);
    %plot the motor calibration data
    if plot_flag
        figure(3);
        hold on
        plot(t,theta,'c','linewidth',1);
        xlabel('time (sec)'); ylabel('angle (rad)');
        title('Pendulum Calibration Data');
    end
end

function [omega_d, sigma, omega_n, zeta] = fit_underdamped(t, theta, plot_flag)
    %runs the custom curve fitting tool
    signal_params = regress_underdamped_response(t,theta);
    %unpack results
    C = signal_params.C; S = signal_params.S; D = signal_params.D;
    omega_d = signal_params.omega_d; sigma = signal_params.sigma;
    omega_n = signal_params.omega_n; zeta = signal_params.zeta;

    if plot_flag
        % create a function for time values
        y_fit = D + (C*cos(omega_d * t) + S*sin(omega_d * t)) .* exp(-sigma * t);

        %plot step response data
        figure(4)
        hold on
        plot(t,theta,'r.','DisplayName','theta');
        plot(t, y_fit, 'b', 'DisplayName','theta fit')
        %create axis labels, title and legend
        title('Measured v. Fit Angle of Pendulum System');
        my_legend = legend();
        set(my_legend,'location','southeast');
        xlabel('t (sec)');
        ylabel('theta (radians)');
    end
end