function main
    plot_flag = true; %choose whether to show plots or not

    [t, y_L, v_L, y_R, v_R] = load_motor_data(plot_flag);
    range = [100:200];
    single_v_L = v_L(range);
    single_t = t(range) - t(range(1));
    [alpha, c] = fitData(single_t, single_v_L, plot_flag)
    [t, theta] = load_pendulum_data(plot_flag);

    % U_step * beta = c
    U_step = 300 %just 1 bc -300,300 is -1,1 in motor space?
    beta = c/U_step
    tau = 1/alpha
    

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
        figure( 'Name', 'untitled fit 1' );
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
        figure();
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
        figure();
        hold on
        plot(t,theta,'c','linewidth',1);
        xlabel('time (sec)'); ylabel('angle (rad)');
        title('Pendulum Calibration Data');
    end
end
