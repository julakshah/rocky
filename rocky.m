function main
    [t, y_L, v_L, y_R, v_R] = load_motor_data()
    range = [100:200]
    single_v_L = v_L(range)
    single_t = t(range) - t(range(1))
    fitData(single_t, single_v_L);
end


function fitData(t, v_L)
    [xData, yData] = prepareCurveData( t, v_L );
    
    % Set up fittype and options.
    ft = fittype( 'c*(1-exp(-a*x))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [1 0.0461713906311539];
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts )
    
    % Plot fit with data.
    figure( 'Name', 'untitled fit 1' );
    h = plot( fitresult, xData, yData );
    legend( h, 'v_L vs. t', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 't', 'Interpreter', 'none' );
    ylabel( 'v_L', 'Interpreter', 'none' );
    grid on
end

function [t, y_L, v_L, y_R, v_R] = load_motor_data()
    %path and file name of data
    fpath = './'; %path (change this!)
    fname_in = 'motor_calib.txt'; %file name (change this!)
    %load the motor calibration data
    motor_data = importdata([fpath,fname_in]);
    %unpack the motor calibration data
    t = motor_data(:,1);
    y_L = motor_data(:,2); v_L = motor_data(:,3);
    y_R = motor_data(:,4); v_R = motor_data(:,5);
    %plot the motor calibration data
    figure(1);
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