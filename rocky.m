function main

end


function fitData
    % Fit: '3b'.
    [xData, yData] = prepareCurveData( t, v_L );
    
    % Set up fittype and options.
    ft = fittype( 'c*(1-exp(-a*x))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [1 0.0461713906311539];
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    
    % Plot fit with data.
    figure( 'Name', 'untitled fit 1' );
    h = plot( fitresult, xData, yData );
    legend( h, 'v_L vs. t', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 't', 'Interpreter', 'none' );
    ylabel( 'v_L', 'Interpreter', 'none' );
    grid on
end