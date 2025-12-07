function SmithchartPlot(input,Gamma,Z0,ZL)
smithplot();
hold on

switch input
    case 'Gamma'
        plot(Gamma, 'ro', 'MarkerSize',10, 'LineWidth',2)
        

    case 'Load'
        zL = ZL/Z0;
        Gamma = (zL - 1) ./ (zL + 1);    % refleksionskoefficient
        smithplot(Gamma, 'ro', 'MarkerSize', 10, 'LineWidth', 2);

end

