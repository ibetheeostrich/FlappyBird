function [] = GraphGood()
    %% Plot parameters
    myorange        = [0.8500 0.3250 0.0980];
    myblue          = [0 0.4470 0.7410];
    myblack         = [0 0 0]/255;
    mygreen         = [0.4660 0.6740 0.1880];
    mycyan          = [0.3010 0.7450 0.9330];
    myyellow        = [0.9290 0.6940 0.1250];
    mypurple        = [0.4940 0.1840 0.5560];
    set(groot,'defaultAxesColorOrder',[mypurple;mygreen;myorange;myblue;myyellow;mycyan;myblack]);
    alw             = 1;                        % AxesLineWidth
    fsz             = 11;                       % Fontsize
    lw              = 1;                        % LineWidth
    msz             = 7;                       % MarkerSize
    set(0,'defaultLineLineWidth',lw);           % set the default line width to lw
    set(0,'defaultLineMarkerSize',msz);         % set the default line marker size to msz
    set(0,'defaultAxesLineWidth',alw);           % set the default line width to lw
    set(0,'defaultAxesFontSize',fsz);         % set the default line marker size to msz
    % set(0,'defaulttextInterpreter','latex')
    % set(0,'defaultAxesTickLabelInterpreter', 'latex');
    % set(0,'defaultLegendInterpreter', 'latex');
    set(0,'defaultFigureColor','w');
    set(0,'defaultAxesColor','w');

    %% How to plot
    % hFig            = figure();
    % hold on; grid on; grid minor; box on;
    % 
    % plot(x,y)
    % 
    % xlabel('xlabel ($s$)')
    % ylabel('ylabel ($\phi$)')
    % 
    % set(gca,'GridLineStyle','-')
    % set(gca,'MinorGridLineStyle','-')
    % set(gca,'GridColor','k')
    % set(gca,'MinorGridColor','k')
    % 
    % hleg = legend('$\phi$','$\theta$','$\psi$');
    % set(hleg,'EdgeColor',hleg.Color);
    % set(hleg,'Location','best');
    % set(hleg,'Interpreter','latex')
    
end