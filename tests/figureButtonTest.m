function [kIs, kIsnt, i] = figureButtonTest

kIs = 1;
kIsnt = 1;
fig = uifigure;
ax  = uiaxes('Parent', fig, 'Units', 'pixels', ...
    'Position', [50 80 400 250]);
buttonIs = uibutton(fig, 'push');
buttonIs.ButtonPushedFcn = @(btn,event) buttonIsPushed(btn,ax);
buttonIs.Position = [450 190 100 22];
buttonIs.Text = 'Yes is artifact';

buttonIsnt = uibutton(fig, 'push');
buttonIsnt.ButtonPushedFcn = @(btn,event) buttonIsntPushed(btn,ax);
buttonIsnt.Position = [450 160 100 22];
buttonIsnt.Text = 'No it isnt';

pause(.5)

i = 1;

waitfor(fig);

    function checkFinish()
        if i > 10
            close(fig)
        end
    end

    function updatePlotting(ax)
       title(ax, i) 
    end

    function buttonIsntPushed(btn,ax)
        kIsnt = kIsnt + 1;
        
        i = i + 1;
        updatePlotting(ax)
        checkFinish()
    end

    function buttonIsPushed(btn,ax)
        kIs = kIs + 1;
        i = i + 1;
        updatePlotting(ax)
        checkFinish
    end

end