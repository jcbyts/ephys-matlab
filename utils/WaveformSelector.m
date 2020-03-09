classdef WaveformSelector < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure              matlab.ui.Figure
        Figure                matlab.ui.Figure
        AxTime                matlab.graphics.axis.Axes
        AxEnergy              matlab.graphics.axis.Axes
        AxWaveforms           matlab.graphics.axis.Axes
        SelectPointsButton    matlab.ui.control.Button
        ClearSelectionButton  matlab.ui.control.Button
        SwitchGroupButton     matlab.ui.control.Button
        FinishSessionButton   matlab.ui.control.Button
        RemoveSelectedButton  matlab.ui.control.Button
        hEnUnselected
        hEnSelected
        hTimeUnselected
        hTimeSelected
        hWFUnselected
        hWFSelected
        SS
        WF
        EN
        NLEN
        selectedid
        unselectedid
    end
    
    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)
            
            opengl hardwarebasic
            
            % Create UIFigure
            app.UIFigure = uifigure;
            app.UIFigure.Position = [1100 200 200 300];
            app.UIFigure.Name = 'Controls';
            
            app.Figure = figure;
            app.Figure.Position = [100 100 1000 500];
            app.Figure.Name = 'Plotting';
            set(app.Figure,'renderer','opengl')

            % Create Axes for plotting
            
            % Timeline
            app.AxTime = axes(app.Figure);
            xlabel(app.AxTime, 'Time');
            ylabel(app.AxTime, 'Energy');
            app.AxTime.Position = [.1 .1 .85 .25];

            % Energy
            app.AxEnergy = axes(app.Figure);
            xlabel(app.AxEnergy, 'Energy');
            ylabel(app.AxEnergy, 'Nonlinear Energy');
            app.AxEnergy.Position = [.1 .5 .3 .4];
            
            % Waveforms
            app.AxWaveforms = axes(app.Figure);
            xlabel(app.AxWaveforms, 'Channel');
            ylabel(app.AxWaveforms, 'mV');
            app.AxWaveforms.Position = [.5 .5 .4 .4];
            
            hold(app.AxEnergy, 'on')
            hold(app.AxTime, 'on')
            hold(app.AxWaveforms, 'on')
            
            
            % Create SelectPointsButton
            app.SelectPointsButton = uibutton(app.UIFigure, 'push');
            app.SelectPointsButton.ButtonPushedFcn = createCallbackFcn(app, @SelectPointsButtonPushed, true);
            app.SelectPointsButton.Position = [50 190 100 22];
            app.SelectPointsButton.Text = 'Select Points';

            % Create ClearSelectionButton
            app.ClearSelectionButton = uibutton(app.UIFigure, 'push');
            app.ClearSelectionButton.ButtonPushedFcn = createCallbackFcn(app, @ClearSelectionButtonPushed, true);
            app.ClearSelectionButton.Position = [50 160 100 22];
            app.ClearSelectionButton.Text = 'Clear Selection';

            % Create SwitchGroupButton
            app.SwitchGroupButton = uibutton(app.UIFigure, 'push');
            app.SwitchGroupButton.ButtonPushedFcn = createCallbackFcn(app, @SwitchGroupButtonPushed, true);
            app.SwitchGroupButton.Position = [50 130 100 22];
            app.SwitchGroupButton.Text = 'Switch Group';

            % Create RemoveWaveformsButton
            app.RemoveSelectedButton = uibutton(app.UIFigure, 'push');
            app.RemoveSelectedButton.ButtonPushedFcn = createCallbackFcn(app, @RemoveSelectedButtonPushed, true);
            app.RemoveSelectedButton.Position = [50 100 100 22];
            app.RemoveSelectedButton.Text = 'Remove Selected';
            
            % Create FinishSessionButton
            app.FinishSessionButton = uibutton(app.UIFigure, 'push');
            app.FinishSessionButton.ButtonPushedFcn = createCallbackFcn(app, @FinishSessionButtonPushed, true);
            app.FinishSessionButton.Position = [50 70 100 22];
            app.FinishSessionButton.Text = 'Finish Session';
        end
    end

    methods (Access = private)
        function UpdatePlotting(app)
            
            if isempty(app.unselectedid)
                app.hEnUnselected.XData = nan;
                app.hEnUnselected.YData = nan;
                
                app.hTimeUnselected.XData = nan;
                app.hTimeUnselected.YData = nan;
                
            else
                app.hEnUnselected.XData = app.EN(app.unselectedid);
                app.hEnUnselected.YData = app.NLEN(app.unselectedid);
                
                app.hTimeUnselected.XData = app.SS(app.unselectedid);
                app.hTimeUnselected.YData = app.EN(app.unselectedid);
                
                mm = mean(app.WF(:,app.unselectedid),2);
                sd = std(app.WF(:,app.unselectedid),[],2);
                n = numel(mm);
                hold(app.AxWaveforms, 'off')
                app.hWFUnselected   = plot(app.AxWaveforms, 1:n, mm, 'k', 1:n, mm+sd, 'k--', 1:n, mm-sd, 'k--');
                hold(app.AxWaveforms, 'on')
            end
            
            if isempty(app.selectedid)
                app.hEnSelected.XData = nan;
                app.hEnSelected.YData = nan;
                
                app.hTimeSelected.XData = nan;
                app.hTimeSelected.YData = nan;
            else
                app.hEnSelected.XData   = app.EN(app.selectedid);
                app.hEnSelected.YData   = app.NLEN(app.selectedid);
                
                app.hTimeSelected.XData = app.SS(app.selectedid);
                app.hTimeSelected.YData = app.EN(app.selectedid);
                
                
                mm = mean(app.WF(:,app.unselectedid),2);
                sd = std(app.WF(:,app.unselectedid),[],2);
                n = numel(mm);
                hold(app.AxWaveforms, 'off')
                app.hWFUnselected   = plot(app.AxWaveforms, app.WF(:,app.selectedid));
            end
            
            
        end
        
        % Button pushed function: FinishSessionButton
        function FinishSessionButtonPushed(app, event)
            close(app.UIFigure, 'force')
            close(app.Figure,'force');
        end
        
        % Button pushed function: RemoveSelectedButton
        function RemoveSelectedButtonPushed(app, event)
            app.SS = app.SS(app.unselectedid);
            app.WF = app.WF(:,app.unselectedid);
            app.EN = app.EN(app.unselectedid);
            app.NLEN = app.NLEN(app.unselectedid);
            
            hold(app.AxWaveforms, 'off')
            app.hWFUnselected   = plot(app.AxWaveforms, app.WF, 'Color', .5*[1 1 1]);
            hold(app.AxWaveforms, 'on')
            
            app.selectedid = [];
            app.unselectedid = 1:numel(app.EN);
            
            app.UpdatePlotting();
        end

        % Button pushed function: SelectPointsButton
        function SelectPointsButtonPushed(app, event)
            axes(app.AxEnergy)
            title(app.AxEnergy, 'Click to create a polygon. Press enter to select', 'Fontweight', 'normal')
%             app.AxEnergy.Title.String = '
            [x,y] = ginput();
            
            in = inpolygon(app.EN,app.NLEN,x,y);
            app.selectedid = union(app.selectedid, find(in));
            app.unselectedid = setdiff(1:numel(app.EN), app.selectedid);
            
            app.UpdatePlotting();
        end

        % Button pushed function: ClearSelectionButton
        function ClearSelectionButtonPushed(app, event)
            app.selectedid = [];
            app.unselectedid = setdiff(1:numel(app.EN), app.selectedid);
            
            app.UpdatePlotting();
        end

        % Button pushed function: SwitchGroupButton
        function SwitchGroupButtonPushed(app, event)
            tmp = app.selectedid;
            app.selectedid = app.unselectedid;
            app.unselectedid = tmp;
            
            app.UpdatePlotting();
        end
    end

    

    methods (Access = public)

        % Construct app
        function app = WaveformSelector(SS,WF)
            
            % Create and configure components
            createComponents(app)
            
            app.SS   = SS;
            app.WF   = WF;
            app.NLEN = mean(abs((WF(2:end-1,:)-WF(1:end-2,:)) .* WF(2:end-1,:) .* WF(3:end,:)));
            app.EN   = mean(WF.^2);
            
            % plotting
            app.hEnUnselected   = plot(app.AxEnergy, app.EN, app.NLEN, '.k');
            app.hEnSelected     = plot(app.AxEnergy, nan, nan, '.r');
            app.hTimeUnselected = plot(app.AxTime, app.SS, app.EN, '.k');
            app.hTimeSelected   = plot(app.AxTime, nan, nan, '.r');
            mm = mean(app.WF,2);
            sd = std(WF,[],2);
            n = numel(mm);
            app.hWFUnselected   = plot(app.AxWaveforms, 1:n, mm, 'k', 1:n, mm+sd, 'k--', 1:n, mm-sd, 'k--');

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end