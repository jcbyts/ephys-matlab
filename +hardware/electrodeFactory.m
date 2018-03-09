function shank = electrodeFactory(name)
% This builds the probe files that we use to handle channel maps and
% filtering corrections

if nargin == 0
    shankList = {'Shank2', ...
        'AcuteZif2Om32', ...
        'MTSingleElectrode', ...
        'MTSingleV132', ...
        'Shank3', ...
        'Shank4', ...
        'Shank5', ...
        'Single202328', ...
        };
    
    name = io.selectFromList(shankList);
    if iscell(name)
        name = name{1};
    end
end

switch name
    case 'Shank2'
        
        shank{1} = hardware.electrode.Shank2;
        shank{1}.headstages{1} = hardware.headstage.intan_RHD2132;
        shank{1}.name = 'V1';
        
    case 'AcuteZif2Om32'
        
        shank{1} = hardware.electrode.AtlasZifOmnDrive_1;
        shank{1}.headstages{1} = hardware.headstage.intan_RHD2132;
        shank{1}.name = name;
        
    case 'MTSingleElectrode'
        
        chanMap = [8 25];
        shank{1} = hardware.electrode.customChannelMap(chanMap);
        shank{1}.name = 'MtBurrHoleMapping';
        
    case 'MTSingleV132'
        
        shank{1} = hardware.electrode.Shank2;
        shank{1}.headstages{1} = hardware.headstage.intan_RHD2132;
        shank{1}.name = 'V1';
        
        chanMap = 4;
        shank{2} = hardware.electrode.customChannelMap(chanMap);
        shank{2}.name = 'MtBurrHoleMapping';
        
    case 'Shank3'
        
        shank{1} = hardware.electrode.Shank2;
        shank{1}.headstages{1} = hardware.headstage.intan_RHD2132;
        shank{1}.name = 'MT';
        
    case 'Shank4'
        
        shank{1} = hardware.electrode.Shank2;
        shank{1}.headstages{1} = hardware.headstage.intan_RHD2132;
        shank{1}.name = 'MT';
        
        shank{2} = hardware.electrode.Shank2;
        shank{2}.headstages{1} = hardware.headstage.intan_RHD2132;
        shank{2}.name = 'MT';
        
    case 'Shank5'
        
        shank{1} = hardware.electrode.Shank2;
        shank{1}.headstages{1} = hardware.headstage.intan_RHD2132;
        shank{1}.name = name;
        
    case 'Single202328'
        chanMap = [20 23 28];
        shank{1} = hardware.electrode.customChannelMap(chanMap);
        shank{1}.name = 'Single202328';
        
        
end