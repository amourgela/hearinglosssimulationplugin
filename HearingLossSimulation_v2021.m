%% Hearing Loss Simulation Plugin
% This plugin simulates 4 perceptual aspects of hearing loss
% High frequency attenuation, spectral smearing, temporal disruption &
% rapid loudness growth. The plugin offers customisability as well as 
% mute and bypass options for each ear.
% Developed by Angeliki Mourgela, Queen Mary University of London 
% Supported by EPSRC & BBC R&D
% For more details on the methods implemented in this code, refer to the
% the accompanying paper titled ''Investigation of a real-time hearing loss 
% simulation for use in audio production'', submitted at the 149th AES
% Convention.
%
%% 2021 Version Updates: 
% Rapid loudness growth is now utilising the envelope expansion method for
% simulation of recruitment by Nejime & Moore (1997).
%%
classdef HearingLossSimulation_v2021 <audioPlugin
properties
    SmearL            = 'Bypass';
    LosstypeL         = 'Bypass';
    LoudGrL           = false;
    SimulL            = false;
    muteL             = false;
    TempJitL          = false;
    
    SmearR            = 'Bypass';
    LosstypeR         = 'Bypass';
    LoudGrR           = false;
    SimulR            = false;
    muteR             = false;
    TempJitR          = false;
    

end 

properties (Access = private, Hidden)
    

    pSmearL    = 10;
    pLosstypeL = 3;
    pSmearR    = 10;
    pLosstypeR = 3;
    mildLossL  
    moderLossL 
    mildLossR  
    moderLossR 
    lpfilt
    FiltBankL
    FiltBankR
    stfL
    istfL
    stfR
    istfR
   
   
end
 properties (Constant) %interface setup 
 PluginInterface = audioPluginInterface(...
 'PluginName','Hearing Loss Simulation','VendorName','Angeliki Mourgela,Queen Mary University of London & BBC R&D',...
 'BackgroundColor',[1,1,1],'OutputChannels',2,...
audioPluginParameter('SimulL','Layout',[1,1],...
'DisplayName','Left Ear Bypass ','Mapping',{'enum','off', 'on'}),...
audioPluginParameter('SimulR','Layout',[1,3],...
'DisplayName','Right Ear Bypass ','Mapping',{'enum','off', 'on'}),...
        audioPluginParameter('muteL','Layout',[2,1],...
'DisplayName','Mute Left Ear','Mapping',{'enum', 'off', 'on'}),...
audioPluginParameter('muteR','Layout',[2,3],...
'DisplayName','Mute Right Ear','Mapping',{'enum', 'off', 'on'}),...
        audioPluginParameter('LosstypeL','Layout',[4,1],...
'DisplayName','High Frequency Attenuation','DisplayNameLocation','Above','Mapping',{'enum','Bypass', 'Mild','Moderate'}),...  
audioPluginParameter('LosstypeR','Layout',[4,3],...
'DisplayName','High Frequency Attenuation','DisplayNameLocation','Above','Mapping',{'enum','Bypass', 'Mild','Moderate'}),...  
        audioPluginParameter('SmearL','Layout',[6,1],...
'DisplayName','Smearing','DisplayNameLocation','Above','Mapping',{'enum','Bypass','Low','High'}),...
audioPluginParameter('SmearR','Layout',[6,3],...
'DisplayName','Smearing','DisplayNameLocation','Above','Mapping',{'enum','Bypass','Low','High'}),...    
            audioPluginParameter('LoudGrL','Layout',[8,1],...
'DisplayName','Rapid Growth','Mapping',{'enum', 'off', 'on'}),...
audioPluginParameter('LoudGrR','Layout',[8,3],...
'DisplayName','Rapid Growth','Mapping',{'enum', 'off', 'on'}),...
audioPluginParameter('TempJitL','Layout',[9,1],...
             'DisplayName','Temporal Disruption','Mapping',{'enum', 'off', 'on'}),...
audioPluginParameter('TempJitR','Layout',[9,3],...
             'DisplayName','Temporal Disruption','Mapping',{'enum', 'off', 'on'}),...
audioPluginGridLayout('RowHeight',[40,40,40,40,40,40,40,40,40],'ColumnWidth',[150,80,150]),'BackgroundImage','background.jpg');  
end
 
methods 

%% constructor section 

function plugin = HearingLossSimulation_v2021()
            WindowLength = 512;
            win =  hamming(WindowLength,'periodic');
            HopLength = 16;
            FFTLength =  WindowLength;
            
plugin.stfL  = dsp.STFT(win,WindowLength-HopLength,FFTLength);
plugin.istfL = dsp.ISTFT(win,WindowLength-HopLength,1,0);

plugin.stfR  = dsp.STFT(win,WindowLength-HopLength,FFTLength);
plugin.istfR = dsp.ISTFT(win,WindowLength-HopLength,1,0);

plugin.FiltBankL = gammatoneFilterBank([20 16000 ],32,getSampleRate(plugin));
plugin.FiltBankR = gammatoneFilterBank([20 16000 ],32,getSampleRate(plugin));

plugin.lpfilt   =  dsp.VariableBandwidthIIRFilter('FilterType','Lowpass','FilterOrder',10,...
    'PassbandFrequency',10);
plugin.mildLossL = multibandParametricEQ('NumEQBands',7,'HasHighShelfFilter',true,'HighShelfCutoff',8000,'HasLowpassFilter',...
    true,'LowpassCutoff',8000,'LowpassSlope',48,'Frequencies',[250,500,...
    1000,2000,4000,6000,8000],'QualityFactors',[2.3,2.3,2.3,2.3,2.3,2.3,2.3],'PeakGains',[0 0 -15 -15 -20 -20 -25]);
plugin.moderLossL = multibandParametricEQ('NumEQBands',7,'HasHighShelfFilter',true,'HighShelfCutoff',8000,'HasLowpassFilter',...
    true,'LowpassCutoff',8000,'LowpassSlope',48,'Frequencies',[250,500,...
    1000,2000,4000,6000,8000],'QualityFactors',[2.3,2.3,2.3,2.3,2.3,2.3,2.3],'PeakGains',[0 0 -20 -20 -30 -30 -45]);
plugin.mildLossR = multibandParametricEQ('NumEQBands',7,'HasHighShelfFilter',true,'HighShelfCutoff',8000,'HasLowpassFilter',...
    true,'LowpassCutoff',8000,'LowpassSlope',48,'Frequencies',[250,500,...
    1000,2000,4000,6000,8000],'QualityFactors',[2.3,2.3,2.3,2.3,2.3,2.3,2.3],'PeakGains',[0 0 -15 -15 -20 -20 -25]);
plugin.moderLossR = multibandParametricEQ('NumEQBands',7,'HasHighShelfFilter',true,'HighShelfCutoff',8000,'HasLowpassFilter',...
    true,'LowpassCutoff',8000,'LowpassSlope',48,'Frequencies',[250,500,...
    1000,2000,4000,6000,8000],'QualityFactors',[2.3,2.3,2.3,2.3,2.3,2.3,2.3],'PeakGains',[0 0 -20 -20 -30 -30 -45]);
end

%% setter & getter functions for tuneable properties
  function set.SmearL(plugin,val)
            validatestring(val,{'Low','High','Bypass'},...
                'set.SmearL','SmearL');
            if strcmp(val,'Bypass')
                plugin.pSmearL = 10;
            elseif strcmp(val,'Low')
                plugin.pSmearL = 100;
           elseif strcmp(val,'High')
                plugin.pSmearL = 200;
            end
         
         smearL = 10;
         
         switch plugin.pSmearL %#ok<*MCSUP>
                case 10
                   smearL = 10;
                case 100
                   smearL = 100;
                case 200
                   smearL = 200;
         end
         
         plugin.lpfilt.PassbandFrequency = smearL;

        end
 
           function val = get.SmearL(plugin)
            x = plugin.pSmearL;
            if     x == 10
                val = 'Bypass';
            elseif x ==100
                val = 'Low';
            elseif x ==200
                val = 'High';
            else
                val = 'Bypass';
        
           end
           end
            
             function set.LosstypeL(plugin,val) 
            validatestring(val,{'Mild','Moderate','Bypass'},...
                'set.LosstypeL','LosstypeL');
            if strcmp(val,'Bypass')
                plugin.pLosstypeL = 3;
           elseif strcmp(val,'Mild')
                plugin.pLosstypeL = 1;
           elseif strcmp(val,'Moderate')
                plugin.pLosstypeL = 2;
           end
             end

             
            function val = get.LosstypeL(plugin) 
            x = plugin.pLosstypeL;
            if     x == 3
                val = 'Bypass';
            elseif x == 2
                val = 'Moderate';
            elseif x == 1
                val = 'Mild';
            else
                val = 'Bypass';
             end
            end
           

function set.SmearR(plugin,val)
            validatestring(val,{'Low','High','Bypass'},...
                'set.SmearR','SmearR');
            if strcmp(val,'Bypass')
                plugin.pSmearR = 10;
            elseif strcmp(val,'Low')
                plugin.pSmearR = 100;
           elseif strcmp(val,'High')
                plugin.pSmearR = 200;
            end
         
         smearR = 10;
         
         switch plugin.pSmearR %#ok<*MCSUP>
                case 10
                   smearR = 10;
                case 100
                   smearR = 100;
                case 200
                   smearR = 200;
         end
         
         plugin.lpfilt.PassbandFrequency = smearR;

        end
 
   function val = get.SmearR(plugin)
            x = plugin.pSmearR;
            if     x == 10
                val = 'Bypass';
            elseif x ==100
                val = 'Low';
            elseif x ==200
                val = 'High';
            else
                val = 'Bypass';
        
           end
    end
            
   function set.LosstypeR(plugin,val) 
            validatestring(val,{'Mild','Moderate','Bypass'},...
                'set.LosstypeR','LosstypeR');
            if strcmp(val,'Bypass')
                plugin.pLosstypeR = 3;
           elseif strcmp(val,'Mild')
                plugin.pLosstypeR = 1;
           elseif strcmp(val,'Moderate')
                plugin.pLosstypeR = 2;
           end
    end

             
   function val = get.LosstypeR(plugin) 
            x = plugin.pLosstypeR;
            if     x == 3
                val = 'Bypass';
            elseif x == 2
                val = 'Moderate';
            elseif x == 1
                val = 'Mild';
            else
                val = 'Bypass';
            end
    end
           
         
 %% main processing          
function out = process(plugin,in)
% Separate left and right channels 

     in_L = in(:,1);
     in_R = in(:,2);
% Create noise signal for smearing      
fs = getSampleRate(plugin);

rng("default");
numOfSamples  = length(in);
xnoise        = randn(numOfSamples,21);
ynoise1       = plugin.lpfilt(xnoise);
ynoise1       = 0.5*ynoise1; %scale it for subtlety


%% Gammatone Filterbank Separation 


hassmearL = plugin.pSmearL;
hassmearR = plugin.pSmearR;


FiltSignal_L = step(plugin.FiltBankL,in_L);
FiltSignal_R = step(plugin.FiltBankR,in_R);

outsm_L_low1 = FiltSignal_L(:,(1:11));
outsm_L_high1= FiltSignal_L(:,(12:32));
outsm_R_low1 = FiltSignal_R(:,(1:11));
outsm_R_high1= FiltSignal_R(:,(12:32));


%% Apply Spectral Smearing 

%left channel 
outsm_L_highS = zeros(size(outsm_L_high1));
for kk = 1:21
if hassmearL ==10
    outsm_L_highS(:,kk) = outsm_L_high1(:,kk);
else   
    outsm_L_highS(:,kk) = 20*outsm_L_high1(:,kk).*ynoise1(:,kk);  %scaled to compensate for loss 
end
end

outsm_L_highSE = sum(outsm_L_highS,2);

 %right channel 
outsm_R_highS = zeros(size(outsm_R_high1));

for gg = 1:21

if hassmearR ==10
    outsm_R_highS(:,gg) = outsm_R_high1(:,gg);
else   
    outsm_R_highS(:,gg) = 20*outsm_R_high1(:,gg).*ynoise1(:,gg);  %scaled to compensate for loss 
end
end

outsm_R_highSE = sum(outsm_R_highS,2);
   
%% Apply High Frequency Attenuation 
%left channel
     
     if plugin.pLosstypeL     == 1
        outAud_L = step(plugin.mildLossL,outsm_L_highSE);
     elseif plugin.pLosstypeL == 2
        outAud_L = step(plugin.moderLossL,outsm_L_highSE);
     elseif plugin.pLosstypeL == 3
        outAud_L = outsm_L_highSE ;
     else 
        outAud_L = outsm_L_highSE ;
     end 
%right channel 
     
     if plugin.pLosstypeR     == 1
        outAud_R = step(plugin.mildLossR,outsm_R_highSE);
     elseif plugin.pLosstypeR == 2
        outAud_R = step(plugin.moderLossR,outsm_R_highSE);
     elseif plugin.pLosstypeR == 3
        outAud_R = outsm_R_highSE;
     else 
        outAud_R = outsm_R_highSE;
     end 


%% Apply Loudness Growth 
%left channel 


if  plugin.LoudGrL==true
   outsm_L_high2 = loudrecruit(outAud_L);

elseif plugin.LoudGrL==false
   outsm_L_high2 = outAud_L;
else
   outsm_L_high2 = outAud_L;
end


if  plugin.LoudGrR==true
    outsm_R_high2 = loudrecruit(outAud_R);
 elseif  plugin.LoudGrR==false
    outsm_R_high2 = outAud_R;
else
   outsm_R_high2 = outAud_R;
end



%% Apply Temporal disruption 

FrameLengthL = length(in_L);
HopLengthL = 16;
numHopsPerFrameL = FrameLengthL/HopLengthL;

%sum the lower portion of the gammatone 

outsm_L_low2 = sum(outsm_L_low1,2);
outsm_R_low2 = sum(outsm_R_low1,2);


if plugin.TempJitL==true 
    
 outsm_L_low3 = zeros(size(in));

for indexL = 1:numHopsPerFrameL        
       
X_L = plugin.stfL(outsm_L_low2((indexL-1)*HopLengthL+1:indexL*HopLengthL));
 
 [L_L,N_L]=size(X_L);
%  
magL = abs(X_L);%original magnitude
phiL = angle(X_L);
rnd_thetaL = -pi/2+pi/2*rand(L_L,N_L); %random phase shifts
mix_thetaL = phiL-rnd_thetaL; 
XL_rand = magL.*exp(1i*mix_thetaL); 
XL_rand = sum(XL_rand,2);
% Convert back to time-domain
outsm_L_low3((indexL-1)*HopLengthL+1:indexL*HopLengthL) = plugin.istfL(XL_rand);        
end    
elseif plugin.TempJitL==false
outsm_L_low3 = outsm_L_low2;
else
outsm_L_low3 = outsm_L_low2;
end

% right channel 

FrameLengthR = length(in_R);
HopLengthR = 16;
numHopsPerFrameR = FrameLengthR/HopLengthR;

if plugin.TempJitR==true 
    
 outsm_R_low3 = zeros(size(in_R));

for indexR = 1:numHopsPerFrameR        
       
X_R = plugin.stfR(outsm_R_low2((indexR-1)*HopLengthR+1:indexR*HopLengthR));
 
 [L_R,N_R]=size(X_R);
%  
magR = abs(X_R);%original magnitude
phiR = angle(X_R);
rnd_thetaR = -pi/2+pi/2*rand(L_R,N_R);%random phase shifts
mix_thetaR = phiR-rnd_thetaR; 
XR_rand = magR.*exp(1i*mix_thetaR); 
XR_rand = sum(XR_rand,2);
% Convert back to time-domain
outsm_R_low3((indexR-1)*HopLengthR+1:indexR*HopLengthR) = plugin.istfR(XR_rand);        
end    
elseif plugin.TempJitR==false
outsm_R_low3 = outsm_R_low2;
else
outsm_R_low3 = outsm_R_low2;
end


%% Sum the channels & Check for Mute and Bypass
LChan = outsm_L_high2+outsm_L_low3;
RChan = outsm_R_high2+outsm_R_low3;

sumL2 = sum(LChan,2);
sumR2 = sum(RChan,2);

% Check for Bypass
if plugin.SimulL==true  
    sumL=in_L;
else
    sumL = sumL2;
end

if plugin.SimulR==true
    sumR=in_R;
else
    sumR = sumR2;
end

%Check for Mute
    if plugin.muteL==true && plugin.muteR==false
out=cat(2,zeros(size(sumL)),sumR);
    elseif plugin.muteR==true && plugin.muteL==false
out=cat(2,sumL,zeros(size(sumR)));
    elseif plugin.muteR==true && plugin.muteL==true
out=zeros(size(in));     
    else
out=cat(2,sumL,sumR);
    end
end
function  reset(plugin)
      
            plugin.FiltBankL.SampleRate = getSampleRate(plugin); 
            plugin.FiltBankR.SampleRate = getSampleRate(plugin);      
            reset(plugin.lpfilt);
            reset(plugin.FiltBankL);
            reset(plugin.FiltBankR);
            reset(plugin.mildLossL);
            reset(plugin.moderLossL);
            reset(plugin.mildLossR);
            reset(plugin.moderLossR);
            reset(plugin.stfL)
            reset(plugin.istfL)
            reset(plugin.stfR)
            reset(plugin.istfR)

end
end
end
