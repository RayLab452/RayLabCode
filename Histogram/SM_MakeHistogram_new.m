function varargout = SM_MakeHistogram_new(varargin)
% SM_MAKEHISTOGRAM_NEW MATLAB code for SM_MakeHistogram_new.fig
%      SM_MAKEHISTOGRAM_NEW, by itself, creates a new SM_MAKEHISTOGRAM_NEW or raises the existing
%      singleton*.
%
%      H = SM_MAKEHISTOGRAM_NEW returns the handle to a new SM_MAKEHISTOGRAM_NEW or the handle to
%      the existing singleton*.
%
%      SM_MAKEHISTOGRAM_NEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SM_MAKEHISTOGRAM_NEW.M with the given input arguments.
%
%      SM_MAKEHISTOGRAM_NEW('Property','Value',...) creates a new SM_MAKEHISTOGRAM_NEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SM_MakeHistogram_new_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SM_MakeHistogram_new_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SM_MakeHistogram_new

% Last Modified by GUIDE v2.5 06-Jun-2022 17:30:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SM_MakeHistogram_new_OpeningFcn, ...
                   'gui_OutputFcn',  @SM_MakeHistogram_new_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before SM_MakeHistogram_new is made visible.
function SM_MakeHistogram_new_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SM_MakeHistogram_new (see VARARGIN)

% Choose default command line output for SM_MakeHistogram_new
handles.output = hObject;
handles.histogramData=[];
handles.fitFlag=0;
handles.perMolHist=[];      % one normalized histogram per molecule (rows)
handles.bobaResults=[];     % output of bobaFRET
handles.spreadResults=[];   % output of moleculeSpread
handles.bobaLabels=[];          % text object handles for the state-population labels
handles.bobaLabelDefaultPos=[]; % [x y] each label sits at with zero manual offset
handles.bobaLabelOffsets=[];    % [dx dy] manual nudge per state, from dragging; label-only
handles.dragLabelIndex=[];      % index of the label currently being dragged (empty = none)
handles.bobaRanges=[];          % patch handles for the per-state shaded ranges
handles.bobaRangeBorders=[];    % line handles for the range border lines

% Bootstrap controls, created here so the .fig file does not need to be edited.
% They sit under the Save button (positions in characters, like the .fig).
figColor = get(hObject,'Color');
handles.nBootLabel = uicontrol(hObject,'Style','text','Units','characters', ...
    'Position',[127 22.5 21 1.1],'String','Bootstrap replicates:', ...
    'HorizontalAlignment','left','BackgroundColor',figColor);
handles.nBoot = uicontrol(hObject,'Style','edit','Units','characters', ...
    'Position',[127 20.4 10 1.7],'String','1000','BackgroundColor','white');
handles.bobaButton = uicontrol(hObject,'Style','pushbutton','Units','characters', ...
    'Position',[149.8 20.4 25 4.7],'String','BOBA bootstrap','FontSize',12, ...
    'Callback',@(src,evt) boba_Callback(src,evt,guidata(src)));
handles.spreadButton = uicontrol(hObject,'Style','pushbutton','Units','characters', ...
    'Position',[149.8 14.4 25 4.7],'String','Spread (molecules)','FontSize',12, ...
    'Callback',@(src,evt) spread_Callback(src,evt,guidata(src)));
handles.showErrorBars = uicontrol(hObject,'Style','checkbox','Units','characters', ...
    'Position',[127 17.4 21 1.7],'String','Show error bars','Value',1, ...
    'BackgroundColor',figColor, ...
    'Callback',@(src,evt) showErrorBars_Callback(src,evt,guidata(src)));
handles.showShading = uicontrol(hObject,'Style','checkbox','Units','characters', ...
    'Position',[109 17.4 17 1.7],'String','Shade band','Value',1, ...
    'BackgroundColor',figColor, ...
    'Callback',@(src,evt) showErrorBars_Callback(src,evt,guidata(src)));
handles.showLabels = uicontrol(hObject,'Style','checkbox','Units','characters', ...
    'Position',[127 15.7 21 1.7],'String','Show state labels','Value',1, ...
    'BackgroundColor',figColor, ...
    'Callback',@(src,evt) showErrorBars_Callback(src,evt,guidata(src)));
% Friendlier front end for the "Range (+/-sigma)" box further down (rangeSigmaMultiplier):
% pick a preset here instead of typing a multiplier by hand. Picking a preset writes
% the corresponding number into that box and redraws; typing directly into that box
% still works too, and switches this popup to "Custom" (see rangeSigmaMultiplier_Callback).
handles.rangeSigmaPreset = uicontrol(hObject,'Style','popupmenu','Units','characters', ...
    'Position',[109 15.7 17 1.7], ...
    'String',{'1 sigma (~68%)','1.5 sigma (~87%)','2 sigma (~95%)','3 sigma (~99.7%)','Custom'}, ...
    'Value',1,'BackgroundColor','white', ...
    'Callback',@(src,evt) rangeSigmaPreset_Callback(src,evt,guidata(src)));
handles.resetLabelsButton = uicontrol(hObject,'Style','pushbutton','Units','characters', ...
    'Position',[127 13.2 21 1.9],'String','Reset label positions','FontSize',9, ...
    'Callback',@(src,evt) resetLabels_Callback(src,evt,guidata(src)));
handles.showRanges = uicontrol(hObject,'Style','checkbox','Units','characters', ...
    'Position',[127 11.0 21 1.7],'String','Show state ranges','Value',1, ...
    'BackgroundColor',figColor, ...
    'Callback',@(src,evt) showErrorBars_Callback(src,evt,guidata(src)));

% The "sigma:" box above, and every fitted width in this program, is the
% Gaussian's standard deviation -- gaussian() below and bobaFRET.m's
% gaussianSigma() are both parameterized directly in sigma, with no FWHM
% conversion anywhere in the fitting math.
%
% How many standard deviations the shaded state-range overlay (drawBobaRanges)
% extends on each side of a state's center. 1 = the classic +/-1 sigma band;
% 1.96 ~= a 95% interval; etc. Positions here may need a small nudge once you
% see them rendered -- this corner of the panel is already fairly crowded.
handles.rangeSigmaLabel = uicontrol(hObject,'Style','text','Units','characters', ...
    'Position',[127 2.0 12 1.7],'String','Range (+/-sigma):', ...
    'HorizontalAlignment','left','BackgroundColor',figColor);
handles.rangeSigmaMultiplier = uicontrol(hObject,'Style','edit','Units','characters', ...
    'Position',[140 2.0 8 1.7],'String','1','BackgroundColor','white', ...
    'Callback',@(src,evt) rangeSigmaMultiplier_Callback(src,evt,guidata(src)));

% States below this population (as a %, of the FULL-DATA fit) are treated as
% noise/non-states and left out of the result: not shown as a label or shaded
% range, not counted in the BOBA bootstrap result popup, and not written to the
% _BOBA.csv file. Default 0.5 means a state fitted at under 0.5% population is
% ignored. This never changes the fit itself (all K Gaussians are still fit and
% still count toward every other state's population) -- it only filters what
% gets reported afterward. Same row as the range-sigma box, reusing the same
% right edge the button column already established.
handles.minPopPercentLabel = uicontrol(hObject,'Style','text','Units','characters', ...
    'Position',[149.8 2.0 16 1.7],'String','Min population (%):', ...
    'HorizontalAlignment','left','BackgroundColor',figColor);
handles.minPopPercent = uicontrol(hObject,'Style','edit','Units','characters', ...
    'Position',[167 2.0 7 1.7],'String','0.5','BackgroundColor','white');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SM_MakeHistogram_new wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SM_MakeHistogram_new_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

handles.workingDir=pwd; % it prints the path of program directory
handles.pathName =[];
handles.fileName=[];
handles.faceColor=[0.9020 0.9020 0.9020];
handles.edgeColor=[0.6510 0.6510 0.6510];
handles.fitColor=[1 0 0];
handles.statement='';
varargout{1} = handles.output;
guidata(hObject, handles);


%% User active buttons...
% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
workingDir="C:\data"; % where program is there
[filelist pathName fi] = uigetfile('*_forvBFRET.dat','Select files on which you wish to perform drift correction.','MultiSelect','on', workingDir);
if iscell(filelist) == 0
    filelist2{1} = filelist;
    filelist = filelist2;
    clear filelist2;
end
handles.pathName=pathName;
handles.filelist=filelist;
guidata(hObject, handles);
[handles.histogramData, handles.perMolHist]=CalculateHistogram(hObject, eventdata, handles);
handles.bobaResults=[];
handles.spreadResults=[];
handles.bobaLabels=[];
handles.bobaLabelDefaultPos=[];
handles.bobaLabelOffsets=[];
handles.fitFlag=0;
guidata(hObject, handles);
plotHistogram(handles.histogramPlot,hObject, eventdata, handles)

% --- Executes on button press in rePlot.
function rePlot_Callback(hObject, eventdata, handles)
redrawPlot(hObject, eventdata, handles);   % keeps error bars / labels / ranges after Replot, not just after BOBA




function [R2_value,parameter]=fitGaussMath(hObject, eventdata, handles)
centers=handles.histogramData(1,:);
norm_T_histCount=handles.histogramData(2,:);
in_Centers=str2double(split(get(handles.in_Centers,'String'),","));
in_sigma=str2double(split(get(handles.in_sigma,'String'),","));   % standard deviation -- gaussian() below takes sigma directly
initialGuesses=[in_Centers in_sigma];
startingGuesses = reshape(initialGuesses', 1, []);
%%%%%%%%%%%%%%%%%%%%%%%%%%
global  c NumTrials TrialError
% 	warning off

% Initializations
NumTrials = 0;  % Track trials
TrialError = 0; % Track errors
% t and y must be row vectors.
tFit = reshape(centers, 1, []);
y = reshape(norm_T_histCount, 1, []);
%-------------------------------------------------------------------------------------------------------------------------------------------
% Perform an iterative fit using the FMINSEARCH function to optimize the height, width and center of the multiple Gaussians.
options = optimset('TolX', 1e-4, 'MaxFunEvals', 10^12);  % Determines how close the model must fit the data
% First, set some options for fminsearch().
options.TolFun = 1e-4;
options.TolX = 1e-4;
options.MaxIter = 100000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEAVY LIFTING DONE RIGHT HERE:
% Run optimization
[parameter, fval, flag, output] = fminsearch(@(lambda)(fitgauss(lambda, tFit, y)), startingGuesses, options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate residuals
means = parameter(1 : 2 : end);
widths = parameter(2 : 2 : end);
numGaussians = length(c);
thisEstimatedCurve=zeros(length(tFit),numGaussians);
for k = 1 : numGaussians
    % Get each curve.
    thisEstimatedCurve(:,k) = c(k) .* gaussian(tFit, means(k), widths(k));
end
% Overall curve estimate is the sum of the component curves.
yhat = sum(thisEstimatedCurve,2);

% residual sum of squares:SS_res=sum((y-yhat)^2)
SS_res= sum((y-yhat').^2);
SS_tot=sum((y-mean(y)).^2);
R2_value = 1-SS_res/SS_tot;
guidata(hObject, handles);
%% -------------------------------------

function CurveFitData=fitGaussPlots(figureHandle,R2_value,parameter,hObject, eventdata, handles)
global  c
numGaussians = length(c);

centers=handles.histogramData(1,:);
tFit = reshape(centers, 1, []);
means = parameter(1 : 2 : end);
widths = parameter(2 : 2 : end);

%% Now plot results.

tFit2 = linspace(tFit(1),tFit(end),1000);
thisEstimatedCurve=zeros(length(tFit2),numGaussians);
CM = jet(numGaussians);
	for k = 1 : numGaussians
		% Get each curve.
		thisEstimatedCurve(:,k) = c(k) .* gaussian(tFit2, means(k), widths(k));
        hold on;
        plot(figureHandle,tFit2, thisEstimatedCurve(:,k),'color',CM(k,:), 'LineWidth', str2double(get(handles.fitThickness,'String'))-1);
	    statement{k}= strcat('peak ', num2str(k),': ',num2str(means(k),1),char(177),num2str(widths(k),1),' \sigma' );
    end
% Overall curve estimate is the sum of the component curves.
yhat2 = sum(thisEstimatedCurve,2);

hold on;
plot(figureHandle, tFit2, yhat2, 'color',handles.fitColor, 'LineWidth', str2double(get(handles.fitThickness,'String')));
hold off;
%%%%%%%%%%%%%%%%%%%%%
uistack(figureHandle, 'top')


CurveFitData=[tFit2' thisEstimatedCurve yhat2];

handles.statement=[strcat('N=',num2str(length(handles.filelist))),statement,{strcat('R2 value: ', num2str(R2_value))}];

textPos=str2double(split(get(handles.textPos,'String'),","));
handles.text=text(figureHandle,textPos(1),textPos(2),handles.statement,'FontSize',str2double(get(handles.axisSize,'String'))-1);
guidata(hObject, handles);



% --- Executes on button press in fitGauss.
function fitGauss_Callback(hObject, eventdata, handles)

plotHistogram(handles.histogramPlot,hObject, eventdata, handles);
[R2_value,parameter]=fitGaussMath(hObject, eventdata, handles);
handles.R2_value=R2_value;
handles.parameter=parameter;
CurveFitData=fitGaussPlots(handles.histogramPlot,R2_value,parameter,hObject, eventdata, handles);
handles.CurveFitData=CurveFitData;
handles.fitFlag=1;
guidata(hObject, handles);

















% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
    dirForSavingHistogram='Histogram';
    % handles.paht comes with '\' at the end
    
    dest_folder=[handles.pathName,dirForSavingHistogram];
    if exist (dest_folder,'dir');
    else mkdir(handles.pathName,dirForSavingHistogram);
    end
    %-----------------------------------------------------
    NamePrefix=get(handles.plotTitle, 'String');
    NamePrefix = regexprep(NamePrefix,'[ ]','');
    NamePrefix = regexprep(NamePrefix,'[,.;:/\~!@#$%^&*()_+]','_');
    
    fName_Picture=[dest_folder '\' NamePrefix];

%% redraw the entire figure again! Stupid!! but that is how you do it.
F=figure();

plotHandle = axes('Parent',F);
plotHistogram(plotHandle,hObject, eventdata, handles);

% Same overlays shown live in the GUI (redrawPlot), baked into the saved image too.
% ranges go on first so the fit curve below draws on top of the shading, and
% persistHandles=false keeps this from touching the live GUI's label/range state.
if ~isempty(handles.bobaResults) && get(handles.showRanges,'Value')==1
    drawBobaRanges(hObject, eventdata, handles, plotHandle, false);
end

 if handles.fitFlag==1
   CurveFitData=fitGaussPlots(plotHandle,handles.R2_value,handles.parameter,hObject, eventdata, handles);
   save(strcat(fName_Picture,'_CurveFitData.dat'), 'CurveFitData', '-ascii');
 end

if ~isempty(handles.bobaResults) && get(handles.showLabels,'Value')==1
    drawBobaLabels(hObject, eventdata, handles, plotHandle, false);
end

 box off;
 %% Now plot results.
saveas(gcf,strcat(fName_Picture,'.jpg'));
saveas(gcf,strcat(fName_Picture,'.eps'));


HistData=handles.histogramData';



save(strcat(fName_Picture,'_Histdata.dat'), 'HistData', '-ascii');

if ~isempty(handles.bobaResults)
    writeBobaSummary(strcat(fName_Picture,'_BOBA.csv'), handles.bobaResults, handles);
end
if ~isempty(handles.spreadResults)
    writeSpreadSummary(strcat(fName_Picture,'_MoleculeSpread.csv'), handles.spreadResults);
end
if ~isempty(handles.bobaResults)
    writeLabelOffsets(strcat(fName_Picture,'_LabelPositions.csv'), handles.bobaLabelOffsets);
end

 close(F);




















%% Bootstrap analysis (uses bobaFRET.m and moleculeSpread.m, which must be in the same folder)
function boba_Callback(hObject, eventdata, handles)
% BOBA: bootstrap over molecules of the Gaussian fit to the mean histogram.
if isempty(handles.perMolHist)
    errordlg('Load data first (Load button).','BOBA bootstrap'); return
end
nBoot=round(str2double(get(handles.nBoot,'String')));
if isnan(nBoot) || nBoot<10
    errordlg('Bootstrap replicates must be a number of at least 10.','BOBA bootstrap'); return
end
% Same starting guesses as the Fit Gaussian button: [center1 sigma1 center2 sigma2 ...]
% (bobaFRET.m, like gaussian() above, is parameterized directly in sigma.)
in_Centers=str2double(split(get(handles.in_Centers,'String'),","));
in_sigma=str2double(split(get(handles.in_sigma,'String'),","));
if numel(in_Centers)~=numel(in_sigma) || any(isnan([in_Centers; in_sigma]))
    errordlg('Centers and sigma must be lists of numbers of equal length.','BOBA bootstrap'); return
end
startGuess=reshape([in_Centers in_sigma]', 1, []);
% If Fit Gaussian has been run, start BOBA from ITS solution (the red curve) instead of the
% raw guess boxes, so the bootstrap's full-data fit begins where the GUI fit ended up.
% bobaFRET sorts states by center itself, so the order here doesn't matter. Falls back to
% the guess boxes if no fit exists yet, or if the number of states no longer matches
% (e.g. the boxes were edited after fitting).
if handles.fitFlag==1 && numel(handles.parameter)==numel(startGuess)
    startGuess=handles.parameter(:)';
    startGuess(2:2:end)=abs(startGuess(2:2:end));   % widths are sigmas, must be positive
end
set(hObject,'Enable','off','String','Running...'); drawnow;
restoreButton=onCleanup(@() set(hObject,'Enable','on','String','BOBA bootstrap'));
try
    handles.bobaResults=bobaFRET(handles.perMolHist, handles.histogramData(1,:), startGuess, nBoot);
catch ME
    errordlg(ME.message,'BOBA bootstrap'); return
end
guidata(hObject, handles);
redrawPlot(hObject, eventdata, handles);   % show the error bars on the histogram
r=handles.bobaResults; s=r.summary;
keep=statesAboveMinPop(r.full.pop, handles);   % states below "Min population (%)" are left out of this report
msg={sprintf('%d molecules, %d/%d replicates converged, fit R^2 = %.3f', r.nMolecules, r.boot.nValid, r.boot.nBoot, r.full.R2)};
rangeMult = str2double(get(handles.rangeSigmaMultiplier,'String'));
if any(~keep)
    msg{end+1}=sprintf('(%d state(s) below %s%% population hidden)', sum(~keep), get(handles.minPopPercent,'String'));
end
if isnan(rangeMult) || rangeMult <=0
    rangeMult = 1; %same fallback as drawBobaRanges
end
for k=1:numel(s.popMean)
    if ~keep(k), continue; end
    mu = r.full.mu(k);  sg = r.full.sigma(k);
    msg{end+1}=sprintf('State %d: center %.3f +/- %.3f, population %.3f +/- %.3f (%g sigma)', k, s.muMean(k), s.muStd(k), s.popMean(k), s.popStd(k));
    % 1) sigma of the fitted Gaussian itself
    msg{end+1}=sprintf('     95%% range (center +/- 1.96 sigma): %.3f to %.3f', mu-rangeMult*sg, mu+rangeMult*sg);
    % 2) bootstrap sigma CI on the center
    m = r.boot.mu(:,k);  m = sort(m(~isnan(m)));
    if numel(m) >= 20
        ci = m([max(1,round(0.025*numel(m))), min(numel(m),round(0.975*numel(m)))]);
        msg{end+1}=sprintf('     center 95%% CI (bootstrap): %.3f to %.3f', ci(1), ci(2));
    end
end
msgbox(msg,'BOBA bootstrap result');

function spread_Callback(hObject, eventdata, handles)
% Spread of per-molecule FRET across molecules, with bootstrap uncertainty.
if isempty(handles.perMolHist)
    errordlg('Load data first (Load button).','Spread'); return
end
nBoot=round(str2double(get(handles.nBoot,'String')));
if isnan(nBoot) || nBoot<10
    errordlg('Bootstrap replicates must be a number of at least 10.','Spread'); return
end
paths=cellfun(@(f) [handles.pathName f], handles.filelist, 'UniformOutput', false);
try
    r=moleculeSpread(paths, nBoot, [], [-0.2 1.2], false);   % same window as the histogram
catch ME
    errordlg(ME.message,'Spread'); return
end
handles.spreadResults=r;
guidata(hObject, handles);
msg={sprintf('%d molecules, %d replicates', r.median.n, nBoot), ...
     sprintf('Median FRET: SD across molecules %.3f (95%% interval %.3f to %.3f)', r.median.sd, r.median.sdCI95(1), r.median.sdCI95(2)), ...
     sprintf('In-range mean FRET: SD across molecules %.3f (95%% interval %.3f to %.3f)', r.mean.sd, r.mean.sdCI95(1), r.mean.sdCI95(2))};
msgbox(msg,'Spread across molecules');

function redrawPlot(hObject, eventdata, handles)
% Redraw the histogram (with bootstrap error bars if available) and, if a Gaussian fit exists, the fit.
% bar() clears the axes each time, so the Gaussian curves and state labels are
% redrawn here too, after every Load / Replot / Fit / Save / checkbox toggle.
plotHistogram(handles.histogramPlot, hObject, eventdata, handles);
if ~isempty(handles.bobaResults) && get(handles.showRanges,'Value')==1
    drawBobaRanges(hObject, eventdata, handles);   % shading first, so the fit curve draws on top of it
    handles=guidata(hObject);
end
if handles.fitFlag==1
    fitGaussPlots(handles.histogramPlot, handles.R2_value, handles.parameter, hObject, eventdata, handles);
end
if ~isempty(handles.bobaResults) && get(handles.showLabels,'Value')==1
    drawBobaLabels(hObject, eventdata, handles);
end

function showErrorBars_Callback(hObject, eventdata, handles)
if isempty(handles.histogramData), return; end
redrawPlot(hObject, eventdata, handles);

function resetLabels_Callback(hObject, eventdata, handles)
% Clears manual label nudges only; does not touch the histogram, the fit, or bobaResults.
if isempty(handles.bobaResults), return; end
handles.bobaLabelOffsets=zeros(size(handles.bobaLabelOffsets));
guidata(hObject, handles);
redrawPlot(hObject, eventdata, handles);

function rangeSigmaPreset_Callback(hObject, eventdata, handles)
% Picking a preset writes the matching multiplier into the "Range (+/-sigma)"
% box and redraws. "Custom" (the last entry) does nothing here -- it's only
% ever selected automatically, when rangeSigmaMultiplier_Callback below sees
% the box was hand-edited to a value that doesn't match a preset.
presets = [1 1.5 2 3];   % sigma multipliers, in the same order as the popup's String
idx = get(hObject,'Value');
if idx <= numel(presets)
    set(handles.rangeSigmaMultiplier,'String', num2str(presets(idx)));
end
guidata(hObject, handles);
if isempty(handles.histogramData), return; end
redrawPlot(hObject, eventdata, handles);

function rangeSigmaMultiplier_Callback(hObject, eventdata, handles)
% Typing directly into the "Range (+/-sigma)" box still works; this just keeps
% the preset popup in sync -- "Custom" if the typed value doesn't match one of
% the presets, otherwise the matching preset (so switching back and forth stays
% consistent either way you edit it).
presets = [1 1.5 2 3];
val = str2double(get(hObject,'String'));
match = find(abs(presets - val) < 1e-9, 1);
if isempty(match)
    set(handles.rangeSigmaPreset,'Value', numel(get(handles.rangeSigmaPreset,'String')));   % "Custom"
else
    set(handles.rangeSigmaPreset,'Value', match);
end
guidata(hObject, handles);
if isempty(handles.histogramData), return; end
redrawPlot(hObject, eventdata, handles);

function keep = statesAboveMinPop(pop, handles)
% True for each state whose FULL-DATA fit population is at or above the
% "Min population (%)" box. pop is r.full.pop (one value per state, from the
% full-data fit) -- the decision of whether a state counts as real is made
% once, from the best available estimate, not re-evaluated per bootstrap
% replicate. Every caller that reports/draws/exports states uses this same
% mask, so "hidden" means hidden everywhere consistently.
minPop = str2double(get(handles.minPopPercent,'String')) / 100;
if isnan(minPop) || minPop < 0
    minPop = 0;   % empty/invalid box = no filtering
end
keep = pop(:)' >= minPop;

function drawBobaLabels(hObject, eventdata, handles, ax, persistHandles)
% Draws (or redraws) the "State k: pop% +/- sd%" labels in a row along the top of the
% plot, each centered over its state; overlapping labels stack downward.
% This only places text on top of the existing plot; it reads handles.bobaResults
% (already computed by boba_Callback) and never recomputes or changes it.
%
% ax             axes to draw into (default: the live GUI plot, handles.histogramPlot).
%                 Pass a different axes (e.g. the Save button's off-screen figure) to
%                 put the same labels on a saved image.
% persistHandles  true (default) stores the new label handles/positions back into
%                 handles via guidata, and makes labels draggable, as in the live GUI.
%                 Pass false for a one-off target (like a Save figure that gets closed
%                 right after) so this call can't overwrite the live GUI's label state.
if nargin<4 || isempty(ax),             ax = handles.histogramPlot; end
if nargin<5 || isempty(persistHandles), persistHandles = true;       end
r=handles.bobaResults;
K=numel(r.full.mu);
keep=statesAboveMinPop(r.full.pop, handles);   % states below "Min population (%)" get no label
if ~isequal(size(handles.bobaLabelOffsets),[K 2])
    handles.bobaLabelOffsets=zeros(K,2);   % number of states changed since the last fit; drop old offsets
end
% Default layout: every label sits in a row at the TOP of the plot, centered over its
% state's center. Where two labels would overlap horizontally, the later (further right)
% one steps down one label-height at a time until it clears. Manual drag offsets are
% added on top of this default, as before.
yl=ylim(ax);  xl=xlim(ax);
yTop=yl(2)-0.02*diff(yl);                  % top edge of the first row of labels
defaultPos=zeros(K,2);
labels=gobjects(K,1);
hold(ax,'on');
for k=1:K
    if ~keep(k), continue; end
    labelStr=sprintf('State %d\n%.1f%% \\pm %.1f%%', k, 100*r.summary.popMean(k), 100*r.summary.popStd(k));
    if persistHandles
        buttonDownFcn = @(src,evt) startDragLabel(src,evt,ancestor(ax,'figure'));
    else
        buttonDownFcn = '';   % a one-off (e.g. saved) copy is never dragged
    end
    % created at the top first, only so its on-screen size (Extent) can be measured
    labels(k)=text(ax, r.full.mu(k), yTop, labelStr, 'HorizontalAlignment','center', ...
        'VerticalAlignment','top', 'FontSize',9, 'FontWeight','bold', 'Color',[0 0 0.6], ...
        'Margin',1, 'ButtonDownFcn', buttonDownFcn, 'UserData', k);
end
placed=zeros(0,4);                         % rectangles already taken: [x1 x2 y1 y2]
for k=find(keep)                           % states are sorted by center, so this runs left to right
    ext=get(labels(k),'Extent');           % [x y width height] in axes (data) units
    w=ext(3);  h=ext(4);
    cx=min(max(r.full.mu(k), xl(1)+w/2), xl(2)-w/2);   % keep the label inside the axes
    yt=yTop;
    while yt-h > yl(1)
        clash=any(placed(:,1) < cx+w/2 & placed(:,2) > cx-w/2 & placed(:,3) < yt & placed(:,4) > yt-h);
        if ~clash, break; end
        yt=yt-1.05*h;                      % step down one label-height
    end
    placed(end+1,:)=[cx-w/2, cx+w/2, yt-h, yt]; %#ok<AGROW>
    defaultPos(k,:)=[cx, yt];
    set(labels(k),'Position',[defaultPos(k,:)+handles.bobaLabelOffsets(k,:), 0]);
end
hold(ax,'off');
if persistHandles
    handles.bobaLabels=labels;
    handles.bobaLabelDefaultPos=defaultPos;
    guidata(hObject, handles);
end

function drawBobaRanges(hObject, eventdata, handles, ax, persistHandles)
% Cosmetic only: shades each state's mean +/- N SD range, limited to the area under
% that state's fitted Gaussian curve (a low-opacity patch), with a dashed edge line at
% each end of the range that rises only as high as the curve there. Reads handles.bobaResults
% (already computed by boba_Callback) and does not recompute or change it.
%
% "SD" here is the fitted state's own width (full.sigma, the native Gaussian
% parameter -- no FWHM conversion involved), i.e. how broad that Gaussian is --
% not the bootstrap uncertainty on its center (summary.muStd), which is a
% different, usually much smaller, number.
% Swap r.full.sigma for r.summary.muStd below if you want the latter instead.
%
% The shaded band extends handles.rangeSigmaMultiplier standard deviations on
% each side of the center (default 1, i.e. the classic +/-1 sigma band); change
% the number in that box to widen/narrow it, e.g. 1.96 for a ~95% interval.
%
% ax, persistHandles -- same meaning as in drawBobaLabels, above.
if nargin<4 || isempty(ax),             ax = handles.histogramPlot; end
if nargin<5 || isempty(persistHandles), persistHandles = true;       end
r=handles.bobaResults;
K=numel(r.full.mu);
yl=ylim(ax);   % lock to the current y-range so patches/lines don't autoscale it

rangeMult = str2double(get(handles.rangeSigmaMultiplier,'String'));
if isnan(rangeMult) || rangeMult<=0
    rangeMult = 1;   % fall back to the classic +/-1 sigma band if the box is empty/invalid
end
keep = statesAboveMinPop(r.full.pop, handles);   % states below "Min population (%)" get no shaded range

% One color per state, cycled if there are more states than colors.
% These are ordinary MATLAB [R G B] triplets (0-1 each). MATLAB patch/line
% colors are always given as R,G,B in that order -- there is no BGR mode to
% switch on. If you want e.g. state 1 to render the way [0.85 0.33 0.10]
% (an orange-red) would look in BGR order, use its reversed triplet
% [0.10 0.33 0.85] (a blue) instead; that is the only way to get "BGR" colors
% out of MATLAB, by swapping the numbers yourself. Edit palette to taste.
palette=[0.85 0.33 0.10;   % state 1
         0.00 0.45 0.74;   % state 2
         0.47 0.67 0.19;   % state 3
         0.93 0.69 0.13;   % state 4
         0.49 0.18 0.56];  % state 5

patches=gobjects(K,1);
borders=gobjects(K,2);
hold(ax,'on');
for k=1:K
    if ~keep(k), continue; end
    col=palette(mod(k-1,size(palette,1))+1,:);
    lo=r.full.mu(k)-rangeMult*r.full.sigma(k);
    hi=r.full.mu(k)+rangeMult*r.full.sigma(k);
    % Shade only the area UNDER this state's fitted Gaussian curve, between lo and hi
    % (not a full-height column): the polygon follows the curve along its top and runs
    % back along the baseline. Uses the same BOBA full-data fit (amp, mu, sigma) as the
    % range itself, and is clipped to the visible x- and y-range.
    xl = xlim(ax);
    xs = linspace(max(lo,xl(1)), min(hi,xl(2)), 200);
    if numel(xs)>1 && xs(end)>xs(1)
        gk = @(x) r.full.amp(k) .* gaussian(x, r.full.mu(k), r.full.sigma(k));
        ys = min(gk(xs), yl(2));
        patches(k)=patch(ax,[xs fliplr(xs)],[ys yl(1)*ones(size(xs))],col, ...
            'FaceAlpha',0.3,'EdgeColor','none','HitTest','off');
        % dashed edge lines rise only as high as the curve at each end of the range
        borders(k,1)=line(ax,[lo lo],[yl(1) min(gk(lo),yl(2))],'Color',col,'LineWidth',1.2,'LineStyle','--','HitTest','off');
        borders(k,2)=line(ax,[hi hi],[yl(1) min(gk(hi),yl(2))],'Color',col,'LineWidth',1.2,'LineStyle','--','HitTest','off');
    end
end
hold(ax,'off');
ylim(ax,yl);
handles.bobaRanges=patches;
handles.bobaRangeBorders=borders;
guidata(hObject, handles);

function startDragLabel(src, evt, figHandle)
% Begins a drag: remember which label, and route mouse movement/release to it.
handles=guidata(figHandle);
handles.dragLabelIndex=get(src,'UserData');
guidata(figHandle, handles);
set(figHandle,'WindowButtonMotionFcn',@(s,e) dragLabelMotion(s,e,figHandle));
set(figHandle,'WindowButtonUpFcn',@(s,e) stopDragLabel(s,e,figHandle));

function dragLabelMotion(src, evt, figHandle)
% Moves only the dragged text object's on-screen position; no data changes.
handles=guidata(figHandle);
if isempty(handles.dragLabelIndex), return; end
cp=get(handles.histogramPlot,'CurrentPoint');
set(handles.bobaLabels(handles.dragLabelIndex),'Position',[cp(1,1) cp(1,2) 0]);

function stopDragLabel(src, evt, figHandle)
% Ends the drag and records the manual offset so it survives the next redraw.
handles=guidata(figHandle);
k=handles.dragLabelIndex;
if ~isempty(k)
    p=get(handles.bobaLabels(k),'Position');
    handles.bobaLabelOffsets(k,:)=p(1:2)-handles.bobaLabelDefaultPos(k,:);
end
handles.dragLabelIndex=[];
set(figHandle,'WindowButtonMotionFcn','');
set(figHandle,'WindowButtonUpFcn','');
guidata(figHandle, handles);

function writeBobaSummary(fileName, r, handles)
% handles is only used to read the "Min population (%)" box (statesAboveMinPop);
% nothing here recomputes or changes the fit itself.
fid=fopen(fileName,'w');
if fid<0, warning('Could not write %s', fileName); return; end
s=r.summary;
keep=statesAboveMinPop(r.full.pop, handles);
fprintf(fid,'molecules,%d,replicates_converged,%d,replicates_total,%d,R2_full_fit,%.4f\n', r.nMolecules, r.boot.nValid, r.boot.nBoot, r.full.R2);
if any(~keep)
    fprintf(fid,'# %d state(s) below %s%% population omitted below (still counted in every other state''s population)\n', ...
        sum(~keep), get(handles.minPopPercent,'String'));
end
% The last four columns (center and sigma 95% intervals) were added after the original
% layout, so they're appended at the end to keep earlier columns in the same positions.
% Results computed before that addition have no muCI95/sigmaCI95 and write NaN there.
if isfield(s,'muCI95'),    muCI=s.muCI95;       else, muCI=nan(2,numel(s.popMean));    end
if isfield(s,'sigmaCI95'), sgCI=s.sigmaCI95;    else, sgCI=nan(2,numel(s.popMean));    end
fprintf(fid,'state,center_mean,center_sd,sigma_mean,sigma_sd,population_mean,population_sd,population_ci95_low,population_ci95_high,center_ci95_low,center_ci95_high,sigma_ci95_low,sigma_ci95_high\n');
for k=1:numel(s.popMean)
    if ~keep(k), continue; end
    fprintf(fid,'%d,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n', k, s.muMean(k), s.muStd(k), s.sigmaMean(k), s.sigmaStd(k), ...
        s.popMean(k), s.popStd(k), s.popCI95(1,k), s.popCI95(2,k), muCI(1,k), muCI(2,k), sgCI(1,k), sgCI(2,k));
end
% per-bin histogram values and their bootstrap error bars
fprintf(fid,'\nbin_center,mean_normalized_count,bootstrap_se,ci95_low,ci95_high\n');
for i=1:numel(r.hist.centers)
    fprintf(fid,'%.5f,%.6f,%.6f,%.6f,%.6f\n', r.hist.centers(i), r.hist.mean(i), r.hist.se(i), r.hist.ci95(1,i), r.hist.ci95(2,i));
end
fclose(fid);

function writeLabelOffsets(fileName, offsets)
% Manual label-position nudges only (display, not data); written so a saved
% figure's labels can be reproduced. Does nothing if no label was dragged.
if isempty(offsets) || ~any(offsets(:))
    return
end
fid=fopen(fileName,'w');
if fid<0, warning('Could not write %s', fileName); return; end
fprintf(fid,'state,dx,dy\n');
for k=1:size(offsets,1)
    fprintf(fid,'%d,%.5f,%.5f\n', k, offsets(k,1), offsets(k,2));
end
fclose(fid);

function writeSpreadSummary(fileName, r)
fid=fopen(fileName,'w');
if fid<0, warning('Could not write %s', fileName); return; end
fprintf(fid,'per_molecule_value,n_molecules,mean,sd_across_molecules,se_of_sd,sd_ci95_low,sd_ci95_high\n');
fprintf(fid,'median_all_frames,%d,%.5f,%.5f,%.5f,%.5f,%.5f\n', r.median.n, r.median.mean, r.median.sd, r.median.sdSE, r.median.sdCI95(1), r.median.sdCI95(2));
fprintf(fid,'mean_in_range,%d,%.5f,%.5f,%.5f,%.5f,%.5f\n', r.mean.n, r.mean.mean, r.mean.sd, r.mean.sdSE, r.mean.sdCI95(1), r.mean.sdCI95(2));
fprintf(fid,'\nfile,n_frames,median_FRET,mean_FRET_in_range,fraction_in_range\n');
for i=1:numel(r.names)
    fprintf(fid,'%s,%d,%.5f,%.5f,%.5f\n', r.names{i}, r.nFrames(i), r.medianE(i), r.meanE(i), r.fracInRange(i));
end
fclose(fid);



%% ========================  helper functions   ====================================

function [histogramData, perMolHist]=CalculateHistogram(hObject, eventdata, handles)
BinSize =str2double(get(handles.binSize,'String'));
edges = -0.2+BinSize/2:BinSize:1.2+BinSize/2;
centers=edges(1)+BinSize/2:BinSize:edges(end)-BinSize/2;
pathName=handles.pathName;
filelist=handles.filelist;
T_histCount=zeros(size(centers));
perMolHist=zeros(length(filelist),length(centers));
for i=1:length(filelist)
data=importdata(strcat(pathName,filelist{i}));
FRET=data(:,2)./(data(:,2)+data(:,1));
histCount = histcounts(FRET,edges);
histCount=histCount./size(data,1);
perMolHist(i,:)=histCount;
T_histCount=T_histCount+histCount;
end
norm_T_histCount=T_histCount./length(filelist);
histogramData=[centers; norm_T_histCount];
%=======================================================================================================================================================
function plotHistogram(plotHandle,hObject, eventdata, handles)
centers=handles.histogramData(1,:);
norm_T_histCount=handles.histogramData(2,:);
axes(plotHandle)
bar(centers, norm_T_histCount, 'FaceColor',handles.faceColor, 'EdgeColor',handles.edgeColor, 'LineWidth',str2double(get(handles.edgeWidth,'String')));
set(plotHandle, 'XLim', [-0.2 1.2], 'XTick', -0.2:0.2:1.2,'XTickLabel', -0.2:0.2:1.2);
plotHandle.FontSize =str2double(get(handles.tickSize,'String'));
ylabel(plotHandle,'Normalized Count', 'FontSize',str2double(get(handles.axisSize,'String')) );
xlabel(plotHandle,'FRET Efficiency', 'FontSize',str2double(get(handles.axisSize,'String')) );
title(plotHandle,get(handles.plotTitle,'String'), 'FontSize', str2double(get(handles.titleSize,'String')));
% Bootstrap error bars (+/- 1 standard error of each bin), once the bootstrap has been run.
% The shaded band and the whisker caps toggle independently (showShading / showErrorBars).
if isfield(handles,'bobaResults') && ~isempty(handles.bobaResults)
    se=handles.bobaResults.hist.se;
    lo=norm_T_histCount-se;
    hi=norm_T_histCount+se;
    hold(plotHandle,'on');
    if get(handles.showShading,'Value')==1
        % One closed polygon: forward along the upper bound, backward along
        % the lower bound (fliplr on both halves). Drawn before the whiskers
        % so the caps land on top of it.
        % NOTE: the bootstrap SE band is often only a few pixels tall, so a
        % low-alpha fill with no outline is effectively invisible (it has no
        % minimum render width the way a line does). A visible edge color
        % traces the upper/lower bounds as actual lines regardless of how
        % thin the fill is, and a higher FaceAlpha keeps the fill itself
        % visible even stacked under the state-range shading.
        fill(plotHandle, [centers fliplr(centers)], [hi fliplr(lo)], [0.25 0.25 0.25], ...
            'FaceAlpha',0.4, 'EdgeColor',[0.25 0.25 0.25], 'LineWidth',1, 'HitTest','off');
    end
    if get(handles.showErrorBars,'Value')==1
        errorbar(plotHandle, centers, norm_T_histCount, se, se, 'LineStyle','none', 'Color',[0 0 0], 'LineWidth',1.2, 'CapSize',4);
    end
    hold(plotHandle,'off');
end
% A normalized count can never be negative, so the y-axis never should be
% either -- MATLAB's default autoscale otherwise pads a bit of blank space
% below zero (visible as a gap under the bars). Keep whatever top MATLAB
% chose (it already accounts for the error bar caps above); just clamp the
% bottom to exactly 0.
yl=ylim(plotHandle);
ylim(plotHandle,[0, yl(2)]);
guidata(hObject, handles);
%=======================================================================================================================================================
function theError = fitgauss(lambda, t, y)
% Fitting function for multiple overlapping Gaussians, with statements
% added (lines 18 and 19) to slow the progress and plot each step along the
% way, for educational purposes.
% Author: T. C. O'Haver, 2006

global c NumTrials TrialError
try
	
	A = zeros(length(t), round(length(lambda) / 2));
	for j = 1 : length(lambda) / 2
		A(:,j) = gaussian(t, lambda(2 * j - 1), lambda(2 * j))';
	end
	
	c = A \ y';
	z = A * c;
	theError = norm(z - y');
	
	% Penalty so that heights don't become negative.
	if sum(c < 0) > 0
		theError = theError + 1000000;
	end
	
	NumTrials = NumTrials + 1;
	TrialError(NumTrials) = theError;
catch ME
	% Some error happened if you get here.
	callStackString = GetCallStack(ME);
	errorMessage = sprintf('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', ...
		mfilename, callStackString, ME.message);
	WarnUser(errorMessage);
end
%=======================================================================================================================================================
function callStackString = GetCallStack(errorObject)
try
	theStack = errorObject.stack;
	callStackString = '';
	stackLength = length(theStack);
	% Get the date of the main, top level function:
	% 	d = dir(theStack(1).file);
	% 	fileDateTime = d.date(1:end-3);
	if stackLength <= 3
		% Some problem in the OpeningFcn
		% Only the first item is useful, so just alert on that.
		[folder, baseFileName, ext] = fileparts(theStack(1).file);
		baseFileName = sprintf('%s%s', baseFileName, ext);	% Tack on extension.
		callStackString = sprintf('%s in file %s, in the function %s, at line %d\n', callStackString, baseFileName, theStack(1).name, theStack(1).line);
	else
		% Got past the OpeningFcn and had a problem in some other function.
		for k = 1 : length(theStack)-3
			[folder, baseFileName, ext] = fileparts(theStack(k).file);
			baseFileName = sprintf('%s%s', baseFileName, ext);	% Tack on extension.
			callStackString = sprintf('%s in file %s, in the function %s, at line %d\n', callStackString, baseFileName, theStack(k).name, theStack(k).line);
		end
	end
catch ME
	errorMessage = sprintf('Error in program %s.\nTraceback (most recent at top):\nError Message:\n%s', ...
		mfilename, ME.message);
	WarnUser(errorMessage);
end
%=======================================================================================================================================================
% Pops up a warning message, and prints the error to the command window.
function WarnUser(warningMessage)
if nargin == 0
	return; % Bail out if they called it without any arguments.
end
try
	fprintf('%s\n', warningMessage);
	uiwait(warndlg(warningMessage));
	% Write the warning message to the log file
	folder = 'C:\Users\Public\Documents\MATLAB Settings';
	if ~exist(folder, 'dir')
		mkdir(folder);
	end
	fullFileName = fullfile(folder, 'Error Log.txt');
	fid = fopen(fullFileName, 'at');
	if fid >= 0
		fprintf(fid, '\nThe error below occurred on %s.\n%s\n', datestr(now), warningMessage);
		fprintf(fid, '-------------------------------------------------------------------------------\n');
		fclose(fid);
	end
catch ME
	message = sprintf('Error in WarnUser():\n%s', ME.message);
	fprintf('%s\n', message);
	uiwait(warndlg(message));
end
%=======================================================================================================================================================
function g = gaussian(x, peakPosition, sigma)
%  gaussian(x,pos,sigma) = gaussian peak centered on pos, with standard
%  deviation sigma. x may be scalar, vector, or matrix; pos and sigma both scalar.
%  Parameterized directly in sigma -- there is no FWHM conversion anywhere in
%  this function or in bobaFRET.m's matching gaussianSigma(). (If you ever need
%  FWHM for an external report, FWHM = 2*sqrt(2*ln(2))*sigma, but nothing in
%  this program computes or relies on that quantity.)
% Examples: gaussian([0 1 2],1,1) gives a peak centered at x=1 with sigma=1.
g = exp(-((x - peakPosition) .^ 2) ./ (2 .* sigma .^ 2));
%=======================================================================================================================================================
function yhat = PlotComponentCurves(x, y, t, c, parameter)
try
	fontSize = 20;
	% Get the means and widths.
	means = parameter(1 : 2 : end);
	widths = parameter(2 : 2 : end);
	% Now plot results.
	hFig2 = figure;
	hFig2.Name = 'Fitted Component Curves';
	% 	plot(x, y, '--', 'LineWidth', 2)
	hold on;
	yhat = zeros(1, length(t));
	numGaussians = length(c);
	legendStrings = cell(numGaussians + 2, 1);
	for k = 1 : numGaussians
		% Get each component curve.
		thisEstimatedCurve = c(k) .* gaussian(t, means(k), widths(k));
		% Plot component curves.
		plot(x, thisEstimatedCurve, '-', 'LineWidth', 2);
		hold on;
		% Overall curve estimate is the sum of the component curves.
		yhat = yhat + thisEstimatedCurve;
		legendStrings{k} = sprintf('Estimated Gaussian %d', k);
	end
	% Plot original summation curve, that is the actual curve.
	plot(x, y, 'r-', 'LineWidth', 1)
	% Plot estimated summation curve, that is the estimate of the curve.
	plot(x, yhat, 'k--', 'LineWidth', 2)
	grid on;
	xlabel('X', 'FontSize', fontSize)
	ylabel('Y', 'FontSize', fontSize)
	caption = sprintf('Estimation of %d Gaussian Curves that will fit data.', numGaussians);
	title(caption, 'FontSize', fontSize, 'Interpreter', 'none');
	grid on
	legendStrings{numGaussians+1} = sprintf('Actual original signal');
	legendStrings{numGaussians+2} = sprintf('Sum of all %d Gaussians', numGaussians);
	legend(legendStrings);
	xlim(sort([x(1) x(end)]));
	hFig2.WindowState = 'maximized';
	drawnow;
	
catch ME
	% Some error happened if you get here.
	callStackString = GetCallStack(ME);
	errorMessage = sprintf('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', ...
		mfilename, callStackString, ME.message);
	WarnUser(errorMessage);
end
% --- Executes on button press in fitColor.
function fitColor_Callback(hObject, eventdata, handles)
handles.fitColor = uisetcolor([0.6 0.8 1]);
guidata(hObject, handles);












%% Functions that are not actively needed....
% --- Executes on button press in faceColor.
function faceColor_Callback(hObject, eventdata, handles)
handles.faceColor = uisetcolor([0.6 0.8 1]);
guidata(hObject, handles);
% --- Executes on button press in edgeColor.
function edgeColor_Callback(hObject, eventdata, handles)
handles.edgeColor = uisetcolor([0.6 0.8 1]);
guidata(hObject, handles);


function titleSize_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function titleSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to titleSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function axisSize_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function axisSize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tickSize_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function tickSize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fitThickness_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function fitThickness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fitThickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function textSize_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function textSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function textPos_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function textPos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function savingfileName_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function savingfileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savingfileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function in_Centers_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function in_Centers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in_Centers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function in_sigma_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function in_sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edgeWidth_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function edgeWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edgeWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plotTitle_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function plotTitle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotTitle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function binSize_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function binSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to binSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
