function varargout = convertor(varargin)
% CONVERTOR MATLAB code for convertor.fig
%      CONVERTOR, by itself, creates a new CONVERTOR or raises the existing
%      singleton*.
%
%      H = CONVERTOR returns the handle to a new CONVERTOR or the handle to
%      the existing singleton*.
%
%      CONVERTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONVERTOR.M with the given input arguments.
%
%      CONVERTOR('Property','Value',...) creates a new CONVERTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before convertor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to convertor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help convertor

% Last Modified by GUIDE v2.5 19-Aug-2019 15:34:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @convertor_OpeningFcn, ...
                   'gui_OutputFcn',  @convertor_OutputFcn, ...
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


% --- Executes just before convertor is made visible.
function convertor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to convertor (see VARARGIN)

% Choose default command line output for convertor
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes convertor wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = convertor_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in edf2asc.
function edf2asc_Callback(hObject, eventdata, handles)
% hObject    handle to edf2asc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edf2asc


% --- Executes on button press in asc2mat.
function asc2mat_Callback(hObject, eventdata, handles)
% hObject    handle to asc2mat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of asc2mat



function filePath1_Callback(hObject, eventdata, handles)
% hObject    handle to filePath1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.edf2ascLog,'String','');
set(handles.asc2matLog,'String','');
% Hints: get(hObject,'String') returns contents of filePath1 as text
%        str2double(get(hObject,'String')) returns contents of filePath1 as a double


% --- Executes during object creation, after setting all properties.
function filePath1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filePath1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filePath2 = get(handles.filePath2,'String');

if isempty(filePath2)
    filePath1 = get(handles.filePath1,'String');
    if isempty(filePath1)
        set(handles.pathLog2,'String','Please input the data path!');
        set(handles.pathLog2,'ForegroundColor',[1 0 0]);
        return
    else
        set(handles.filePath2,'String',filePath1);
        set(handles.pathLog2,'String','Is this way right?');
        set(handles.pathLog2,'ForegroundColor',[1 0 0]);
        return
    end
else
    functionPath = pwd;
    try
        cd(filePath2);
        set(handles.pathLog2,'String','Start reading files.')
        set(handles.pathLog2,'ForegroundColor',[0.2 0.9 0]);
    catch
        set(handles.pathLog2,'String','Path error');
        set(handles.pathLog2,'ForegroundColor',[1 0 0]);
        return
    end
end

dataFile = dir(fullfile(filePath2,'Converted_*.mat'));
fileName = cell(length(dataFile),1);
for i = 1:length(dataFile)
    fileName{i} = dataFile(i).name;
end
set(handles.fileList,'String',fileName)
cd(functionPath);
    

% --- Executes on selection change in fileList.
function fileList_Callback(hObject, eventdata, handles)
% hObject    handle to fileList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fileList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fileList
filePath2 = get(handles.filePath2,'String');
if isempty(filePath2)
    filePath1 = get(handles.filePath1,'String');
    if isempty(filePath1)
        set(handles.pathLog2,'String','Please input the data path!');
        set(handles.pathLog2,'ForegroundColor',[1 0 0]);
        return
    else
        set(handles.filePath2,'String',filePath1);
        set(handles.pathLog2,'String','Will work on displayed pathway.');
        set(handles.pathLog2,'ForegroundColor',[1 0 0]);
    end
else
    set(handles.pathLog2,'String',' ');
    set(handles.pathLog2,'ForegroundColor',[0.2 0.9 0]);
end

fileList = get(handles.fileList,'String');
fileNum = get(handles.fileList,'Value');
fileName = fullfile(filePath2,fileList{fileNum});

clear data
data = load(fileName);

if ~isempty(data.eyeData)
    figure(1)
    cla;
    
    plot(data.eyeData(:,2),'r');
    hold on
    plot(data.eyeData(:,3),'g');
    plot(data.eyeData(:,4),'b');
    
    if ~isempty(data.movingStart)
        msIndex = ones(length(data.movingStart(:,1)),1);
        for i = 1:length(data.movingStart(:,1))
            msIndex(i) = find(data.movingStart(i,1) <= data.eyeData(:,1),1);
        end
        plot([msIndex msIndex],[0 8000],'--m');
    end
    
    if ~isempty(data.startChoice)
        scIndex = ones(length(data.startChoice(:,1)),1);
        for i = 1:length(data.startChoice(:,1))
            scIndex(i) = find(data.startChoice(i,1) <= data.eyeData(:,1),1);
        end
        plot([scIndex scIndex],[0 8000],'--k');
    end
else
    set(handles.pathLog2,'String','No eye data in this file.');
    set(handles.pathLog2,'ForegroundColor',[1 0 0]);
end
    






% --- Executes during object creation, after setting all properties.
function fileList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filePath1 = get(handles.filePath1,'String');
edf2asc = get(handles.edf2asc,'Value');

if isempty(filePath1)
    set(handles.pathLog1,'String','You need to input file path before start!')
    set(handles.pathLog1,'ForegroundColor',[1 0 0]);
    return
else
    functionPath = pwd;
    exePath = fullfile(functionPath,'edf2asc.exe');
    exeCmd = [exePath ' -ntime_chcek'];
    
    try
        cd(filePath1);
        set(handles.pathLog1,'String','Convert in processing...')
        set(handles.pathLog1,'ForegroundColor',[0.2 0.9 0]);
    catch
        set(handles.pathLog1,'String','Path error');
        set(handles.pathLog1,'ForegroundColor',[1 0 0]);
        return
    end
end

if edf2asc
    edfFiles = dir('*.edf');
    edfFileNum = length(edfFiles);
    
    state = 0;
    set(handles.edf2ascLog,'String','Please Wait... 0%');
    set(handles.edf2ascLog,'ForegroundColor',[0.9 0.7 0]);
    WaitSecs(0.2);
    
    for i = 1:edfFileNum
        edfNamei = edfFiles(i).name;
        ascNamei = strrep(edfNamei,'.edf','.asc');
        edfPathi = fullfile(filePath1,edfNamei);
        if exist(ascNamei,'file')
            delete(ascNamei)
        end
        cmd = [exeCmd 32 edfPathi];
        system(cmd);
        
        if i < edfFileNum/4
            curState = 0;
        elseif i < edfFileNum/2 && i > edfFileNum/4
            curState = 1;
        elseif i < edfFileNum*3/4 && i > edfFileNum/2
            curState = 2;
        elseif i > edfFileNum*3/4 && i < edfFileNum
            curState = 3;
        elseif i == edfFileNum
            curState = 4;
        end
        
        if curState > state
            if curState == 1
                set(handles.edf2ascLog,'String','Please Wait... 25%');
                set(handles.edf2ascLog,'ForegroundColor',[0.9 0.7 0]);
                WaitSecs(0.2);
            elseif  curState == 2
                set(handles.edf2ascLog,'String','Please Wait... 50%');
                set(handles.edf2ascLog,'ForegroundColor',[0.9 0.7 0]);
                WaitSecs(0.2);
            elseif curState == 3
                set(handles.edf2ascLog,'String','Please Wait... 75%');
                set(handles.edf2ascLog,'ForegroundColor',[0.9 0.7 0]);
                WaitSecs(0.2);
            elseif curState == 4
                set(handles.edf2ascLog,'String','Mission success.');
                set(handles.edf2ascLog,'ForegroundColor',[0.2 0.9 0]);
                WaitSecs(0.2);
            end
            state = curState;
        end
    end
end

asc2mat = get(handles.asc2mat,'Value');
if asc2mat
    ascFile = dir('*.asc');
    ascFileNum = length(ascFile);
    
    state = 0;
    set(handles.asc2matLog,'String','Please Wait... 0%');
    set(handles.asc2matLog,'ForegroundColor',[0.9 0.7 0]);
    WaitSecs(0.5);
    
    for i = 1:ascFileNum
        tic;
        ascNamei = ascFile(i).name;
        saveNamei = ['Converted_' strrep(ascNamei,'.asc','.mat')];
        
        if exist(saveNamei,'file')
            delete(saveNamei)
        end
        
        fid = fopen(ascNamei);
        fseek(fid,0,'eof');
        numline=ftell(fid);
        fclose(fid);
        rawData=importdata(ascNamei,' ',numline);
        
        msStr = rawData(contains(rawData,'Moving Start','IgnoreCase',true));
        movingStart = nan(length(msStr),2);
        for j = 1:length(msStr)
            movingStart(j,:) = cell2mat(cellfun(@str2num,regexp(msStr{j},'\d*\.?\d*','match'),'UniformOutput',false));
        end
        delIndex = diff(movingStart(:,2))==0;
        movingStart(delIndex,:) = [];
        
        scStr = rawData(contains(rawData,'Start choice','IgnoreCase',true));
        startChoice = nan(length(scStr),2);
        for j = 1:length(scStr)
            startChoice(j,:) = cell2mat(cellfun(@str2num,regexp(scStr{j},'\d*\.?\d*','match'),'UniformOutput',false));
        end
        delIndex = diff(startChoice(:,2))==0;
        startChoice(delIndex,:) = [];
        
        dmStr = rawData(contains(rawData,'Decision made','IgnoreCase',true));
        decisionMade = nan(length(dmStr),2);
        for j = 1:length(dmStr)
            decisionMade(j,:) = cell2mat(cellfun(@str2num,regexp(dmStr{j},'\d*\.?\d*','match'),'UniformOutput',false));
        end
        delIndex = diff(decisionMade(:,2))==0;
        decisionMade(delIndex,:) = [];
        
        tcStr = rawData(contains(rawData,'Trial complete','IgnoreCase',true));
        trialComplete = nan(length(tcStr),2);
        for j = 1:length(tcStr)
            trialComplete(j,:) = cell2mat(cellfun(@str2num,regexp(tcStr{j},'\d*\.?\d*','match'),'UniformOutput',false));
        end
        delIndex = diff(trialComplete(:,2))==0;
        trialComplete(delIndex,:) = [];
        
        trStr = rawData(contains(rawData,'Trial repeat','IgnoreCase',true));
        trialRepeat = nan(length(trStr),2);
        for j = 1:length(trStr)
            trialRepeat(j,:) = cell2mat(cellfun(@str2num,regexp(trStr{j},'\d*\.?\d*','match'),'UniformOutput',false));
        end
        trialRepeat;
        
        eyeDataInd = contains(rawData,'...');
        eyeData = rawData(eyeDataInd);
        eyeData = cell2mat(cellfun(@str2num,strrep(eyeData,'...',''),'UniformOutput',false));
        
        save(saveNamei,'eyeData','movingStart','startChoice','decisionMade','trialComplete','trialRepeat')
        
        if i < ascFileNum/4
            curState = 0;
        elseif i < ascFileNum/2 && i > ascFileNum/4
            curState = 1;
        elseif i < ascFileNum*3/4 && i > ascFileNum/2
            curState = 2;
        elseif i > ascFileNum*3/4 && i < ascFileNum
            curState = 3;
        elseif i == ascFileNum
            curState = 4;
        end
        
        if curState > state
            if curState == 1
                set(handles.asc2matLog,'String','Please Wait... 25%');
                set(handles.asc2matLog,'ForegroundColor',[0.9 0.7 0]);
                WaitSecs(0.5);
            elseif  curState == 2
                set(handles.asc2matLog,'String','Please Wait... 50%');
                set(handles.asc2matLog,'ForegroundColor',[0.9 0.7 0]);
                WaitSecs(0.5);
            elseif curState == 3
                set(handles.asc2matLog,'String','Please Wait... 75%');
                set(handles.asc2matLog,'ForegroundColor',[0.9 0.7 0]);
                WaitSecs(0.5);
            elseif curState == 4
                set(handles.asc2matLog,'String','Mission success.');
                set(handles.asc2matLog,'ForegroundColor',[0.2 0.9 0]);
                WaitSecs(0.5);
            end
            state = curState;
        end
    end
end
set(handles.pathLog1,'String','Mission complete.')
set(handles.pathLog1,'ForegroundColor',[0.2 0.9 0]);
set(handles.filePath2,'String',filePath1);
cd(functionPath);
toc



function filePath2_Callback(hObject, eventdata, handles)
% hObject    handle to filePath2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filePath2 as text
%        str2double(get(hObject,'String')) returns contents of filePath2 as a double


% --- Executes during object creation, after setting all properties.
function filePath2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filePath2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
