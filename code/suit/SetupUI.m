classdef SetupUI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                     matlab.ui.Figure
        PROCESSFLAGSPanel            matlab.ui.container.Panel
        CONNSTEPSPanel               matlab.ui.container.Panel
        DropDown                     matlab.ui.control.DropDown
        AddROICheckBoxF              matlab.ui.control.CheckBox
        ndlevelCheckBoxF             matlab.ui.control.CheckBox
        stlevelCheckBoxF             matlab.ui.control.CheckBox
        DenoisingCheckBoxF           matlab.ui.control.CheckBox
        PreprocessingCheckBoxF       matlab.ui.control.CheckBox
        SetupCheckBoxF               matlab.ui.control.CheckBox
        MRISetupCheckBoxF            matlab.ui.control.CheckBox
        OverWriteCheckBoxF           matlab.ui.control.CheckBox
        TextArea                     matlab.ui.control.TextArea
        FunctionalProgressPanel      matlab.ui.container.Panel
        PlotflapmapCheckBox          matlab.ui.control.CheckBox
        ReslicetoSUITCheckBox        matlab.ui.control.CheckBox
        analysisCheckBox             matlab.ui.control.CheckBox
        ndlevelCheckBox              matlab.ui.control.CheckBox
        SmoothCheckBox               matlab.ui.control.CheckBox
        contrastCheckBox             matlab.ui.control.CheckBox
        NegCheckBox                  matlab.ui.control.CheckBox
        PosCheckBox                  matlab.ui.control.CheckBox
        stlevelCheckBox              matlab.ui.control.CheckBox
        CorigistertoT1CheckBox       matlab.ui.control.CheckBox
        SliceTimingCheckBox          matlab.ui.control.CheckBox
        RealignmentCheckBox          matlab.ui.control.CheckBox
        RunCheckBoxFunc              matlab.ui.control.CheckBox
        AnatomicalProgressPanel      matlab.ui.container.Panel
        NormalizeCheckBox            matlab.ui.control.CheckBox
        SegmentandIsolateCheckBox    matlab.ui.control.CheckBox
        RunCheckBoxAna               matlab.ui.control.CheckBox
        DataprogressPanel            matlab.ui.container.Panel
        subidentifyEditField         matlab.ui.control.EditField
        subidentifyEditFieldLabel    matlab.ui.control.Label
        GetDataCheckBox              matlab.ui.control.CheckBox
        RawDataCheckBox              matlab.ui.control.CheckBox
        progressDropDown             matlab.ui.control.DropDown
        progressDropDownLabel        matlab.ui.control.Label
        behaveCheckBox               matlab.ui.control.CheckBox
        localCheckBox                matlab.ui.control.CheckBox
        RawTabGroup                  matlab.ui.container.TabGroup
        InitialTab                   matlab.ui.container.Tab
        condNameEditField            matlab.ui.control.EditField
        condNameEditFieldLabel       matlab.ui.control.Label
        ChckImgOpEditField           matlab.ui.control.EditField
        ChckImgOpEditField_2Label    matlab.ui.control.Label
        FolderNestEditField          matlab.ui.control.EditField
        FolderNestEditField_3Label   matlab.ui.control.Label
        subpathEditField             matlab.ui.control.EditField
        subpathEditField_3Label      matlab.ui.control.Label
        TextArea_3                   matlab.ui.control.TextArea
        BehaveFileinfoPanel          matlab.ui.container.Panel
        ExtentionEditFieldBEH        matlab.ui.control.EditField
        ExtentionEditField_2Label    matlab.ui.control.Label
        FolderNestEditFieldBEH       matlab.ui.control.EditField
        FolderNestEditField_4Label   matlab.ui.control.Label
        subpathEditFieldBEH          matlab.ui.control.EditField
        subpathEditField_4Label      matlab.ui.control.Label
        UITableRAW                   matlab.ui.control.Table
        FTPinfo                      matlab.ui.control.EditField
        FTPServerDropDown            matlab.ui.control.DropDown
        FTPServerDropDownLabel       matlab.ui.control.Label
        ConnTabGroup                 matlab.ui.container.TabGroup
        CONNInitialTab               matlab.ui.container.Tab
        PrepStepsEditField           matlab.ui.control.EditField
        PrepStepsEditFieldLabel      matlab.ui.control.Label
        ROIcorTable                  matlab.ui.control.Table
        OldprojectCheckBox           matlab.ui.control.CheckBox
        ResultFoldEditField          matlab.ui.control.EditField
        ResultFoldEditFieldLabel     matlab.ui.control.Label
        sliceorderEditField          matlab.ui.control.EditField
        sliceorderEditFieldLabel     matlab.ui.control.Label
        FiltBandHEditField           matlab.ui.control.NumericEditField
        FiltBandHEditFieldLabel      matlab.ui.control.Label
        AnaSourceEditField           matlab.ui.control.EditField
        AnaSourceEditFieldLabel      matlab.ui.control.Label
        AnaNameEditField             matlab.ui.control.EditField
        AnaNameEditFieldLabel        matlab.ui.control.Label
        ROIpathEditField             matlab.ui.control.EditField
        ROIpathEditFieldLabel        matlab.ui.control.Label
        DataPathEditField            matlab.ui.control.EditField
        DataPathEditFieldLabel       matlab.ui.control.Label
        NameEditField                matlab.ui.control.EditField
        NameEditField_2Label         matlab.ui.control.Label
        FiltBandLEditField           matlab.ui.control.NumericEditField
        FiltBandLEditFieldLabel      matlab.ui.control.Label
        FuncVresEditField            matlab.ui.control.NumericEditField
        FuncVresEditField_2Label     matlab.ui.control.Label
        StrucVresEditField           matlab.ui.control.NumericEditField
        StrucVresEditField_2Label    matlab.ui.control.Label
        fwhmEditField                matlab.ui.control.NumericEditField
        fwhmEditField_2Label         matlab.ui.control.Label
        TREditField                  matlab.ui.control.NumericEditField
        TREditField_2Label           matlab.ui.control.Label
        ProjectPathEditField         matlab.ui.control.EditField
        ProjectPathEditField_2Label  matlab.ui.control.Label
        CONNTab                      matlab.ui.control.Table
    end

    
    methods (Access = private)
        function InitialVar(app,varargin)
            f = []; o = [];
            if nargin > 1
                f = varargin{1};
            end
            if nargin > 2
                o = varargin{2};
            end
            if ~isempty(f)
                if f.conn
                    app.progressDropDown.Value = "CONN";
                else
                    app.progressDropDown.Value = "Raw";
                end
                app.localCheckBox.Value = f.local;
                app.behaveCheckBox.Value = f.beh;
                app.RawDataCheckBox.Value = f.raw;
                app.GetDataCheckBox.Value = f.GetData;
                % anatomical flags
                app.RunCheckBoxAna.Value = f.ana.ana;
                app.SegmentandIsolateCheckBox.Value = f.ana.SegIso;
                app.NormalizeCheckBox.Value = f.ana.Norm;
                
                % functional flags
                app.RunCheckBoxFunc.Value = f.func.func;
                app.RealignmentCheckBox.Value = f.func.relign;
                app.SliceTimingCheckBox.Value = f.func.sliTim;
                app.CorigistertoT1CheckBox.Value = f.func.corig;
                app.stlevelCheckBox.Value = f.func.fstL;
                app.PosCheckBox.Value = f.func.fstconP;
                app.NegCheckBox.Value = f.func.fstconN;
                app.contrastCheckBox.Value = f.func.fstcon;
                app.ReslicetoSUITCheckBox.Value = f.func.resliceSUIT;
                app.SmoothCheckBox.Value = f.func.smooth;
                app.ndlevelCheckBox.Value = f.gfunc.func;
                app.analysisCheckBox.Value = f.gfunc.sndL;
                app.PlotflapmapCheckBox.Value = f.gfunc.flapmapPlot;
                app.OldprojectCheckBox.Value = f.old_connPrj;
 
            else
                app.progressDropDown.Value = "Raw";
                app.localCheckBox.Value = false;
                app.behaveCheckBox.Value = false;
                app.RawDataCheckBox.Value = false;
                app.GetDataCheckBox.Value = false;
                % anatomical flags
                app.RunCheckBoxAna.Value = false;
                app.SegmentandIsolateCheckBox.Value = false;
                app.NormalizeCheckBox.Value = false;
                
                % functional flags
                app.RunCheckBoxFunc.Value = false;
                app.RealignmentCheckBox.Value = false;
                app.SliceTimingCheckBox.Value = false;
                app.CorigistertoT1CheckBox.Value = false;
                app.stlevelCheckBox.Value = false;
                app.PosCheckBox.Value = false;
                app.NegCheckBox.Value = false;
                app.NegCheckBox.Value = false;
                app.ReslicetoSUITCheckBox.Value = false;
                app.SmoothCheckBox.Value = false;
                app.ndlevelCheckBox.Value = false;
                app.analysisCheckBox.Value = false;
                app.PlotflapmapCheckBox.Value = false;
                app.OldprojectCheckBox.Value = false;
            end

            if ~isempty(o)
                %% SPM VARIABLES
                % ftpServer Data
                app.FTPServerDropDown.UserData = o.ftpServer;
                app.FTPinfo.Value = app.FTPServerDropDown.UserData.(app.FTPServerDropDown.Value);
                
                % subject identifier
                app.subidentifyEditField.Value = o.subidntfr;
                
                % subject file info
                app.subpathEditField.Value = o.subpath;
                app.FolderNestEditField.Value = cell2mat(strcat(o.FolderNest,filesep));
                app.ChckImgOpEditField.Value = o.opath;
                
                % subject name 
                rownum = max([length(o.sub),length(o.round),length(o.sess)]);
                o.sub = cat(1,o.sub,cell(rownum-length(o.sub),1));
                o.round = cat(1,o.round',cell(rownum-length(o.round),1));
                o.sess = cat(1,o.sess',cell(rownum-length(o.sess),1));
                app.UITableRAW.Data =  [o.sub,o.sess,o.round];
        
                % overwrite
                app.OverWriteCheckBoxF.Value = f.ow;

                % BEH file
                app.ExtentionEditFieldBEH.Value = o.BEH.Word.ext;
                app.FolderNestEditFieldBEH.Value = o.BEH.Word.folder;
                app.subpathEditFieldBEH.Value = o.BEH.Word.subpath;
                
                % condition name R
                app.condNameEditField.Value = o.R;

                %% CONN VARIABLES
                % CONN STEPS
                app.SetupCheckBoxF.Value = o.conn_steps.Setup;
                app.PreprocessingCheckBoxF.Value = o.conn_steps.Preprocessing;
                app.DenoisingCheckBoxF.Value = o.conn_steps.Denoising;
                app.stlevelCheckBoxF.Value = o.conn_steps.fst_Analysis;
                app.ndlevelCheckBoxF.Value = o.conn_steps.snd_Analysis;
                app.AddROICheckBoxF.Value = o.conn_steps.Add_Roi || o.conn_steps.Add_RoiResult;
                if o.conn_steps.Add_Roi
                    app.DropDown.Value = "Add roi";
                elseif o.conn_steps.Add_RoiResult
                    app.DropDown.Value = "Result";
                end
                

                % CONN var
                % project name and path
                app.ProjectPathEditField.Value = o.conn_prjPath;
                app.NameEditField.Value = o.conn_prjName;
                app.DataPathEditField.Value = o.subpath;

                % subject session round
                mL = max([length(o.conn_sub),length(o.conn_ses),length(o.conn_round)]);
                app.CONNTab.Data = [[o.conn_sub;cell(mL-length(o.conn_sub),1)], ...
                                     [o.conn_ses';cell(mL-length(o.conn_ses),1)], ...
                                     [o.conn_round';cell(mL-length(o.conn_round),1)]];
                % preprocessing
                app.TREditField.Value = o.conn_TR;
                app.fwhmEditField.Value = o.conn_fwhm;
                app.StrucVresEditField.Value = o.conn_strVres;
                app.FuncVresEditField.Value = o.conn_funcVres;
                app.sliceorderEditField.Value = num2str(o.conn_slice_order);
                % denoising
                app.FiltBandLEditField.Value = o.conn_filtBand(1);
                app.FiltBandHEditField.Value = o.conn_filtBand(2);
                % analysis
                app.AnaNameEditField.Value = o.conn_AnalysisName;
                app.AnaSourceEditField.Value = strjoin(o.conn_AnaSource,', ');
                % roi path
                app.ROIpathEditField.Value = o.conn_ROIpath;
                % Result folder
                app.ResultFoldEditField.Value = o.conn_ResultF;

                % roi corordinate
                Nroi = cell(length(cell2mat(o.ROIcor')),4);
                for i = 1:length(o.ROIname)
                    for j = 1:size(o.ROIcor{i},1)
                        Nroi{(i-1)*size(o.ROIcor{i},1)+j,1}  = [o.ROIname{i},'_',num2str(j)];
                        Nroi{(i-1)*size(o.ROIcor{i},1)+j,2} = o.ROIcor{i}(j,1);
                        Nroi{(i-1)*size(o.ROIcor{i},1)+j,3} = o.ROIcor{i}(j,2);
                        Nroi{(i-1)*size(o.ROIcor{i},1)+j,4} = o.ROIcor{i}(j,3);
                    end
                end
                app.ROIcorTable.Data = Nroi;
                app.PrepStepsEditField.Value = o.conn_prepsteps;

            else
              %% SPM VARIABLES
                % ftpServer Data
                app.FTPServerDropDown.UserData = struct('ip','','account','','password','','infolder','','outfolder','');
                app.FTPinfo.Value = '';
                
                % subject file info
                app.subpathEditField.Value = '';
                app.FolderNestEditField.Value = '';
                app.ChckImgOpEditField.Value = '';

                % subject identifier
                app.subidentifyEditField.Value = '';
                
                
                % subject name 
                app.UITableRAW.Data = cell(100,3);

                % overwrite
                app.OverWriteCheckBoxF.Value = false;

                % BEH file
                app.ExtentionEditFieldBEH.Value = '';
                app.FolderNestEditFieldBEH.Value = '';
                app.subpathEditFieldBEH.Value = '';

                % condition name R
                app.condNameEditField.Value = '';

                %% CONN VARIABLES
                % CONN STEPS
                app.SetupCheckBoxF.Value = false;
                app.PreprocessingCheckBoxF.Value = false;
                app.DenoisingCheckBoxF.Value = false;
                app.stlevelCheckBoxF.Value = false;
                app.ndlevelCheckBoxF.Value = false;
                app.AddROICheckBoxF.Value = false;

                % CONN var
                % subject session round
                app.CONNTab.Data = cell(100,3);

                % preprocessing
                app.TREditField.Value = 0;
                app.fwhmEditField.Value = 0;
                app.StrucVresEditField.Value = 0;
                app.FuncVresEditField.Value = 0;
                app.sliceorderEditField.Value = '';
                % Denoising
                app.FiltBandLEditField.Value = 0;
                app.FiltBandHEditField.Value = 0;
                % analysis
                app.AnaNameEditField.Value = '';
                app.AnaSourceEditField.Value = '';
                % roi path
                app.ROIpathEditField.Value = '';
                % result folder
                app.ResultFoldEditField.Value = '';

                % roi corordinate
                app.ROIcorTable.Data = cell(100,3);

                app.PrepStepsEditField.Value = '';
            end
        end
        function InitialVisible(app)
            % progress panel
            stlevelCheckBoxValueChanged(app);
            ndlevelCheckBoxValueChanged(app);
            SetupCheckBoxFValueChanged(app);
            PreprocessingCheckBoxFValueChanged(app);
            stlevelCheckBoxFValueChanged(app);
            DenoisingCheckBoxFValueChanged(app);
            localCheckBoxValueChanged(app);
            FTPServerDropDownOpening(app);
            FTPServerDropDownValueChanged(app);
            behaveCheckBoxValueChanged(app);
            progressDropDownValueChanged(app);
            RunCheckBoxFuncValueChanged(app);
            RunCheckBoxAnaValueChanged(app);
        end
        
        function [f,o] = SetupVar(app)
            value = app.progressDropDown.Value;
            if value == "CONN"
                f.conn = true;
            else
                f.conn = false;
            end
            f.local = app.localCheckBox.Value;
            f.beh = app.behaveCheckBox.Value;
            f.raw = app.RawDataCheckBox.Value;
            f.GetData = app.GetDataCheckBox.Value;
            f.ana.ana = app.RunCheckBoxAna.Value;
            f.ana.SegIso = app.SegmentandIsolateCheckBox.Value;
            f.ana.Norm = app.NormalizeCheckBox.Value;
            f.func.func = app.RunCheckBoxFunc.Value;
            f.func.relign = app.RealignmentCheckBox.Value;
            f.func.sliTim = app.SliceTimingCheckBox.Value;
            f.func.corig = app.CorigistertoT1CheckBox.Value;
            f.func.fstL = app.stlevelCheckBox.Value;
            f.ana.wholemask = false;
            f.func.fstconP = app.PosCheckBox.Value;
            f.func.fstconN = app.NegCheckBox.Value;
            f.func.fstcon = app.contrastCheckBox.Value;
            f.func.resliceSUIT = app.ReslicetoSUITCheckBox.Value;
            f.func.smooth = app.SmoothCheckBox.Value;
            f.gfunc.func = app.ndlevelCheckBox.Value;
            f.gfunc.sndL = app.analysisCheckBox.Value;
            f.gfunc.flapmapPlot = app.PlotflapmapCheckBox.Value;
            f.old_connPrj = app.OldprojectCheckBox.Value;

            o.ftpServer = app.FTPServerDropDown.UserData;
            o.subpath = app.subpathEditField.Value;
            value = app.FolderNestEditField.Value;
            value = split(value,filesep);
            value(cellfun(@isempty,value)) = [];
            o.FolderNest = value;
            o.opath = app.ChckImgOpEditField.Value;
            o.subidntfr = app.subidentifyEditField.Value;

            value = app.UITableRAW.Data;
            o.sub = value(:,1);
            o.sub(cellfun(@isempty,o.sub)) = [];
            o.sess = value(:,2)';         
            o.sess(cellfun(@isempty,o.sess)) = [];
            if isempty(o.sess), o.sess = {''}; end
            o.round = value(:,3)';
            o.round(cellfun(@isempty,o.round)) = [];
            if isempty(o.round), o.round = {''}; end
            f.ow = app.OverWriteCheckBoxF.Value;

            % BEH file
            o.BEH.Word.ext = app.ExtentionEditFieldBEH.Value;
            o.BEH.Word.folder = app.FolderNestEditFieldBEH.Value;
            o.BEH.Word.subpath = app.subpathEditFieldBEH.Value;

            % condition name R
            o.R = app.condNameEditField.Value;

            %% CONN VARIABLES
            % CONN STEPS
            o.conn_steps.Setup = app.SetupCheckBoxF.Value;
            o.conn_steps.Denoising = app.DenoisingCheckBoxF.Value;
            o.conn_steps.fst_Analysis = app.stlevelCheckBoxF.Value;
            o.conn_steps.snd_Analysis = app.ndlevelCheckBoxF.Value;

            % CONN Var
            % project name and path
            o.conn_prjPath = app.ProjectPathEditField.Value;
            o.conn_prjName = app.NameEditField.Value;
            o.subpath = app.DataPathEditField.Value;
            o.conn_rawdata = false;
            o.conn_prepsteps = app.PrepStepsEditField.Value;
            if isempty(o.conn_prepsteps), o.conn_prepsteps = ''; end

            % subject session round
            if app.AddROICheckBoxF.Value
                value = app.DropDown.Value;
                if value == "Result"
                    o.conn_steps.Add_Roi = false;
                    o.conn_steps.Add_RoiResult = true;
                elseif value == "Add Roi"
                    o.conn_steps.Add_Roi = true;
                    o.conn_steps.Add_RoiResult = false;
                end
            else
                o.conn_steps.Add_Roi = false;
                o.conn_steps.Add_RoiResult = false;
            end
            value = app.CONNTab.Data;
            o.conn_sub = value(:,1);
            o.conn_sub(cellfun(@isempty,o.conn_sub)) = [];
            o.conn_ses = value(:,2)';
            o.conn_ses(cellfun(@isempty,o.conn_ses)) = [];
            if isempty(o.conn_ses), o.conn_ses = ''; end
            o.conn_round = value(:,3)';
            o.conn_round(cellfun(@isempty,o.conn_round)) = [];
            if isempty(o.conn_round), o.conn_round = {''}; end
            o.conn_conditon = o.conn_round;

            % preprocessing
            o.conn_TR = app.TREditField.Value;
            o.conn_fwhm = app.fwhmEditField.Value;
            o.conn_strVres = app.StrucVresEditField.Value;
            o.conn_funcVres = app.FuncVresEditField.Value;
            o.conn_slice_order = str2num(app.sliceorderEditField.Value);
            % denoising
            o.conn_filtBand(1) = app.FiltBandLEditField.Value;
            o.conn_filtBand(2) = app.FiltBandHEditField.Value;
            % analysis
            o.conn_AnalysisName = app.AnaNameEditField.Value;
            o.conn_AnaSource = split(app.AnaSourceEditField.Value,', ')';
            % roi path
            o.conn_ROIpath = app.ROIpathEditField.Value;
            o.conn_ResultF = app.ResultFoldEditField.Value;
            % roi corordinate
            value = app.ROIcorTable.Data;
            value(cellfun(@isempty, value)) = [];
            o.ROIname = value(:,1);
            o.ROIname(cellfun(@isempty,o.ROIname)) = [];
            o.ROIname = split(o.ROIname,'_');
            if size(o.ROIname,2) == 1
                o.ROIname = o.ROIname';
            else
                o.ROIname(:,end) = [];
            end
            o.ROIname = unique(strcat(strcat(o.ROIname(:,1:end-1),'_'),o.ROIname(:,end)))';
            o.ROIname(cellfun(@isempty,o.ROIname)) = [];
            o.ROIcor = cell(size(o.ROIname));
            for i = 1:length(o.ROIname)
                idx = contains(string(value(:,1)),string(o.ROIname(i)));
                o.ROIcor{i} = [cell2mat(value(idx,2)),cell2mat(value(idx,3)),cell2mat(value(idx,4))];
            end
            
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Value changed function: stlevelCheckBox
        function stlevelCheckBoxValueChanged(app, event)
            value = app.stlevelCheckBox.Value;
            if value
                app.contrastCheckBox.Visible = "on";
                app.PosCheckBox.Visible = "on";
                app.NegCheckBox.Visible = "on";
            else
                app.contrastCheckBox.Visible = "off";
                app.PosCheckBox.Visible = "off";
                app.NegCheckBox.Visible = "off";
            end
        end

        % Value changed function: ndlevelCheckBox
        function ndlevelCheckBoxValueChanged(app, event)
            value = app.ndlevelCheckBox.Value;
            if value
                app.analysisCheckBox.Visible = "on";
                app.PlotflapmapCheckBox.Visible = "on";
            else
                app.analysisCheckBox.Visible = "off";
                app.PlotflapmapCheckBox.Visible = "off";
            end
        end

        % Value changed function: RunCheckBoxFunc
        function RunCheckBoxFuncValueChanged(app, event)
            value = app.RunCheckBoxFunc.Value;
            if value
                app.analysisCheckBox.Visible = "on";
                app.ndlevelCheckBox.Visible = "on";
                app.PlotflapmapCheckBox.Visible = "on";
                app.SmoothCheckBox.Visible = "on";
                app.contrastCheckBox.Visible = "on";
                app.NegCheckBox.Visible = "on";
                app.PosCheckBox.Visible = "on";
                app.stlevelCheckBox.Visible = "on";
                app.CorigistertoT1CheckBox.Visible = "on";
                app.SliceTimingCheckBox.Visible = "on";
                app.RealignmentCheckBox.Visible = "on";
                app.ReslicetoSUITCheckBox.Visible = "on";
            else
                app.analysisCheckBox.Visible = "off";
                app.ndlevelCheckBox.Visible = "off";
                app.PlotflapmapCheckBox.Visible = "off";
                app.SmoothCheckBox.Visible = "off";
                app.contrastCheckBox.Visible = "off";
                app.NegCheckBox.Visible = "off";
                app.PosCheckBox.Visible = "off";
                app.stlevelCheckBox.Visible = "off";
                app.CorigistertoT1CheckBox.Visible = "off";
                app.SliceTimingCheckBox.Visible = "off";
                app.RealignmentCheckBox.Visible = "off";
                app.ReslicetoSUITCheckBox.Visible = "off";
            end
        end

        % Value changed function: RunCheckBoxAna
        function RunCheckBoxAnaValueChanged(app, event)
            value = app.RunCheckBoxAna.Value;
            if value
                app.NormalizeCheckBox.Visible = "on";
                app.SegmentandIsolateCheckBox.Visible = "on";
            else
                app.NormalizeCheckBox.Visible = "off";
                app.SegmentandIsolateCheckBox.Visible = "off";
            end
        end

        % Value changed function: progressDropDown
        function progressDropDownValueChanged(app, event)
            value = app.progressDropDown.Value;
            switch value
                case "Raw"
                    app.RawTabGroup.Visible = "on";
                    app.ConnTabGroup.Visible = "off";
                    app.localCheckBox.Visible = "on";
                    app.behaveCheckBox.Visible = "on";
                    app.AnatomicalProgressPanel.Visible = "on";
                    app.FunctionalProgressPanel.Visible = "on";
                    app.MRISetupCheckBoxF.Visible = "on";
                    app.OverWriteCheckBoxF.Visible = "on";
                    app.CONNSTEPSPanel.Visible = "off";
                case "CONN"
                    app.RawTabGroup.Visible = "off";
                    app.ConnTabGroup.Visible = "on";
                    app.localCheckBox.Visible = "off";
                    app.behaveCheckBox.Visible = "off";
                    app.AnatomicalProgressPanel.Visible = "off";
                    app.FunctionalProgressPanel.Visible = "off";
                    app.MRISetupCheckBoxF.Visible = "off";
                    app.OverWriteCheckBoxF.Visible = "off";
                    app.CONNSTEPSPanel.Visible = "on";
            end
        end

        % Value changed function: SetupCheckBoxF
        function SetupCheckBoxFValueChanged(app, event)
            value = app.SetupCheckBoxF.Value;
            if value
                app.TREditField.Visible = "on";
                app.TREditField_2Label.Visible = "on";
            else
                app.TREditField.Visible = "off";
                app.TREditField_2Label.Visible = "off";
            end
        end

        % Value changed function: PreprocessingCheckBoxF
        function PreprocessingCheckBoxFValueChanged(app, event)
            value = app.PreprocessingCheckBoxF.Value;
            if value
                app.sliceorderEditField.Visible = "on";
                app.sliceorderEditFieldLabel.Visible = "on";
                app.fwhmEditField.Visible = "on";
                app.fwhmEditField_2Label.Visible = "on";
                app.StrucVresEditField.Visible = "on";
                app.StrucVresEditField_2Label.Visible = "on";
                app.FuncVresEditField.Visible = "on";
                app.FuncVresEditField_2Label.Visible = "on";
                app.StrucVresEditField.Visible = "on";
                app.StrucVresEditField_2Label.Visible = "on";
                app.PrepStepsEditField.Visible = "on";
                app.PrepStepsEditFieldLabel.Visible = "on";
            else
                app.sliceorderEditField.Visible = "off";
                app.sliceorderEditFieldLabel.Visible = "off";
                app.fwhmEditField.Visible = "off";
                app.fwhmEditField_2Label.Visible = "off";
                app.StrucVresEditField.Visible = "off";
                app.StrucVresEditField_2Label.Visible = "off";
                app.FuncVresEditField.Visible = "off";
                app.FuncVresEditField_2Label.Visible = "off";
                app.StrucVresEditField.Visible = "off";
                app.StrucVresEditField_2Label.Visible = "off";
                app.PrepStepsEditField.Visible = "off";
                app.PrepStepsEditFieldLabel.Visible = "off";
            end
        end

        % Value changed function: stlevelCheckBoxF
        function stlevelCheckBoxFValueChanged(app, event)
            value = app.stlevelCheckBoxF.Value;
            if value
                app.AnaNameEditField.Visible = "on";
                app.AnaNameEditFieldLabel.Visible = "on";
                app.AnaSourceEditField.Visible = "on";
                app.AnaSourceEditFieldLabel.Visible = "on";
            else
                app.AnaNameEditField.Visible = "off";
                app.AnaNameEditFieldLabel.Visible = "off";
                app.AnaSourceEditField.Visible = "off";
                app.AnaSourceEditFieldLabel.Visible = "off";
            end
        end

        % Value changed function: DenoisingCheckBoxF
        function DenoisingCheckBoxFValueChanged(app, event)
            value = app.DenoisingCheckBoxF.Value;
            if value
                app.FiltBandLEditField.Visible = "on";
                app.FiltBandLEditFieldLabel.Visible = "on";
                app.FiltBandHEditField.Visible = "on";
                app.FiltBandHEditFieldLabel.Visible = "on";
            else
                app.FiltBandLEditField.Visible = "off";
                app.FiltBandLEditFieldLabel.Visible = "off";
                app.FiltBandHEditField.Visible = "off";
                app.FiltBandHEditFieldLabel.Visible = "off";
            end
        end

        % Value changed function: localCheckBox
        function localCheckBoxValueChanged(app, event)
            value = app.localCheckBox.Value;
            if value 
                app.FTPServerDropDown.Visible = "off";
                app.FTPServerDropDownLabel.Visible = "off";
                app.FTPinfo.Visible = "off";
            else
                app.FTPServerDropDown.Visible = "on";
                app.FTPServerDropDownLabel.Visible = "on";
                app.FTPinfo.Visible = "on";
            end
        end

        % Drop down opening function: FTPServerDropDown
        function FTPServerDropDownOpening(app, event)
            value = app.FTPinfo.Value;
            field = app.FTPServerDropDown.Value;
            app.FTPServerDropDown.UserData.(field) = value;
        end

        % Value changed function: FTPServerDropDown
        function FTPServerDropDownValueChanged(app, event)
            value = app.FTPServerDropDown.Value;
            app.FTPinfo.Value = app.FTPServerDropDown.UserData.(value);
        end

        % Value changed function: behaveCheckBox
        function behaveCheckBoxValueChanged(app, event)
            value = app.behaveCheckBox.Value;
            if value 
                app.BehaveFileinfoPanel.Visible = "on";
            else
                app.BehaveFileinfoPanel.Visible = "off";
            end
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            [f,o] = SetupVar(app);
            assignin('base',"f",f);
            assignin('base',"o",o);
            assignin('base',"UIFlag",true)
            % put value to UIFigure UserData
            delete(app)
        end

        % Value changed function: RawDataCheckBox
        function RawDataCheckBoxValueChanged(app, event)
            value = app.RawDataCheckBox.Value;
            app.localCheckBox.Value = ~value;
            localCheckBoxValueChanged(app);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 711 472];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create ConnTabGroup
            app.ConnTabGroup = uitabgroup(app.UIFigure);
            app.ConnTabGroup.Position = [318 6 389 464];

            % Create CONNInitialTab
            app.CONNInitialTab = uitab(app.ConnTabGroup);
            app.CONNInitialTab.Title = 'CONN Initial';

            % Create CONNTab
            app.CONNTab = uitable(app.CONNInitialTab);
            app.CONNTab.ColumnName = {'Subject'; 'Session'; 'Round'};
            app.CONNTab.RowName = {};
            app.CONNTab.ColumnEditable = true;
            app.CONNTab.Position = [182 200 190 181];

            % Create ProjectPathEditField_2Label
            app.ProjectPathEditField_2Label = uilabel(app.CONNInitialTab);
            app.ProjectPathEditField_2Label.HorizontalAlignment = 'right';
            app.ProjectPathEditField_2Label.Position = [1 412 70 22];
            app.ProjectPathEditField_2Label.Text = 'Project Path';

            % Create ProjectPathEditField
            app.ProjectPathEditField = uieditfield(app.CONNInitialTab, 'text');
            app.ProjectPathEditField.Position = [74 412 144 22];

            % Create TREditField_2Label
            app.TREditField_2Label = uilabel(app.CONNInitialTab);
            app.TREditField_2Label.HorizontalAlignment = 'right';
            app.TREditField_2Label.Position = [41 359 25 22];
            app.TREditField_2Label.Text = 'TR';

            % Create TREditField
            app.TREditField = uieditfield(app.CONNInitialTab, 'numeric');
            app.TREditField.Position = [73 359 100 22];

            % Create fwhmEditField_2Label
            app.fwhmEditField_2Label = uilabel(app.CONNInitialTab);
            app.fwhmEditField_2Label.HorizontalAlignment = 'right';
            app.fwhmEditField_2Label.Position = [33 306 34 22];
            app.fwhmEditField_2Label.Text = 'fwhm';

            % Create fwhmEditField
            app.fwhmEditField = uieditfield(app.CONNInitialTab, 'numeric');
            app.fwhmEditField.Position = [73 306 100 22];

            % Create StrucVresEditField_2Label
            app.StrucVresEditField_2Label = uilabel(app.CONNInitialTab);
            app.StrucVresEditField_2Label.HorizontalAlignment = 'right';
            app.StrucVresEditField_2Label.Position = [3 280 64 22];
            app.StrucVresEditField_2Label.Text = 'Struc. Vres';

            % Create StrucVresEditField
            app.StrucVresEditField = uieditfield(app.CONNInitialTab, 'numeric');
            app.StrucVresEditField.Position = [73 280 100 22];

            % Create FuncVresEditField_2Label
            app.FuncVresEditField_2Label = uilabel(app.CONNInitialTab);
            app.FuncVresEditField_2Label.HorizontalAlignment = 'right';
            app.FuncVresEditField_2Label.Position = [4 253 63 22];
            app.FuncVresEditField_2Label.Text = 'Func. Vres';

            % Create FuncVresEditField
            app.FuncVresEditField = uieditfield(app.CONNInitialTab, 'numeric');
            app.FuncVresEditField.Position = [73 253 100 22];

            % Create FiltBandLEditFieldLabel
            app.FiltBandLEditFieldLabel = uilabel(app.CONNInitialTab);
            app.FiltBandLEditFieldLabel.HorizontalAlignment = 'right';
            app.FiltBandLEditFieldLabel.Position = [0 226 70 22];
            app.FiltBandLEditFieldLabel.Text = 'Filt. Band(L)';

            % Create FiltBandLEditField
            app.FiltBandLEditField = uieditfield(app.CONNInitialTab, 'numeric');
            app.FiltBandLEditField.Position = [73 226 100 22];

            % Create NameEditField_2Label
            app.NameEditField_2Label = uilabel(app.CONNInitialTab);
            app.NameEditField_2Label.HorizontalAlignment = 'right';
            app.NameEditField_2Label.Position = [218 413 40 22];
            app.NameEditField_2Label.Text = ' Name';

            % Create NameEditField
            app.NameEditField = uieditfield(app.CONNInitialTab, 'text');
            app.NameEditField.Position = [261 413 116 22];

            % Create DataPathEditFieldLabel
            app.DataPathEditFieldLabel = uilabel(app.CONNInitialTab);
            app.DataPathEditFieldLabel.HorizontalAlignment = 'right';
            app.DataPathEditFieldLabel.Position = [16 386 55 22];
            app.DataPathEditFieldLabel.Text = 'DataPath';

            % Create DataPathEditField
            app.DataPathEditField = uieditfield(app.CONNInitialTab, 'text');
            app.DataPathEditField.Position = [75 386 302 22];

            % Create ROIpathEditFieldLabel
            app.ROIpathEditFieldLabel = uilabel(app.CONNInitialTab);
            app.ROIpathEditFieldLabel.HorizontalAlignment = 'right';
            app.ROIpathEditFieldLabel.Position = [13 122 53 22];
            app.ROIpathEditFieldLabel.Text = 'ROI path';

            % Create ROIpathEditField
            app.ROIpathEditField = uieditfield(app.CONNInitialTab, 'text');
            app.ROIpathEditField.Position = [73 122 100 21];

            % Create AnaNameEditFieldLabel
            app.AnaNameEditFieldLabel = uilabel(app.CONNInitialTab);
            app.AnaNameEditFieldLabel.HorizontalAlignment = 'right';
            app.AnaNameEditFieldLabel.Position = [4 174 62 22];
            app.AnaNameEditFieldLabel.Text = 'Ana Name';

            % Create AnaNameEditField
            app.AnaNameEditField = uieditfield(app.CONNInitialTab, 'text');
            app.AnaNameEditField.Position = [73 174 100 21];

            % Create AnaSourceEditFieldLabel
            app.AnaSourceEditFieldLabel = uilabel(app.CONNInitialTab);
            app.AnaSourceEditFieldLabel.HorizontalAlignment = 'right';
            app.AnaSourceEditFieldLabel.Position = [0 148 68 22];
            app.AnaSourceEditFieldLabel.Text = 'Ana Source';

            % Create AnaSourceEditField
            app.AnaSourceEditField = uieditfield(app.CONNInitialTab, 'text');
            app.AnaSourceEditField.Position = [73 148 100 22];

            % Create FiltBandHEditFieldLabel
            app.FiltBandHEditFieldLabel = uilabel(app.CONNInitialTab);
            app.FiltBandHEditFieldLabel.HorizontalAlignment = 'right';
            app.FiltBandHEditFieldLabel.Position = [-1 200 72 22];
            app.FiltBandHEditFieldLabel.Text = 'Filt. Band(H)';

            % Create FiltBandHEditField
            app.FiltBandHEditField = uieditfield(app.CONNInitialTab, 'numeric');
            app.FiltBandHEditField.Position = [73 200 100 22];

            % Create sliceorderEditFieldLabel
            app.sliceorderEditFieldLabel = uilabel(app.CONNInitialTab);
            app.sliceorderEditFieldLabel.HorizontalAlignment = 'right';
            app.sliceorderEditFieldLabel.Position = [6 333 60 22];
            app.sliceorderEditFieldLabel.Text = 'slice order';

            % Create sliceorderEditField
            app.sliceorderEditField = uieditfield(app.CONNInitialTab, 'text');
            app.sliceorderEditField.Position = [73 333 100 22];

            % Create ResultFoldEditFieldLabel
            app.ResultFoldEditFieldLabel = uilabel(app.CONNInitialTab);
            app.ResultFoldEditFieldLabel.HorizontalAlignment = 'right';
            app.ResultFoldEditFieldLabel.Position = [3 96 66 22];
            app.ResultFoldEditFieldLabel.Text = 'Result Fold';

            % Create ResultFoldEditField
            app.ResultFoldEditField = uieditfield(app.CONNInitialTab, 'text');
            app.ResultFoldEditField.Position = [73 96 100 22];

            % Create OldprojectCheckBox
            app.OldprojectCheckBox = uicheckbox(app.CONNInitialTab);
            app.OldprojectCheckBox.Text = 'Old project';
            app.OldprojectCheckBox.Position = [238 167 80 22];

            % Create ROIcorTable
            app.ROIcorTable = uitable(app.CONNInitialTab);
            app.ROIcorTable.ColumnName = {'name'; 'Coord_X'; 'Coord_Y'; 'Coord_Z'};
            app.ROIcorTable.RowName = {};
            app.ROIcorTable.ColumnEditable = true;
            app.ROIcorTable.Position = [183 22 190 141];

            % Create PrepStepsEditFieldLabel
            app.PrepStepsEditFieldLabel = uilabel(app.CONNInitialTab);
            app.PrepStepsEditFieldLabel.HorizontalAlignment = 'right';
            app.PrepStepsEditFieldLabel.Position = [6 70 64 22];
            app.PrepStepsEditFieldLabel.Text = 'Prep.Steps';

            % Create PrepStepsEditField
            app.PrepStepsEditField = uieditfield(app.CONNInitialTab, 'text');
            app.PrepStepsEditField.Position = [73 70 100 22];

            % Create RawTabGroup
            app.RawTabGroup = uitabgroup(app.UIFigure);
            app.RawTabGroup.Position = [319 6 389 464];

            % Create InitialTab
            app.InitialTab = uitab(app.RawTabGroup);
            app.InitialTab.Title = 'Initial';

            % Create FTPServerDropDownLabel
            app.FTPServerDropDownLabel = uilabel(app.InitialTab);
            app.FTPServerDropDownLabel.HorizontalAlignment = 'right';
            app.FTPServerDropDownLabel.Position = [13 413 63 22];
            app.FTPServerDropDownLabel.Text = 'FTPServer';

            % Create FTPServerDropDown
            app.FTPServerDropDown = uidropdown(app.InitialTab);
            app.FTPServerDropDown.Items = {'ip', 'account', 'password', 'infolder', 'outfolder'};
            app.FTPServerDropDown.DropDownOpeningFcn = createCallbackFcn(app, @FTPServerDropDownOpening, true);
            app.FTPServerDropDown.ValueChangedFcn = createCallbackFcn(app, @FTPServerDropDownValueChanged, true);
            app.FTPServerDropDown.Position = [80 412 90 24];
            app.FTPServerDropDown.Value = 'password';

            % Create FTPinfo
            app.FTPinfo = uieditfield(app.InitialTab, 'text');
            app.FTPinfo.Position = [181 414 196 22];

            % Create UITableRAW
            app.UITableRAW = uitable(app.InitialTab);
            app.UITableRAW.ColumnName = {'Subject'; 'Session'; 'Round'};
            app.UITableRAW.RowName = {};
            app.UITableRAW.ColumnEditable = true;
            app.UITableRAW.Position = [188 128 190 280];

            % Create BehaveFileinfoPanel
            app.BehaveFileinfoPanel = uipanel(app.InitialTab);
            app.BehaveFileinfoPanel.Title = 'Behave File info';
            app.BehaveFileinfoPanel.Position = [7 193 178 105];

            % Create subpathEditField_4Label
            app.subpathEditField_4Label = uilabel(app.BehaveFileinfoPanel);
            app.subpathEditField_4Label.HorizontalAlignment = 'right';
            app.subpathEditField_4Label.Position = [18 58 48 22];
            app.subpathEditField_4Label.Text = 'subpath';

            % Create subpathEditFieldBEH
            app.subpathEditFieldBEH = uieditfield(app.BehaveFileinfoPanel, 'text');
            app.subpathEditFieldBEH.Position = [70 58 100 22];

            % Create FolderNestEditField_4Label
            app.FolderNestEditField_4Label = uilabel(app.BehaveFileinfoPanel);
            app.FolderNestEditField_4Label.HorizontalAlignment = 'right';
            app.FolderNestEditField_4Label.Position = [1 32 64 22];
            app.FolderNestEditField_4Label.Text = 'FolderNest';

            % Create FolderNestEditFieldBEH
            app.FolderNestEditFieldBEH = uieditfield(app.BehaveFileinfoPanel, 'text');
            app.FolderNestEditFieldBEH.Position = [70 32 100 22];

            % Create ExtentionEditField_2Label
            app.ExtentionEditField_2Label = uilabel(app.BehaveFileinfoPanel);
            app.ExtentionEditField_2Label.HorizontalAlignment = 'right';
            app.ExtentionEditField_2Label.Position = [10 7 55 22];
            app.ExtentionEditField_2Label.Text = 'Extention';

            % Create ExtentionEditFieldBEH
            app.ExtentionEditFieldBEH = uieditfield(app.BehaveFileinfoPanel, 'text');
            app.ExtentionEditFieldBEH.Position = [70 7 100 22];

            % Create TextArea_3
            app.TextArea_3 = uitextarea(app.InitialTab);
            app.TextArea_3.Position = [5 6 372 116];

            % Create subpathEditField_3Label
            app.subpathEditField_3Label = uilabel(app.InitialTab);
            app.subpathEditField_3Label.HorizontalAlignment = 'right';
            app.subpathEditField_3Label.Position = [24 384 48 22];
            app.subpathEditField_3Label.Text = 'subpath';

            % Create subpathEditField
            app.subpathEditField = uieditfield(app.InitialTab, 'text');
            app.subpathEditField.Position = [75 384 101 22];

            % Create FolderNestEditField_3Label
            app.FolderNestEditField_3Label = uilabel(app.InitialTab);
            app.FolderNestEditField_3Label.HorizontalAlignment = 'right';
            app.FolderNestEditField_3Label.Position = [6 358 64 22];
            app.FolderNestEditField_3Label.Text = 'FolderNest';

            % Create FolderNestEditField
            app.FolderNestEditField = uieditfield(app.InitialTab, 'text');
            app.FolderNestEditField.Position = [75 358 101 22];

            % Create ChckImgOpEditField_2Label
            app.ChckImgOpEditField_2Label = uilabel(app.InitialTab);
            app.ChckImgOpEditField_2Label.HorizontalAlignment = 'right';
            app.ChckImgOpEditField_2Label.Position = [3 332 68 22];
            app.ChckImgOpEditField_2Label.Text = 'ChckImgOp';

            % Create ChckImgOpEditField
            app.ChckImgOpEditField = uieditfield(app.InitialTab, 'text');
            app.ChckImgOpEditField.Position = [75 332 101 22];

            % Create condNameEditFieldLabel
            app.condNameEditFieldLabel = uilabel(app.InitialTab);
            app.condNameEditFieldLabel.HorizontalAlignment = 'right';
            app.condNameEditFieldLabel.Position = [8 305 63 22];
            app.condNameEditFieldLabel.Text = 'condName';

            % Create condNameEditField
            app.condNameEditField = uieditfield(app.InitialTab, 'text');
            app.condNameEditField.Position = [76 305 100 22];

            % Create PROCESSFLAGSPanel
            app.PROCESSFLAGSPanel = uipanel(app.UIFigure);
            app.PROCESSFLAGSPanel.Title = 'PROCESS FLAGS';
            app.PROCESSFLAGSPanel.Position = [2 7 312 465];

            % Create DataprogressPanel
            app.DataprogressPanel = uipanel(app.PROCESSFLAGSPanel);
            app.DataprogressPanel.Title = 'Data progress';
            app.DataprogressPanel.Position = [7 369 294 69];

            % Create localCheckBox
            app.localCheckBox = uicheckbox(app.DataprogressPanel);
            app.localCheckBox.ValueChangedFcn = createCallbackFcn(app, @localCheckBoxValueChanged, true);
            app.localCheckBox.Text = 'local';
            app.localCheckBox.Position = [97 25 78 22];

            % Create behaveCheckBox
            app.behaveCheckBox = uicheckbox(app.DataprogressPanel);
            app.behaveCheckBox.ValueChangedFcn = createCallbackFcn(app, @behaveCheckBoxValueChanged, true);
            app.behaveCheckBox.Text = 'behave';
            app.behaveCheckBox.Position = [149 25 78 22];

            % Create progressDropDownLabel
            app.progressDropDownLabel = uilabel(app.DataprogressPanel);
            app.progressDropDownLabel.HorizontalAlignment = 'right';
            app.progressDropDownLabel.Position = [11 28 60 22];
            app.progressDropDownLabel.Text = 'progress';

            % Create progressDropDown
            app.progressDropDown = uidropdown(app.DataprogressPanel);
            app.progressDropDown.Items = {'Raw', 'CONN'};
            app.progressDropDown.ValueChangedFcn = createCallbackFcn(app, @progressDropDownValueChanged, true);
            app.progressDropDown.Position = [14 6 74 22];
            app.progressDropDown.Value = 'Raw';

            % Create RawDataCheckBox
            app.RawDataCheckBox = uicheckbox(app.DataprogressPanel);
            app.RawDataCheckBox.ValueChangedFcn = createCallbackFcn(app, @RawDataCheckBoxValueChanged, true);
            app.RawDataCheckBox.Text = 'RawData';
            app.RawDataCheckBox.Position = [212 26 78 22];

            % Create GetDataCheckBox
            app.GetDataCheckBox = uicheckbox(app.DataprogressPanel);
            app.GetDataCheckBox.Text = 'GetData';
            app.GetDataCheckBox.Position = [212 4 67 22];

            % Create subidentifyEditFieldLabel
            app.subidentifyEditFieldLabel = uilabel(app.DataprogressPanel);
            app.subidentifyEditFieldLabel.HorizontalAlignment = 'right';
            app.subidentifyEditFieldLabel.Position = [89 2 66 22];
            app.subidentifyEditFieldLabel.Text = 'sub identify';

            % Create subidentifyEditField
            app.subidentifyEditField = uieditfield(app.DataprogressPanel, 'text');
            app.subidentifyEditField.Position = [159 3 48 22];

            % Create AnatomicalProgressPanel
            app.AnatomicalProgressPanel = uipanel(app.PROCESSFLAGSPanel);
            app.AnatomicalProgressPanel.Title = 'Anatomical Progress';
            app.AnatomicalProgressPanel.Position = [7 281 141 79];

            % Create RunCheckBoxAna
            app.RunCheckBoxAna = uicheckbox(app.AnatomicalProgressPanel);
            app.RunCheckBoxAna.ValueChangedFcn = createCallbackFcn(app, @RunCheckBoxAnaValueChanged, true);
            app.RunCheckBoxAna.Text = 'Run';
            app.RunCheckBoxAna.Position = [5 35 44 22];

            % Create SegmentandIsolateCheckBox
            app.SegmentandIsolateCheckBox = uicheckbox(app.AnatomicalProgressPanel);
            app.SegmentandIsolateCheckBox.Text = 'Segment and Isolate';
            app.SegmentandIsolateCheckBox.Position = [5 19 132 22];

            % Create NormalizeCheckBox
            app.NormalizeCheckBox = uicheckbox(app.AnatomicalProgressPanel);
            app.NormalizeCheckBox.Text = 'Normalize';
            app.NormalizeCheckBox.Position = [5 2 76 22];

            % Create FunctionalProgressPanel
            app.FunctionalProgressPanel = uipanel(app.PROCESSFLAGSPanel);
            app.FunctionalProgressPanel.Title = 'Functional Progress';
            app.FunctionalProgressPanel.Position = [155 114 141 246];

            % Create RunCheckBoxFunc
            app.RunCheckBoxFunc = uicheckbox(app.FunctionalProgressPanel);
            app.RunCheckBoxFunc.ValueChangedFcn = createCallbackFcn(app, @RunCheckBoxFuncValueChanged, true);
            app.RunCheckBoxFunc.Text = 'Run';
            app.RunCheckBoxFunc.Position = [5 203 44 22];

            % Create RealignmentCheckBox
            app.RealignmentCheckBox = uicheckbox(app.FunctionalProgressPanel);
            app.RealignmentCheckBox.Text = 'Realignment';
            app.RealignmentCheckBox.Position = [5 186 89 22];

            % Create SliceTimingCheckBox
            app.SliceTimingCheckBox = uicheckbox(app.FunctionalProgressPanel);
            app.SliceTimingCheckBox.Text = 'Slice Timing';
            app.SliceTimingCheckBox.Position = [5 169 87 22];

            % Create CorigistertoT1CheckBox
            app.CorigistertoT1CheckBox = uicheckbox(app.FunctionalProgressPanel);
            app.CorigistertoT1CheckBox.Text = 'Corigister to T1';
            app.CorigistertoT1CheckBox.Position = [5 152 104 22];

            % Create stlevelCheckBox
            app.stlevelCheckBox = uicheckbox(app.FunctionalProgressPanel);
            app.stlevelCheckBox.ValueChangedFcn = createCallbackFcn(app, @stlevelCheckBoxValueChanged, true);
            app.stlevelCheckBox.Text = '1st level';
            app.stlevelCheckBox.Position = [5 135 66 22];

            % Create PosCheckBox
            app.PosCheckBox = uicheckbox(app.FunctionalProgressPanel);
            app.PosCheckBox.Visible = 'off';
            app.PosCheckBox.Text = 'Pos';
            app.PosCheckBox.Position = [27 103 43 22];

            % Create NegCheckBox
            app.NegCheckBox = uicheckbox(app.FunctionalProgressPanel);
            app.NegCheckBox.Visible = 'off';
            app.NegCheckBox.Text = 'Neg';
            app.NegCheckBox.Position = [27 86 44 22];

            % Create contrastCheckBox
            app.contrastCheckBox = uicheckbox(app.FunctionalProgressPanel);
            app.contrastCheckBox.Visible = 'off';
            app.contrastCheckBox.Text = 'contrast';
            app.contrastCheckBox.Position = [27 120 65 22];

            % Create SmoothCheckBox
            app.SmoothCheckBox = uicheckbox(app.FunctionalProgressPanel);
            app.SmoothCheckBox.Text = 'Smooth';
            app.SmoothCheckBox.Position = [6 52 63 22];

            % Create ndlevelCheckBox
            app.ndlevelCheckBox = uicheckbox(app.FunctionalProgressPanel);
            app.ndlevelCheckBox.ValueChangedFcn = createCallbackFcn(app, @ndlevelCheckBoxValueChanged, true);
            app.ndlevelCheckBox.Text = '2nd level';
            app.ndlevelCheckBox.Position = [6 35 70 22];

            % Create analysisCheckBox
            app.analysisCheckBox = uicheckbox(app.FunctionalProgressPanel);
            app.analysisCheckBox.Visible = 'off';
            app.analysisCheckBox.Text = 'analysis';
            app.analysisCheckBox.Position = [25 18 65 22];

            % Create ReslicetoSUITCheckBox
            app.ReslicetoSUITCheckBox = uicheckbox(app.FunctionalProgressPanel);
            app.ReslicetoSUITCheckBox.Text = 'Reslice to SUIT';
            app.ReslicetoSUITCheckBox.Position = [6 69 105 22];

            % Create PlotflapmapCheckBox
            app.PlotflapmapCheckBox = uicheckbox(app.FunctionalProgressPanel);
            app.PlotflapmapCheckBox.Visible = 'off';
            app.PlotflapmapCheckBox.Text = 'Plot flapmap';
            app.PlotflapmapCheckBox.Position = [25 0 89 22];

            % Create TextArea
            app.TextArea = uitextarea(app.PROCESSFLAGSPanel);
            app.TextArea.Position = [12 9 289 98];

            % Create OverWriteCheckBoxF
            app.OverWriteCheckBoxF = uicheckbox(app.PROCESSFLAGSPanel);
            app.OverWriteCheckBoxF.Text = 'OverWrite';
            app.OverWriteCheckBoxF.Position = [72 131 76 22];

            % Create MRISetupCheckBoxF
            app.MRISetupCheckBoxF = uicheckbox(app.PROCESSFLAGSPanel);
            app.MRISetupCheckBoxF.Text = 'MRI Setup';
            app.MRISetupCheckBoxF.Position = [72 152 79 22];

            % Create CONNSTEPSPanel
            app.CONNSTEPSPanel = uipanel(app.PROCESSFLAGSPanel);
            app.CONNSTEPSPanel.Title = 'CONN STEPS';
            app.CONNSTEPSPanel.Position = [9 166 139 192];

            % Create SetupCheckBoxF
            app.SetupCheckBoxF = uicheckbox(app.CONNSTEPSPanel);
            app.SetupCheckBoxF.ValueChangedFcn = createCallbackFcn(app, @SetupCheckBoxFValueChanged, true);
            app.SetupCheckBoxF.Text = 'Setup';
            app.SetupCheckBoxF.Position = [42 144 53 22];

            % Create PreprocessingCheckBoxF
            app.PreprocessingCheckBoxF = uicheckbox(app.CONNSTEPSPanel);
            app.PreprocessingCheckBoxF.ValueChangedFcn = createCallbackFcn(app, @PreprocessingCheckBoxFValueChanged, true);
            app.PreprocessingCheckBoxF.Text = 'Preprocessing';
            app.PreprocessingCheckBoxF.Position = [20 121 99 22];

            % Create DenoisingCheckBoxF
            app.DenoisingCheckBoxF = uicheckbox(app.CONNSTEPSPanel);
            app.DenoisingCheckBoxF.ValueChangedFcn = createCallbackFcn(app, @DenoisingCheckBoxFValueChanged, true);
            app.DenoisingCheckBoxF.Text = 'Denoising';
            app.DenoisingCheckBoxF.Position = [32 99 75 22];

            % Create stlevelCheckBoxF
            app.stlevelCheckBoxF = uicheckbox(app.CONNSTEPSPanel);
            app.stlevelCheckBoxF.ValueChangedFcn = createCallbackFcn(app, @stlevelCheckBoxFValueChanged, true);
            app.stlevelCheckBoxF.Text = '1st level';
            app.stlevelCheckBoxF.Position = [37 75 66 22];

            % Create ndlevelCheckBoxF
            app.ndlevelCheckBoxF = uicheckbox(app.CONNSTEPSPanel);
            app.ndlevelCheckBoxF.Text = '2nd level';
            app.ndlevelCheckBoxF.Position = [34 52 70 22];

            % Create AddROICheckBoxF
            app.AddROICheckBoxF = uicheckbox(app.CONNSTEPSPanel);
            app.AddROICheckBoxF.Text = 'Add ROI';
            app.AddROICheckBoxF.Position = [35 27 68 22];

            % Create DropDown
            app.DropDown = uidropdown(app.CONNSTEPSPanel);
            app.DropDown.Items = {'Result', 'Add roi'};
            app.DropDown.Position = [58 5 78 22];
            app.DropDown.Value = 'Result';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = SetupUI(varargin)
            if ~nargin
                f = [];
                o = [];
            end
            if nargin > 0
                f = varargin{1};
            end
            if nargin > 1
                o = varargin{2};
            end
           
            % Create UIFigure and components
            createComponents(app)
            
            % Initial components value
            InitialVar(app,f,o);
            InitialVisible(app)
            % Register the app with App Designer
            registerApp(app, app.UIFigure)

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