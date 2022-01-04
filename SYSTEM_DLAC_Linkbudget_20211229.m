%% 2020-02-28, Image Rejection calculation ?
%% 2020-02-29, Apply Receiver BPF before downconversion
%% 2020-02-29, QEC Estimation and Correction
%% 2020-02-28, Phase Channel Offset correct: R.2.1.1. and R.2.1.2.
%% 2020-03-03, Apply Mutli-Stage of FIR filter
%% 2020-03-03, How to configure the FIR paramters by Norder, bw tolerance, Window and MStage ?
%% 2020-03-03, Why the filitering by 'conv' will Not generate the time shift ?
%% 2020-03-03, Analysis receiver parameter with perfect txwaveform ?
%% 2020-03-03, Add DC or IM2 to waveform
%% 2020-03-04, IQ waveform QEC correct for DL part ?
%% 2020-03-04, Important !!!, 'IFFT' will add the NoiseFloor out-of-band
%% 2020-03-04, Important !!!, Place the FIR Filter after signal did the 'IFFT'. (%% B.3.1. Apply FIR Filter)
%% 2020-03-05, Evaluate the DC impact in receiver ?
%% 2020-03-05, Issue1, large fs %% B.1.0. Reduce the Lengths of Carrier
%% 2020-03-05, TOPIC1, The Linearity Relationship between fs and FIR filter Norder !!! ?
%% 2020-03-10, Issue2, If without NCO shift firstly, the QEC estimation will Not correct ?
%% 2020-03-10, TOPIC2, The Comparsion the fNCO(1,5,10,20MHz...) VS QEC result ?
%% 2020-03-10, Issue3, offset_mea is issue to EVM calculation ?
%% 2020-03-11, T.5.4. QEC for the Imbalance LO
%% 2020-03-12, Choose the HighFrequencyBand(fLO+fBB) for QEC to avoid the spectrum alising in LowFrequencyBand
%% 2020-03-12, Imbalance LO will cause the ImageCarrier, this ImageCarrier will impact the ACLR performance
%% 2020-03-12, Imbalance LO will cause the ImageCarrier(located at outside of band), when flag_NCO ON, this ImageCarrier will impact the ACLR performance
%% 2020-03-12, Imbalance LO will cause the ImageCarrier(located at inside of band), when flag_NCO OFF, this ImageCarrier will impact the EVM performance
%% 2020-03-14, If PN_LOU.PN_ThetaDeg==0, no-assignment
%% 2020-03-14, Add PhaseDriftDeg for AM/PM, in LO function ver. g3
%% 2020-03-15, Issue4, Ipwr relative between I and Q changed after Decimation ?
%% 2020-03-15, T.4.1. Generate Multi-Carriers, check impact to ACLR and EVM ?
%% 2020-03-16, Add NoiseFigure for Receiver ?
%% 2020-03-18, Update Gain_dB for Multi-Carriers
%% 2020-03-19, Update Ipwr_dB_inband_NoiseFloor, the 10*log10(BW) has already included the 10*log10(df)
%% 2020-03-26, Compare EVM result of the filtering application by 1. time-domain: conv(x,y) and 2. Freq-domain: ifft(fft(x).fft(y)) -> NO Big Different
%% 2020-03-26, choice the filtering application by Freq-domation method, 'FD'
%% 2020-05-22, filter by freq. domail will NOT introduce the delay from FIR, the second advantage NOT cause the edge of bandwidth grow to impact ACLR
%% 2020-03-26, Change all Output results to be Column type
%% 2020-03-29, Issue5, Add NOISE FIGURE
%% 2020-03-29, T.7.0.0. Add AM-PM to Gain_IM3_g3 ?
%% 2020-04-23, fNCO should based on raster df
%% 2020-04-23, Issue6, CFR will cause the NOISE FLOOR and bad ACLR
%% 2020-04-23, Issue6, For AC, compare the phase accuracy between low SNR wo CFR and higher SNR wi CFR ???
%% 2020-04-23, Issue6, compare the CFR output between different carrier type, BW/NCarriers/fNCO... ???
%% 2020-04-28, Nsamps=307200, to reduce EVM calculation time
%% 2020-05-13, Sampling Rate should larger than IBW before NCO and MultiCarrier Combination
%% 2020-05-19, Compare results of different Target power
%% 2020-05-28, %% T.5.1. Apply BPF after Upconversion and before PA
%% 2020-07-24, Add ADC performance for RX: SFDR/IIP3/NF ???
%% 2020-07-24, Add Blocking for RX
%% 2020-07-29, Extend Data Rate by Repmat Samples --> Remove!!!
%% 2020-07-31, T4e. LO with Unlinearity
%% 2020-08-01, TEST N*InterModulation, 'off', default
%% 2020-08-03, R3a. Apply IFBPF Filter after Mixer or before ADC
%% 2021-04-01, RXQEC Method Comparsion: trx_QEC_g and dsp_iqc_widely_g
%% 2021-04-09, TEST, flagT2b_TEST_ADC, ADC output, Quantization Error/Dither/NoiseShape
%% 2021-04-19, TEST, flagR3o_TEST_ADC, ADC output, Quantization Error/Dither/NoiseShape
%% 2021-04-19, TEST, flagR3o_TEST_ADC, Add ADC Quantization Error to estimate the QEC ability
%% 2021-04-22, TEST, flagC4_TEST_Delay, Apply Phase Delay to waveform

%% AC ================================================================================================================================
%% 2020-04-12, flag_Waveform_Analysis: 'DL' or 'AC'
%% 2020-04-13, C0c. Generate AC waveform
%% 2020-04-13, Method is NO impact to phase accuracy, Method will has different Spectrum results ???
%% 2020-04-13, NsampsEVM: to reduce the samples of EVM calculation for AC analysis
%% 2020-04-14, AC phase accuracy, There is NO impact for filter application by FD or td
%% 2020-04-14, AC phase accuracy, LO Imbalance has impact
%% 2020-04-14, AC phase accuracy, Imbalance mix with TX and RX will Large impact the AC phase accuracy
%% 2020-04-14, For EVM accuracy with Large offset_CH, the calculation length from ref and mea should be sufficient, how to choice the Length ???
%% 2020-04-14, Y = CIRCSHIFT(X,K,DIM)
%% 2020-04-15, AC phase correction, flag_cor_OffsetChPhs=1 will directly correct the phase difference from each branch
%% 2020-04-15, AC phase accuracy, What the LO performance impact the phase estimation accuracy ? LO IMB/PMN/SPURS ?
%% 2020-05-22, DLAC mode, DL signal is NOT for demodulation, the role is only for Interference to AC signal
%% 2021-03-10, TEST, flagD7c_TEST_ACBrGainVariation: waveformC0_AC with Gain variations for each Branch
%% 2021-03-15, Check AC Combination: Ipwr of carrier and noise floor
%% 2021-05-14, TEST, flagC6_TEST_ACBrGainVariation: waveformC0_AC with Gain variations for each Branch


if ~exist('NsampsCA','var')&&~exist('NsampsBB','var')&&~exist('NsampsTX','var')&&~exist('NsampsRX','var')
    clear all
    clc
    close all
    warning off
    format short g
end

%% 2020-05-25, Excel Read
xlsFile='Matlab_Linkbudget.xlsx';
xlsSheetName = 'AC_TXNoiseFloor';
% xlsSheetName = 'QEC';
xlsRangeArrayInput ='A1:B550';
% xlsRangeArrayOutput ='J1:J110';
xlsReadShift = [0 0];
xlsWriteShift = [0 0];
fnum_dir = 'C:\Users\egranli\Documents\MATLAB\grant\SAVE_fig';
fnum_dir = [];

xlsRangeArrayInput ='A1:E550';
xlsReadShift = [0 1];
xlsWriteShift = [0 1];

xlsRangeArrayInput ='F1:I550';
xlsReadShift = [0 0];
xlsWriteShift = [0 0];



% case1:
xlsSheetName = 'Input1';
xlsRangeArrayInput ='B1:E550';
xlsReadShift = [0 0];
xlsWriteShift = [0 0];

% case2:
xlsSheetName = 'AC_TXNoiseFloor';
xlsRangeArrayInput ='K1:N550';
xlsReadShift = [0 0];
xlsWriteShift = [0 0];

% initialization
% clear NsampsTX NsampsRX NsampsDP
% clear waveformD7_ACPhsAlign
% close all

%% Carrier Configuration ========================================================================================================================
if ~exist('NsampsCA','var') || (exist('waveformD7_ACPhsAlign','var'))
    
    %% 2020-04-12, flag_AnalysisWaveform: 'DL'\'UL'\'DL2UL'\'AC'\'DLAC'\
    %% 2020-04-12, DL: TX ONLY\UL: RX ONLY\DL2UL: FROM DL to UL\AC: AC source only\DLAC: AC+DL source
    flagC0_AnalysisWF = ExcelRead('flagC0_AnalysisWF',1,[],[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    MOD = ExcelRead('MOD',1,[],'char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift); %; % Modulation
    %% 2020-03-05, Reduce the Lengths of Carrier
    %         flagC0_SamplesDMC = 'off'; % on/off
    flagC0_SamplesDMC = ExcelRead('flagC0_SamplesDMC',1,[],'char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift); % on/off
    Nbr_DL = ExcelRead('Nbr_DL',1,[],'scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    %% 2020-05-22, filter by freq. domail will NOT introduce the delay from FIR, the second advantage NOT cause the edge of bandwidth grow to impact ACLR
    mehtod_filter = ExcelRead('mehtod_filter',1,[],'char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    % mehtod_filter = 'td';
    %% C0. Carrier settings
    if strcmp(flagC0_AnalysisWF,'DLAC')
        NCarrierDL = ExcelRead('NCarrierDL',[],'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        bw_Channelcell_DL = ExcelRead('bw_Channelcell_DL',[NCarrierDL],'column','cell',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        NCarrierAC = ExcelRead('NCarrierAC',[],'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        bw_Channelcell_AC = ExcelRead('bw_Channelcell_AC',[NCarrierAC],'column','cell',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        
        % table
        if ~isrow(bw_Channelcell_DL)
            bw_Channelcell_DLAll=bw_Channelcell_DL.';
        else
            bw_Channelcell_DLAll=bw_Channelcell_DL;
        end
        tableInput_Carrier =table(flagC0_AnalysisWF, MOD, flagC0_SamplesDMC, bw_Channelcell_DLAll, bw_Channelcell_AC);
        
    elseif strcmp(flagC0_AnalysisWF,'DL2UL')||strcmp(flagC0_AnalysisWF,'DL')
        NCarrierDL = ExcelRead('NCarrierDL',[],'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        bw_Channelcell_DL = ExcelRead('bw_Channelcell_DL',[NCarrierDL],'column','cell',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        NCarrierAC = 0;
        
        % table
        bw_Channelcell_DLAll=bw_Channelcell_DL.';
        tableInput_Carrier =table(flagC0_AnalysisWF, MOD, flagC0_SamplesDMC, bw_Channelcell_DLAll);
        
    elseif strcmp(flagC0_AnalysisWF,'AC')
        NCarrierAC = ExcelRead('NCarrierAC',[],'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        bw_Channelcell_AC = ExcelRead('bw_Channelcell_AC',[NCarrierAC],'column','cell',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        NCarrierDL = 0;
        
        % table
        tableInput_Carrier =table(flagC0_AnalysisWF, flagC0_SamplesDMC, MOD, bw_Channelcell_AC);
    end
    if NCarrierDL>=2 && ~strcmp(flagC0_SamplesDMC,'off')
        for k=1:nchoosek(NCarrierDL,2)
            bw_ChannelCompare(k) = strcmp(bw_Channelcell_DL{k,:},bw_Channelcell_DL{k+1,:});
            if k+1>nchoosek(NCarrierDL,2)
                break
            end
        end
        
        if any(bw_ChannelCompare)==0
            flagC0_SamplesDMC=4;
        end
    end
    
    %% C0a. RF settings
    IBW = ExcelRead('IBW',1,[],'scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift); %
    fRF = ExcelRead({'fRF_start','fRF_stop'},[2 1],'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift); %
    orderHD = ExcelRead('orderHD',1,[],'scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift); %
    fnum = 050601;
    % table
    IBW_MHz=IBW/1e6;
    fRF_MHz=fRF/1e6;
    tableInput_RF = table(IBW_MHz,fRF_MHz,orderHD);
    tableC0Input = table(tableInput_Carrier, tableInput_RF)
    
    disp([flagC0_AnalysisWF, ' Analysis ===================================================='])
    
    if strcmp(flagC0_AnalysisWF,'DL2UL')||strcmp(flagC0_AnalysisWF,'DL')||strcmp(flagC0_AnalysisWF,'DLAC')
        disp('Generate Carrier ====================================================')
        
        if size(bw_Channelcell_DL,1)<NCarrierDL && size(bw_Channelcell_DL,1)==1
            bw_Channelcell_DL = repmat(bw_Channelcell_DL,NCarrierDL,1);
        elseif size(bw_Channelcell_DL,1)<NCarrierDL
            error('!')
        end
        
        NCarriers = NCarrierDL+NCarrierAC;
        for idCarDL=1:NCarrierDL
            
            %% C0b. Generate DL waveform
            bw_Channelcell_DL_idC = bw_Channelcell_DL(idCarDL,:);
            tableInput(idCarDL,:) = table(idCarDL, bw_Channelcell_DL_idC, MOD, flagC0_SamplesDMC);
            
            for idBR=1:Nbr_DL
                [waveformC0_OFDM(:,idBR), bbdata_grid, DLconfig(idCarDL,idBR)] = OFDM_Mod_g(bw_Channelcell_DL_idC{1,:},MOD,flagC0_SamplesDMC,'LTE',[fnum, NCarriers, 1, NCarrierAC+idCarDL],fnum_dir);                
            end
            % export to cell
            waveformC0cell_DL{idCarDL,:}=waveformC0_OFDM;
            NameC0cell_DL{idCarDL,:}='DL';
            [~,DIMFFT] = max(size(waveformC0_OFDM));
            
            % export:
            fsDL(idCarDL,:) = DLconfig(idCarDL,:).fs;
            bwCarrierDL(idCarDL,:) = DLconfig(idCarDL,:).bwCarrier;
            bwInbandDL = bwCarrierDL/2*[-1 1];
            bwChDL(idCarDL,:) = DLconfig(idCarDL,:).bwChannel;
            NsampsC0_DL(idCarDL,:) = length(waveformC0_OFDM);
            df_DL(idCarDL,:) = fsDL(idCarDL,:)/NsampsC0_DL(idCarDL);
            %             NbrDL(idCarDL,:) = size(waveformC0_OFDM,2);
            title(['DL Carrier',num2str(idCarDL),', fs:',num2str(fsDL(idCarDL,:)/1e6),'MHz'])
            waveformC0_OFDM = [];
        end
    else
        waveformC0cell_DL=[];
        NsampsC0_DL=[];
        fsDL=[];
        bwInbandDL=[];
        NameC0cell_DL=[];
        df_DL=[];
        bwChDL=[];
        NCarriers = NCarrierDL+NCarrierAC;
        
    end % strcmp(flag_AnalysisWF,'DL2UL')||strcmp(flag_AnalysisWF,'DL')||strcmp(flag_AnalysisWF,'DLAC')
    
    %% C0c. Generate AC waveform
    if strcmp(flagC0_AnalysisWF,'AC')||strcmp(flagC0_AnalysisWF,'DLAC')
        
        if strcmp(bw_Channelcell_AC,'20MHz')
            bwChAC = 20e6;
            fsAC = 30.72e6;
            NsampsC0_AC = 19200*1;
        elseif strcmp(bw_Channelcell_AC,'100MHz')
            bwChAC = 100e6;
            fsAC = 122.88e6;
            NsampsC0_AC = 19200*5;
        end
        if exist('NsampsC0_DL','var')&&~isempty(NsampsC0_DL)
            %             if strcmp(bw_Channelcell_AC,'20MHz')
            %                 bwChAC = 20e6;
            %                 fsAC = 30.72e6;
            %             elseif strcmp(bw_Channelcell_AC,'100MHz')
            %                 bwChAC = 100e6;
            %                 fsAC = 122.88e6;
            %             end
            NsampsC0_AC = max(NsampsC0_DL)*fsAC/max(fsDL);
        elseif ~exist('NsampsC0_AC','var')||isempty(NsampsC0_AC)
            NsampsC0_AC = [];
        end
        %                     NsampsC0_AC = 19200/4; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                     NsampsC0_AC = 6144+1;
        
        %         subplot(NCarriers,1,1)
        if ~exist('waveformD7_ACPhsAlign','var')||isempty(waveformD7_ACPhsAlign)
            [waveformC0_AC, ACconfig] = AntCal_genACSource_g100(bw_Channelcell_AC{1,:}, [], NsampsC0_AC, [], [], [fnum, NCarriers, 1, 1] ,fnum_dir);
            title(['AC Carrier',', fs:',num2str(ACconfig.fs/1e6),'MHz'])
        elseif exist('ACconfig','var')&&~isempty(ACconfig)
            waveformC0_AC = waveformD7_ACPhsAlign;
            if length(waveformC0_AC)~=NsampsC0_AC
                error('!')
            end
            PLOT_FFT_dB_g(waveformD7_ACPhsAlign(:,:), ACconfig.fs, length(waveformD7_ACPhsAlign), ['ACPhsAlign'], 'df', 'full');
            title(['AC Carrier with Phase Alignment'])
        else
            error('!')
        end
        
        Debug_ACDemodulation = 'off';
        if strcmp(Debug_ACDemodulation,'on')
            
            %% T7. Add PhaseShift for each branch
            ACPhsShiftNbrDeg_T7 = ExcelRead('ACPhsShiftNbrDeg_T7',4,'row','vector',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
            waveformT7_addPhsShift = exp(1i*ACPhsShiftNbrDeg_T7./180*pi).*waveformC0_AC;
            
            %% T8. Combine all branches for AC
            waveformT8_CombACNbr = sum(waveformT7_addPhsShift,2); % all branch combination
            
            %% D7b. AC DeModulation
            ACDomod_Nshift = 0;
            fnum=fnum+1;
            [ACdmod_t0Mean, ACdmod_p0DegMean, ACdmod_phEstDeg, ACdmod_phEstDegDrift, ACdmod_dataCapCor, ACdmod_SNR, ACdmod_dataCapNbrwoPD] = AntCal_phaseDemodulateApp_g100(waveformT8_CombACNbr,ACconfig,ACDomod_Nshift,fnum);
            ACPhsEstDegError = ACdmod_p0DegMean-ACPhsShiftNbrDeg_T7;
        end
        
        % export to cell
        waveformC0cell_AC{1}=waveformC0_AC;
        NsampsC0_AC = length(waveformC0_AC);
        fsAC = ACconfig.fs;
        %         NfftAC = ACconfig.Nfft;
        %         dfSubcarrier = fsAC/NfftAC;
        bwInbandAC = ACconfig.bwInband/2*[-1 1];
        Nbr_AC = size(waveformC0_AC,2);
        df_AC = fsAC/NsampsC0_AC;
        NameC0cell_AC{1}='AC';
        [~,DIMFFT] = max(size(waveformC0_AC));
        
        if exist('waveformD7_ACPhsAlign','var')&&~isempty(waveformD7_ACPhsAlign)
            waveformC0_AC = waveformD7_ACPhsAlign;
        end
        
        
    else
        waveformC0cell_AC=[];
        NsampsC0_AC=[];
        fsAC=[];
        bwInbandAC=[];
        NameC0cell_AC=[];
        df_AC=[];
        bwChAC=[];
    end % strcmp(flag_AnalysisWF,'AC')||strcmp(flag_AnalysisWF,'DLAC')
    
    % export
    waveformC0cell=[waveformC0cell_AC;waveformC0cell_DL]; % AC signal is priority 1
    NsampsC0=[NsampsC0_AC;NsampsC0_DL];
    fsC0=[fsAC;fsDL];
    bwInbandC0=[bwInbandAC;bwInbandDL];
    NameC0cell=[NameC0cell_AC;NameC0cell_DL];
    bwChNC=[bwChAC;bwChDL];
    bwInband=bwInbandC0;
    NameCell=NameC0cell;
    df=[df_AC;df_DL];
    NCarriers = size(waveformC0cell,1);
    fs = fsC0;
    waveformCell = waveformC0cell;
    Nsamps=max(NsampsC0);
    
    % excel
    ExcelWrite(NameCell,'CarriersC0',size(NameCell),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    ExcelWrite(NCarriers,'NCarriersC0',size(NCarriers),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    ExcelWrite(bwInband,'bwInbandC0',size(bwInband),'r1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    
    %% C1. UPSampling, input to FPGA sampling rate
    flagC1_UPSamp = ExcelRead('flagC1_UPSamp',1,'row','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift); % UPSamp/off
    
    if ~strcmp(flagC1_UPSamp,'off')
        fsC1_UPS = ExcelRead('fsC1_UPS',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        ratioC1_UPS = fsC1_UPS./fs;
        % excel
        ExcelWrite(ratioC1_UPS,'ratioC1_UPS',size(ratioC1_UPS),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        
        %% C1. FIRC1UPS
        flagC1_UPSFIR = ExcelRead('flagC1_UPSFIR',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        if strcmp(flagC1_UPSFIR,'off')
            FIRc1_UPS ={1};
        else
            % FIR input
            FIRc1_Wtype = ExcelRead('FIRc1_Wtype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
            FIRc1_Ftype = ExcelRead('FIRc1_Ftype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
            FIRc1_Order = ExcelRead('FIRc1_Order',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
            FIRc1_fTolerance = ExcelRead('FIRc1_fTolerance',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
            FIRc1_K_AttdB = ExcelRead('FIRc1_K_AttdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
            FIRc1_K_fdelta = ExcelRead('FIRc1_K_fdelta',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
            FIRc1_fcutoffL = ExcelRead('FIRc1_fcutoffL',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
            FIRc1_fcutoffH = ExcelRead('FIRc1_fcutoffH',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
            
            FIRc1_UPS = SYM_FIRApp(FIRc1_Wtype,FIRc1_Ftype,FIRc1_Order,FIRc1_K_AttdB,FIRc1_K_fdelta,FIRc1_fTolerance,FIRc1_fcutoffL,FIRc1_fcutoffH,df,bwInband,NCarriers);
        end
        
        fnum = fnum+1;
        NameCellflag=repmat([{flagC1_UPSamp},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        [waveformC1cell_UPS, fsC1_UPS, NsampsC1_UPS, EVMc1_UPS, tableC1_UPS] = SYM_UPSampApp(waveformCell, ratioC1_UPS, fs, [FIRc1_UPS], mehtod_filter, bwInband, [], [], fnum, NameCelltitle, fnum_dir);
        
        % export
        waveformCell= waveformC1cell_UPS;
        fs=fsC1_UPS;
        Nsamps=NsampsC1_UPS;
        %         df=fs./Nsamps;
        
        % excel
        %         ExcelWrite(tableC1_UPS.tableFilterInput.FIR_Order,'FIRc1_Order_out',size(tableC1_UPS.tableFilterInput.FIR_Order),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(EVMc1_UPS,2),'EVMc1_UPS',size(EVMc1_UPS),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(fsC1_UPS/1e6,2),'fsC1_UPS_MHz',size(fsC1_UPS),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(Nsamps,'NsampsC1_UPS',size(NsampsC1_UPS),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        
        %         case {'ETNSamp'}
        %             %% 2020-07-29, Extend Data Rate by Repmat Samples
        %             %% C1b. Extend Data Rate
        %             ratioC1b_fsExtend = fsC1_UPS./fsC0;
        %             fnum = fnum+1;
        %             for idC=1:NCarriers
        %                 waveformC1cell_ETN{idC,:} = repmat(waveformC0cell{idC,:}, ratioC1b_fsExtend(idC,:),1);
        %                 figure(fnum)
        %                 subplot(NCarriers,1,idC)
        %                 PLOT_FFT_dB_g(waveformC1cell_ETN{idC,:}, fsC1_UPS, length(waveformC0cell{idC,:})/1, [], 'df', 'full', 'pwr', [fnum,NCarriers,1,idC]);
        %                 [IpwrdB_C0{idC}, ~, ~] = Pwr_Inband_g(fft(waveformC0cell{idC,:}, length(waveformC0cell{idC,:}), DIMFFT), fs(idC,:), bwInband(idC,:), [], 'full', 0);
        %                 [IpwrdB_C1b{idC}, ~, ~] = Pwr_Inband_g(fft(waveformC1cell_ETN{idC,:}, length(waveformC1cell_ETN{idC,:}), DIMFFT), fs(idC,:), bwInband(idC,:), [], 'full', 0);
        %                 title([NameCell{idC},newline,'Extend to fs:',num2str(fsC1_UPS/1e6),'MHz'])
        %             end
        %
        %             % export
        %             waveformCell= waveformC1cell_ETN;
        %             fs=fsC1_UPS;
        %             Nsamps=length(waveformC1cell_ETN{idC,:});
        %             df=fs./Nsamps;
        %
        %             % excel
        % %             ExcelWrite(tableC1_UPS.tableFilterInput.FIR_Order,'FIRc1_Order_out',size(tableC1_UPS.tableFilterInput.FIR_Order),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        % %             ExcelWrite(round(EVMc1_UPS,2),'EVMc1_UPS',size(EVMc1_UPS),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        %             ExcelWrite(round(fsC1_UPS/1e6,2),'fsC1_UPS_MHz',size(fsC1_UPS),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        %             ExcelWrite(Nsamps,'NsampsC1_UPS',size(NsampsC1_UPS),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        
    end % flagC1_UPSFIR
    
    %% C2. Add Digital Gain
    flagC2_addDigGain = ExcelRead('flagC2_addDigGain',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagC2_addDigGain,'off')
        % input
        GaindB_C2 = ExcelRead('GaindB_C2',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        IpwrdB_C2Target = ExcelRead('IpwrdB_C2Target',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        bwACLROffsetC2 = ExcelRead('bwACLROffsetC2',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        if fs<bwACLROffsetC2*2
            bwACLROffsetC2=[];
            flagC2_bwACLROffset=0;
        else
            flagC2_bwACLROffset=1;
        end
        %         bwACLROffsetC2 = sum(bwChNC);
        
        fnum=fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagC2_addDigGain},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        %         [waveformB3cell_wiGain, EVMb3_Gain,  tableB3_Gain] = SYM_AMPApp(waveformCell, GaindB_C2, NFdB_C2, Poip3dB_C2, Pim2dB_C2, PDCdB_C2, AMtoPMDegDrift_C2, [], fs, IpwrdB_C2Target, bwInband, bwACLROffsetC2, fnum, NameCelltitle);
        [waveformC2cell_wiGain, EVMc2_Gain,  tableC2_Gain] = SYM_AMPApp(waveformCell, GaindB_C2, [], [], [], [], [], [], fs, IpwrdB_C2Target, bwInband, bwACLROffsetC2, fnum, NameCelltitle, fnum_dir, 'Digital');
        
        for idC=1:NCarriers
            [IpwrdB_C1, ~, ~] = Pwr_Inband_g(fft(waveformCell{idC,:}, length(waveformCell{idC,:}), DIMFFT), fs(idC), bwInband(idC,:), [], 'full', 0);
            [IpwrdB_C2, ~, ~] = Pwr_Inband_g(fft(waveformC2cell_wiGain{idC,:}, length(waveformC2cell_wiGain{idC,:}), DIMFFT), fs(idC), bwInband(idC,:), [], 'full', 0);
            IpwrdB_C1_mean(idC,:)=round(mean(IpwrdB_C1),2);
            IpwrdB_C2_mean(idC,:)=round(mean(IpwrdB_C2),2);
        end
        
        % export
        waveformCell = waveformC2cell_wiGain;
        
        % excel
        ExcelWrite(round(EVMc2_Gain,2),'EVMc2_Gain',size(EVMc2_Gain),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(IpwrdB_C1_mean,'IpwrdB_C1_mean',size(IpwrdB_C1_mean),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(IpwrdB_C2_mean,'IpwrdB_C2_mean',size(IpwrdB_C2_mean),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        if flagC2_bwACLROffset
            ACLRdBLeft_C2Gain_min = fix(tableC2_Gain.tableACLR.ACLRdBLeft_min);
            ACLRdBRight_C2Gain_min = fix(tableC2_Gain.tableACLR.ACLRdBRight_min);
            ExcelWrite(ACLRdBLeft_C2Gain_min,'ACLRdBLeft_C2Gain_min',size(ACLRdBLeft_C2Gain_min),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
            ExcelWrite(ACLRdBRight_C2Gain_min,'ACLRdBRight_C2Gain_min',size(ACLRdBRight_C2Gain_min),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        end
    end
    
    %% C3. Export to waveformCA
    waveformCell_CA = waveformCell;
    fsCA = fs;
    NsampsCA = Nsamps;
    bwInbandCA = bwInband;
    clear NsampsBB NsampsTX NsampsRX NsampsDP
    
end

%% C4. 2021-04-22, TEST, flagC4_TEST_Delay, Apply Phase Delay to waveform
flagC4_TEST_Delay = 0
if flagC4_TEST_Delay
    npointOfDelay = 1104.3
    n = 1
    waveformCA = waveformCell{n};
    waveformCA_Delay = dsp_delay_cs_g(waveformCA, fs, npointOfDelay/fs);
    fnum_debug = 20210423
    PLOT_FFT_dB_g(waveformCA, fs, Nsamps, ['Original'], 'dt', 'full', 'amp', [fnum_debug]);
    PLOT_FFT_dB_g(waveformCA_Delay, fs, Nsamps, ['Delay'], 'dt', 'full', 'amp', [fnum_debug]);
    opts.fnum = fnum_debug*10
    opts.fs = fs
    opts.fs = []
    [mea_waveformCA_Delay_Cor, ref_waveformCA_Delay, DelayTotal] = dsp_Delay_EstCor_g(waveformCA, waveformCA_Delay, opts);
    [evmC4_RefCor, delayC4a] = dsp_evm_timexcorr_inband_g(ref_waveformCA_Delay, waveformCA_Delay, fs, bwInband, [], []);
    [evmC4_MeaCor, delayC4b] = dsp_evm_timexcorr_inband_g(waveformCA, mea_waveformCA_Delay_Cor, fs, bwInband, [], []);
    
    [evmC4c, delayC4c] = dsp_evm_timexcorr_inband_g(waveformCA, waveformCA_Delay, fs, bwInband, 'on', []);
    
end

%% C5. 2021-05-02, TEST, flagC5_TEST_UpsFilterConv, Estimate the delay contributed by FIR CONV method
flagC5_TEST_UpsFilterConv = 0
if flagC5_TEST_UpsFilterConv
    n = 1
    waveformC5 = waveformCell_CA{n};
    ratioC5_UPS = 4
    fnum_debug = 20210502
    
    % filter parameters
    FIRc5_Ftype = 'LPF'
    FIRc5_Wtype = 'Remez'
    FIRc5_K_AttdB = [1, 20]
    FIRc5_Wtype = 'Kaiser'
    FIRc5_K_AttdB = [50]
    FIRc5_K_fdelta = 1e6;
    FIRc5_fTolerance = 5e6
    FIRc5_Order = []
    FIRc5_fcutoffL = []
    FIRc5_fcutoffH = []
    bwInbandC5 = bwInbandCA
    fsC5_Ups = ratioC5_UPS*fs
    NsampsC5 = ratioC5_UPS*Nsamps
    [FIRc5_CHF, b_c5Test] = SYM_FIRApp(FIRc5_Wtype,FIRc5_Ftype,FIRc5_Order,FIRc5_K_AttdB,FIRc5_K_fdelta,FIRc5_fTolerance,...
        FIRc5_fcutoffL,FIRc5_fcutoffH,df,bwInbandC5,NCarriers,fsC5_Ups,fnum_debug);
    
    % FIR TimeDomain Convolution+Interploation
    waveformC5_Ups = zeros(ratioC5_UPS, Nsamps);
    waveformC5_Ups(1,:) = waveformC5;
    waveformC5_Ups = waveformC5_Ups(:);
    
    waveformC5_Ups_tdFIR = DSP_filter_g(FIRc5_CHF{:}, waveformC5_Ups, 'td');
    PLOT_FFT_dB_g(waveformC5_Ups_tdFIR, fsC5_Ups, NsampsC5, ['waveformC5 Ups tdFIR'], 'df', 'full', 'pwr', [fnum_debug]);
    
    waveformC5_Ups_td2FIR = DSP_filter_g(FIRc5_CHF{:}, waveformC5_Ups, 'td2');
    PLOT_FFT_dB_g(waveformC5_Ups_td2FIR, fsC5_Ups, NsampsC5, ['waveformC5 Ups td2FIR'], 'df', 'full', 'pwr', [fnum_debug]);
    
    % FIR FreqDomain Convolution+Interploation
    waveformC5_Ups_FDFIR = DSP_filter_g(FIRc5_CHF{:}, waveformC5_Ups, 'FD');
    PLOT_FFT_dB_g(waveformC5_Ups_FDFIR, fsC5_Ups, NsampsC5, ['waveformC5 Ups FDFIR'], 'df', 'full', 'pwr', [fnum_debug]);
    
    % PolyPhase+Interploation
    waveformC5_PolyPhaseUps = FIR_PolyPhase_g(waveformC5, FIRc5_CHF{:}, ratioC5_UPS, 0);
    PLOT_FFT_dB_g(waveformC5_PolyPhaseUps, fsC5_Ups, NsampsC5, ['waveformC5 Ups PolyPhase'], 'df', 'full', 'pwr', [fnum_debug]);
    
    % Delay corect, waveformC5_PolyPhaseUps
    opts.fnum = fnum_debug*10
    opts.fs = fsC5_Ups
    opts.fs = []
    [waveformC5_PolyPhaseUps_MeaDelayCor, waveformC5_RefDelay, DelayTotal] = dsp_Delay_EstCor_g(waveformC5_Ups_FDFIR, waveformC5_PolyPhaseUps, opts);
    PLOT_FFT_dB_g(waveformC5_PolyPhaseUps_MeaDelayCor, fsC5_Ups, NsampsC5, ['waveformC5 PolyPhaseUps DelayCor'], 'df', 'full', 'pwr', [fnum_debug]);
    
    [evmC5_FDAndtD, delayC5a] = dsp_evm_timexcorr_inband_g(waveformC5_Ups_FDFIR, waveformC5_Ups_tdFIR, fsC5_Ups, bwInbandC5, 'on', 2)
    [evmC5_FDAndPP, delayC5b] = dsp_evm_timexcorr_inband_g(waveformC5_Ups_FDFIR, waveformC5_PolyPhaseUps, fsC5_Ups, bwInbandC5, 'on', 2)
    [evmC5_FDAndtD2, delayC5c] = dsp_evm_timexcorr_inband_g(waveformC5_Ups_FDFIR, waveformC5_Ups_td2FIR, fsC5_Ups, bwInbandC5, 'on', 2)
    [evmC5_FDAndPPDelayCor, delayC5d] = dsp_evm_timexcorr_inband_g(waveformC5_Ups_FDFIR, waveformC5_PolyPhaseUps_MeaDelayCor, fsC5_Ups, bwInbandC5, 'on', 2)
    
    %% C6. 2021-05-04, TEST, flagC6_TEST_DnsPolyPhaseFIR
    fnum_debug = 20210505
    flagC6_TEST_DnsPolyPhaseFIR = 0
    FIRc6_Ftype = 'LPF'
    FIRc6_Wtype = 'Remez'
    FIRc6_K_AttdB = [1, 20]
    FIRc6_Wtype = 'Kaiser'
    FIRc6_K_AttdB = [50]
    FIRc6_K_fdelta = 2e6;
    FIRc6_fTolerance = 2e6
    FIRc6_Order = []
    FIRc6_fcutoffL = []
    FIRc6_fcutoffH = []
    bwInbandC6 = bwInbandCA
    fsC6_Dns = 1*fsC5_Ups
    ratioC6_DNS = 1/ratioC5_UPS
    NsampsC6 = NsampsC5/ratioC5_UPS
    [FIRc6_CHF, b_c6Test] = SYM_FIRApp(FIRc6_Wtype,FIRc6_Ftype,FIRc6_Order,FIRc6_K_AttdB,FIRc6_K_fdelta,FIRc6_fTolerance,...
        FIRc6_fcutoffL,FIRc6_fcutoffH,df,bwInbandC6,NCarriers,fsC6_Dns,fnum_debug);
    
    % PolyPhase Decimation
    waveformC6_PolyPhaseDns = FIR_PolyPhase_g(waveformC5_PolyPhaseUps, FIRc6_CHF{:}, ratioC6_DNS, 1);
    PLOT_FFT_dB_g(waveformC6_PolyPhaseDns, fsC6_Dns, NsampsC6, ['waveformC6 Dns PolyPhase'], 'df', 'full', 'pwr', [fnum_debug]);
    
    % FIR FreqDomain Convolution + Decimation
    waveformC6_FDFIR = DSP_filter_g(FIRc6_CHF{:}, waveformC5_PolyPhaseUps, 'td');
    waveformC6_FDFIR_Dns = waveformC6_FDFIR(1:1/ratioC6_DNS:end);
    PLOT_FFT_dB_g(waveformC6_FDFIR_Dns, fsC6_Dns, NsampsC6, ['waveformC6 Dns FDFIR'], 'df', 'full', 'pwr', [fnum_debug]);
    PLOT_FFT_dB_g(waveformC6_PolyPhaseDns, fsC6_Dns, NsampsC6, ['waveformC6 Dns PolyPhase'], 'dt', 'full', 'amp', [fnum_debug*2]);
    PLOT_FFT_dB_g(waveformC6_FDFIR_Dns, fsC6_Dns, NsampsC6, ['waveformC6 Dns FDFIR'], 'dt', 'full', 'amp', [fnum_debug*2]);
    
    opts.fnum = fnum_debug*10
    opts.fs = fsC6_Dns
    opts.fs = []
    [waveformC6_PolyPhaseDns_MeaDelayCor, ref_waveformC6_RefDelay, DelayTotal] = dsp_Delay_EstCor_g(waveformC6_FDFIR_Dns, waveformC6_PolyPhaseDns, opts);
    PLOT_FFT_dB_g(waveformC6_PolyPhaseDns_MeaDelayCor, fsC6_Dns, NsampsC6, ['waveformC6 Dns PolyPhase DelayCor'], 'dt', 'full', 'amp', [fnum_debug*2]);
    
    [evmC6_FDAndPP, delayC6a] = dsp_evm_timexcorr_inband_g(waveformC6_FDFIR_Dns, waveformC6_PolyPhaseDns, fsC6_Dns, bwInbandC6, 'on', 2);
    [evmC6_FDAndPPDelayCor, delayC6b] = dsp_evm_timexcorr_inband_g(waveformC6_FDFIR_Dns, waveformC6_PolyPhaseDns_MeaDelayCor, fsC6_Dns, bwInbandC6, 'on', 2);
    
end

%% C6, 2021-05-14, TEST, flagC6_TEST_ACBrGainVariation: waveformC0_AC with Gain variations for each Branch
flagC6_TEST_ACBrGainVariation = 0;
if flagC6_TEST_ACBrGainVariation
    fnum_debug = 20210514
    % Gain vs Branch
    PGaindB_C6 = 40*[0, 1, -1, 0, 1, 0, -1, 0];
    PGaindB_C6 = 60*[0, 1, -1, 0, 1, 0, -1, 0];
    PGaindB_C6 = 1*[3, -3, -1, 2, 1, 10, -7, 4];
    
    % Add Gain variations
    waveformC6_AC_GainVsBr = waveformC0_AC.*10.^(PGaindB_C6/20);
    [IpwrdB_C6_ACaddGain, ~, ~] = Pwr_Inband_g(fft(waveformC6_AC_GainVsBr), fs, [], 0, 'full', []);
    for idBR=1:numel(IpwrdB_waveformD7c_AC_GainVsBr)
        PLOT_FFT_dB_g(waveformC6_AC_GainVsBr(:,idBR), fsAC, NsampsC0_AC*2^4, ['wfAC +Gain Br',num2str(idBR),', pwrdB:',num2str(round(IpwrdB_C6_ACaddGain(idBR),2))], 'df', 'full', 'pwr', fnum_debug+1);
    end
    % Add Phase shifts
    ACPhsShiftNbrDeg_T7 = ExcelRead('ACPhsShiftNbrDeg_T7',4,'row','vector',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    waveformC6_AC_GainVsBr_PhShift = exp(1i*ACPhsShiftNbrDeg_T7./180*pi).*waveformC6_AC_GainVsBr;
    
    % Add Phase shifts2 
    for idBR = 1:numel(IpwrdB_waveformD7c_AC_GainVsBr)
        waveformC6_AC_GainVsBr_PhShift(:,idBR) = circshift(waveformC6_AC_GainVsBr(:,idBR),ACPhsShiftNbrDeg_T7(idBR),1);
    end
    
    % Sum of AC waveform
    waveformC6_AC_GainVsBr_PhShift_Sum = sum(waveformC6_AC_GainVsBr_PhShift,2);
    PLOT_FFT_dB_g(waveformC6_AC_GainVsBr_PhShift_Sum, fsAC, NsampsC0_AC*2^4, 'wfAC +PhShift +SumAllBr', 'df', 'full', 'pwr', fnum_debug+2);
    [IpwrdB_C6_ACSum, ~, ~] = Pwr_Inband_g(fft(waveformC6_AC_GainVsBr_PhShift_Sum), fs, [], 0, 'full', []);
    
    % AC DeModulation
    ACDomod_Nshift = 0;
    [ACdmod_t0Mean_C6, ACdmod_p0DegMean_C6, ACdmod_phEstDeg_C6, ACdmod_phEstDegDrift_C6, ACdmod_dataCapCor_C6, ACdmod_SNR_C6, ACdmod_dataCapNbrwoPD_C6] = ...
        AntCal_phaseDemodulateApp_g100(waveformC6_AC_GainVsBr_PhShift_Sum,ACconfig,ACDomod_Nshift,fnum_debug+3, fnum_dir);
    
%     % AC EQ
%     [waveformC6_ACPhsAlign,ACEQ_taps_C6,ACEQ_B_C6] = AntCal_EQphsAlignApp(waveformC0_AC,ACconfig,[], ACdmod_phEstDeg_C6,ACEQ_numFirTaps,ACEQ_NoOfAlignBr,fnum, fnum_dir);
    
end

%% Baseband Processing ==========================================================================================================================
if ~exist('NsampsBB','var')
    disp('Baseband Block ====================================================')
    %% B0. Carriers import from Baseband block
    waveformCell = waveformCell_CA;
    fs = fsCA;
    Nsamps = NsampsCA;
    bwInband = bwInbandCA;
    NCarriers = size(waveformCell_CA,1);
    NameCell = NameC0cell;
    
    
    %% B1. Apply Channel Filter
    flagB1_CHFilter = ExcelRead('flagB1_CHFilter',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagB1_CHFilter,'off')
        % Input
        FIRb1_Wtype = ExcelRead('FIRb1_Wtype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRb1_Ftype = ExcelRead('FIRb1_Ftype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRb1_Order = ExcelRead('FIRb1_Order',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRb1_fTolerance = ExcelRead('FIRb1_fTolerance',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRb1_fTolerance = 0.8e6
        FIRb1_K_AttdB = ExcelRead('FIRb1_K_AttdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRb1_K_fdelta = ExcelRead('FIRb1_K_fdelta',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRb1_fcutoffL = ExcelRead('FIRb1_fcutoffL',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRb1_fcutoffH = ExcelRead('FIRb1_fcutoffH',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        
        FIRb1_Wtype = 'Remez'
        FIRb1_K_AttdB = [1, 20]
        FIRb1_Wtype = 'Kaiser'
        FIRb1_K_AttdB = [50]
        FIRb1_K_fdelta = 1e6;
        
        [FIRb1_CHF, b_B1CHF] = SYM_FIRApp(FIRb1_Wtype,FIRb1_Ftype,FIRb1_Order,FIRb1_K_AttdB,FIRb1_K_fdelta,FIRb1_fTolerance,FIRb1_fcutoffL,FIRb1_fcutoffH,df,bwInband,NCarriers);
        
        %         FIRb1_Wtype = 'Remez';
        %         FIRb1_K_AttdB = [2 50];
        %         [FIRb1_CHF, b_B1CHF] = SYM_FIRApp(FIRb1_Wtype,FIRb1_Ftype,FIRb1_Order,FIRb1_K_AttdB,FIRb1_K_fdelta,FIRb1_fTolerance,FIRb1_fcutoffL,FIRb1_fcutoffH,df,bwInband,NCarriers,fs,[]);
        
        fnum = fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagB1_CHFilter},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        [waveformB1cell_CHF, b_B1CHF, EVMb1_CHF, tableB1_CHF] = SYM_FilterApp(waveformCell, FIRb1_CHF, fs, mehtod_filter, bwInband, fnum, NameCelltitle, fnum_dir);
        
        % export
        waveformCell= waveformB1cell_CHF;
        
        % excel
        ExcelWrite(tableB1_CHF.tableInput.FIR_Order,'FIRb1_Order_out',size(tableB1_CHF.tableInput.FIR_Order),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(EVMb1_CHF,2),'EVMb1_CHF',size(EVMb1_CHF),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        
    end
    
    %% B2. Add RippleEQ
    flagB2_addEQRipple = ExcelRead('flagB2_addEQRipple',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagB2_addEQRipple,'off')
        % input
        RippledBMAX_B2 = ExcelRead('RippledBMAX_B2',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        
        fnum=fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagB2_addEQRipple},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        [waveformB2cell_EQRipple, EVMb2_EQRipple, tableB2_EQRipple] = SYM_MagRippleEQApp(waveformCell, RippledBMAX_B2, fs, bwInband, bwInband, fnum, NameCelltitle, fnum_dir);
        
        % export waveform
        waveformCell = waveformB2cell_EQRipple;
        
        % excel
        ExcelWrite(round(EVMb2_EQRipple,2),'EVMb2_EQRipple',size(EVMb2_EQRipple),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    %% B2b. 2021-04-27, TEST, OFDM DeModulation to evaluate the Ripple effect
    flagB2b_DemodForRipple = 0
    
    if flagB2b_DemodForRipple
        waveformB2Bcell_ref = waveformCell_CA
        [rxGridB2Bcell, EVMB2B_DemodulationDL] = SYM_DeMODApp(waveformB2cell_EQRipple, waveformB2Bcell_ref, DLconfig, fs, bwInbandC0, fnum, fnum_dir);
    end
    
    
    %% B3. Add Digital Gain
    flagB3_addDigGain = ExcelRead('flagB3_addDigGain',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagB3_addDigGain,'off')
        % input
        GaindB_B3 = ExcelRead('GaindB_B3',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        %         IpwrdB_B3Target = ExcelRead('IpwrdB_B3Target',[NCarriers],'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        IpwrdB_B3Target = ExcelRead('IpwrdB_B3Target',[NCarriers],'column','cell',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        bwACLROffsetB3 = ExcelRead('bwACLROffsetB3',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        %         bwACLROffsetB3 = sum(bwChNC);
        if fs<2*bwACLROffsetB3
            bwACLROffsetB3 = [];
            flagB3_bwACLROffset = 0;
        else
            flagB3_bwACLROffset = 1;
        end
        
        fnum=fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagB3_addDigGain},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        [waveformB3cell_wiGain, EVMb3_Gain,  tableB3_Gain] = SYM_AMPApp(waveformCell, GaindB_B3, [], [], [], [], [], [], fs, IpwrdB_B3Target, bwInband, bwACLROffsetB3, fnum, NameCelltitle, fnum_dir, 'Digital');
        
        for idC=1:NCarriers
            [IpwrdB_B2, ~, ~] = Pwr_Inband_g(fft(waveformCell{idC,:}, length(waveformCell{idC,:}), DIMFFT), fs(idC,:), bwInband(idC,:), [], 'full', 0);
            [IpwrdB_B3, ~, ~] = Pwr_Inband_g(fft(waveformB3cell_wiGain{idC,:}, length(waveformB3cell_wiGain{idC,:}), DIMFFT), fs(idC,:), bwInband(idC,:), [], 'full', 0);
            IpwrdB_B2_mean(idC,:)=round(mean(IpwrdB_B2),2);
            IpwrdB_B3_mean(idC,:)=round(mean(IpwrdB_B3),2);
        end
        
        % export
        waveformCell = waveformB3cell_wiGain;
        
        % excel
        ExcelWrite(round(EVMb3_Gain,2),'EVMb3_Gain',size(EVMb3_Gain),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(IpwrdB_B2_mean,'IpwrdB_B2_mean',size(IpwrdB_B2_mean),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(IpwrdB_B3_mean,'IpwrdB_B3_mean',size(IpwrdB_B3_mean),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        
        if flagB3_bwACLROffset
            ACLRdBLeft_B3Gain_min = fix(tableB3_Gain.tableACLR.ACLRdBLeft_min);
            ACLRdBRight_B3Gain_min = fix(tableB3_Gain.tableACLR.ACLRdBRight_min);
            ExcelWrite(ACLRdBLeft_B3Gain_min,'ACLRdBLeft_B3Gain_min',size(ACLRdBLeft_B3Gain_min),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
            ExcelWrite(ACLRdBRight_B3Gain_min,'ACLRdBRight_B3Gain_min',size(ACLRdBRight_B3Gain_min),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        end
    end
    
    %% B4. DUC stage B4
    flagB4_DUC = ExcelRead('flagB4_DUC',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagB4_DUC,'off')
        ratioB4_DUC = ExcelRead('ratioB4_DUC',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        fsB4_DUC = ratioB4_DUC*fs;
        bwChCA_2times = 2*sum(diff(bwInbandCA,[],2),1);
        if fsB4_DUC/2 < bwChCA_2times
            ratioB4_DUC = ceil(bwChCA_2times/max(fsCA/2))
            fsB4_DUC = ratioB4_DUC*fs;
            ExcelWrite(ratioB4_DUC,'ratioB4_DUC',size(ratioB4_DUC),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        end
        
        % excel
        %         ExcelWrite(round(fsB4_DUC/1e6,2),'fsB4_DUC_MHz',size(fsB4_DUC),'r1col',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(fsB4_DUC/1e6,2),'fsB4_DUC_MHz',size(fsB4_DUC),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        
        %% B4. FIRb4DUC
        FIRb4_Wtype = ExcelRead('FIRb4_Wtype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRb4_Ftype = ExcelRead('FIRb4_Ftype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRb4_Order = ExcelRead('FIRb4_Order',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRb4_fTolerance = ExcelRead('FIRb4_fTolerance',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRb4_K_AttdB = ExcelRead('FIRb4_K_AttdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRb4_K_fdelta = ExcelRead('FIRb4_K_fdelta',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRb4_fcutoffL = ExcelRead('FIRb4_fcutoffL',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRb4_fcutoffH = ExcelRead('FIRb4_fcutoffH',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        
        FIRb4_Wtype = 'Remez'
        FIRb4_Ftype = 'LPF'
        FIRb4_fTolerance = 2e6
        FIRb4_K_AttdB = [0.1, 50]
        FIRb4_K_fdelta = 2e6
        FIRb4_fcutoffL = max(bwInband,[],2)
        FIRb4_DUC = SYM_FIRApp(FIRb4_Wtype,FIRb4_Ftype,FIRb4_Order,FIRb4_K_AttdB,FIRb4_K_fdelta,FIRb4_fTolerance,FIRb4_fcutoffL,FIRb4_fcutoffH,df,bwInband,NCarriers,fs*ratioB4_DUC,[]);
        if strcmpi(FIRb4_Wtype,'PM_lpf')
            devs_B4 = [(10^(FIRb4_K_AttdB(1)/20)-1)/(10^(FIRb4_K_AttdB(1)/20)+1) 10^(-FIRb4_K_AttdB(2)/20)];
            for k=1:size(FIRb4_fcutoffL,1)
                [n,fo,ao,w] = firpmord([FIRb4_fcutoffL(k,:) FIRb4_fcutoffL(k,:)+FIRb4_K_fdelta], [1 0], [devs_B4], fs(k,:)*ratioB4_DUC);
                FIRb4_DUC{k,:} = firpm(n,fo,ao,w);
                PLOT_FFT_dB_g(cell2mat(FIRb4_DUC(k,:))*Nsamps, fs(k,:)*ratioB4_DUC, Nsamps, ['Length:',num2str(length(cell2mat(FIRb4_DUC(k,:))))], 'df', 'full', 'pwr', [20210419], [], []);
                
            end
        end
        
        fnum=fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagB4_DUC},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        [waveformB4cell_DUC, fsB4_DUC, NsampsB4_DUC, EVMb4_DUC, tableB4_DUC] = SYM_UPSampApp(waveformCell, ratioB4_DUC, fs, [FIRb4_DUC], mehtod_filter, bwInbandC0, [], [], fnum, NameCelltitle, fnum_dir);
        
        % export
        waveformCell= waveformB4cell_DUC;
        fs=fsB4_DUC;
        Nsamps=NsampsB4_DUC;
        
        % excel
        ExcelWrite(tableB4_DUC.tableFilterInput.FIR_Order,'FIRb4_Order_out',size(tableB4_DUC.tableFilterInput.FIR_Order),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(EVMb4_DUC,2),'EVMb4_DUC',size(EVMb4_DUC),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(fs/1e6),'fsB4_DUC_MHz',size(fs),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(Nsamps,'NsampsB4_DUC',size(Nsamps),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    %% B5. Add Channel Offset
    flagB5_addChOffset = ExcelRead('flagB5_addChOffset',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    offsetCH_B5 = ExcelRead('offsetCH_B5',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagB5_addChOffset,'off') % all of offset_CH are ZERO
        tableB5Input_offsetCH = table(offsetCH_B5)
        %     for idBR=1:size(offsetCH_B5,mod(DIMFFT,2)+1)
        for idC=1:NCarriers
            waveformB5cell_offsetCH{idC,:} = circshift(waveformCell{idC},offsetCH_B5,DIMFFT);
        end
        
        % export
        waveformCell = waveformB5cell_offsetCH;
    end
    
    %% B6. NCO shift and Combination
    flagB6_addNCO = ExcelRead('flagB6_addNCO',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if flagB6_addNCO>0
        fNCO = ExcelRead('fNCO',[flagB6_addNCO],'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        flagB6_addNCO_MultiCarrier = ExcelRead('flagB6_addNCO_MultiCarriers',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        NCarriersB6input = size(waveformCell,1);
        
        if fs/2 < (max(abs(fNCO))+abs(diff(bwInband)))
            %% 2020-05-13, Sampling Rate should larger than IBW before NCO and MultiCarrier Combination
            error('the Sampling Rate should larger than IBW before NCO and MultiCarrier Combination!')
        end
        fnum=fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{'B6'},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        idC=1;
        [waveformB6cell_NCO, waveformB6cell_NCOComb, bwInbandB6, EVMb6_NCOComb, tableB6_NCO] = SYM_NCOApp(waveformCell(idC:end,:), fix(fNCO(idC:end,:)/df)*df, fs, bwInband(idC:end,:), flagB6_addNCO_MultiCarrier, fnum, NameCelltitle(idC:end,:), fnum_dir);
        
        % export
        if exist('waveformB6cell_NCOComb','var')&&~isempty(waveformB6cell_NCOComb)
            waveformCell = waveformB6cell_NCOComb;
            bwInbandB6_NCOComb = [min(bwInbandB6(:)) max(bwInbandB6(:))];
            bwInband = bwInbandB6_NCOComb;
            NameC6cell = {['CarriersSUM']};
            if all(diff(df)==0)
                df = df(1);
            end
        else
            waveformCell = waveformB6cell_NCO;
            for idC=1:NCarriers
                [IpwrdB_B6(idC,:), ~, ~] = Pwr_Inband_g(fft(waveformB6cell_NCO{idC,:}, length(waveformB6cell_NCO{idC,:}), DIMFFT), fs, bwInbandB6(idC,:), [], 'full', 0);
            end
            bwInband=bwInbandB6;
            bwInbandB6_NCOComb = NaN;
            EVMb6_NCOComb = NaN;
            NameC6cell = NameCell;
            
        end
        % export
        NCarriersB6_NCOMultiCarrier=size(waveformCell,1);
        NCarriers = NCarriersB6_NCOMultiCarrier;
        NameCell = NameC6cell;
        
        % excel
        ExcelWrite(bwInbandB6,'bwInbandB6',size(bwInbandB6),'r1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(EVMb6_NCOComb,2),'EVMb6_NCO',size(EVMb6_NCOComb,1),'r1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(NameCell,'CarriersB6_NCO',size(NameCell),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(NCarriers,'NCarriersB6_NCO',size(NCarriers),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(bwInbandB6_NCOComb,'bwInbandB6_NCOComb',size(bwInbandB6_NCOComb,1),'r1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        
    end
    
    %% B7. CCDF and Clipping
    flagB7_Clipping = ExcelRead('flagB7_Clipping',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    %% 2020-04-23, Issue6, by flag_ACLR could imporve the NOISE FLOOR and ACLR, but EVM get worse !!!
    if ~strcmp(flagB7_Clipping,'off')
        % input
        thresholddB_B7 = ExcelRead('thresholddB_B7',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        bwACLROffsetB7 = sum(bwChNC);
        flagB7_CFR_ACLRLimiter = ExcelRead('flagB7_CFR_ACLRLimiter',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        flagB7_PEX_DmcItpRatio = ExcelRead('flagB7_PEX_DmcItp',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        
        % input
        bwACLROffsetB7 = sum(bwChNC);
        %         bwACLROffsetB7 = diff(bwInband,[],DIMFFT+1)
        
        if strcmpi(flagB7_CFR_ACLRLimiter, 'aclr')
            ACLRdBLimit_B7 = ExcelRead('ACLRdBLimit_B7',1,'column','cell',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
            FIRb7_CFR_PEX_Ch = ACLRdBLimit_B7;
            
        elseif strcmpi(flagB7_CFR_ACLRLimiter, 'PEXFIR')
            if ~strcmp(flagB7_PEX_DmcItpRatio, 'off') || isnumeric(flagB7_PEX_DmcItpRatio)
                FIRb7_CFR_PEX_DmcItp_Wtype = ExcelRead('FIRb7_CFR_PEX_DmcItp_Wtype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
                FIRb7_CFR_PEX_DmcItp_Wtype = 'Remez';
                FIRb7_CFR_PEX_DmcItp_Wtype = 'Kaiser';
                FIRb7_CFR_PEX_DmcItp_Ftype = 'LPF';
                FIRb7_CFR_PEX_DmcItp_Order = ExcelRead('FIRb7_CFR_PEX_DmcItp_Order',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
                FIRb7_CFR_PEX_DmcItp_fTolerance = ExcelRead('FIRb7_CFR_PEX_DmcItp_fTolerance',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
                FIRb7_CFR_PEX_DmcItp_K_AttdB = ExcelRead('FIRb7_CFR_PEX_DmcItp_K_AttdB',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
                FIRb7_CFR_PEX_DmcItp_K_AttdB = [1 40];
                FIRb7_CFR_PEX_DmcItp_K_AttdB = [20];
                FIRb7_CFR_PEX_DmcItp_K_fdelta = ExcelRead('FIRb7_CFR_PEX_DmcItp_K_fdelta',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
                FIRb7_CFR_PEX_DmcItp_fcutoffL = ExcelRead('FIRb7_CFR_PEX_DmcItp_fcutoffL',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
                FIRb7_CFR_PEX_DmcItp_fcutoffH = ExcelRead('FIRb7_CFR_PEX_DmcItp_fcutoffH',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
                
                PEX_DmcItpRatio = fs./fsC0;
                PEX_NCO = -fNCO;
                PEX_bwInband = bwInbandC0;
                PEX_NCarriers = NCarriersB6input;
                PEX_DmcIpt_fs = fs;
                PEX_Ch_fs = fsC0;
                
                fnum_FIRb7 = 20210417
                % generate FIRb7_CFR_PEX_DmcIpt
                [FIRb7_CFR_PEX_DmcIpt, ~, tableB7_PEXDmcIntFIR] = SYM_FIRApp(FIRb7_CFR_PEX_DmcItp_Wtype,FIRb7_CFR_PEX_DmcItp_Ftype,...
                    FIRb7_CFR_PEX_DmcItp_Order,FIRb7_CFR_PEX_DmcItp_K_AttdB,FIRb7_CFR_PEX_DmcItp_K_fdelta,FIRb7_CFR_PEX_DmcItp_fTolerance,...
                    FIRb7_CFR_PEX_DmcItp_fcutoffL,FIRb7_CFR_PEX_DmcItp_fcutoffH,df,PEX_bwInband,PEX_NCarriers,PEX_DmcIpt_fs,fnum_FIRb7);
                
            else
                PEX_DmcItpRatio = [];
                PEX_NCO = [];
                PEX_bwInband = [];
                PEX_bwInband = bwInband;
                PEX_NCarriers = NCarriersB6_NCOMultiCarrier;
                PEX_DmcIpt_fs = [];
                PEX_Ch_fs = fs;
                
                FIRb7_CFR_PEX_DmcIpt = [];
            end
            
            FIRb7_CFR_PEX_Ch_Wtype = ExcelRead('FIRb7_CFR_PEX_Ch_Wtype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
            FIRb7_CFR_PEX_Ch_Ftype = 'LPF';
            FIRb7_CFR_PEX_Ch_Order = ExcelRead('FIRb7_CFR_PEX_Ch_Order',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
            FIRb7_CFR_PEX_Ch_fTolerance = ExcelRead('FIRb7_CFR_PEX_Ch_fTolerance',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
            FIRb7_CFR_PEX_Ch_K_AttdB = ExcelRead('FIRb7_CFR_PEX_Ch_K_AttdB',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
            FIRb7_CFR_PEX_Ch_K_fdelta = ExcelRead('FIRb7_CFR_PEX_Ch_K_fdelta',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
            FIRb7_CFR_PEX_Ch_fcutoffL = ExcelRead('FIRb7_CFR_PEX_Ch_fcutoffL',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
            FIRb7_CFR_PEX_Ch_fcutoffH = ExcelRead('FIRb7_CFR_PEX_Ch_fcutoffH',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
            
            if strcmp(FIRb7_CFR_PEX_Ch_Wtype,'Remez')
                FIRb7_CFR_PEX_Ch_Wtype = 'Remez';
                FIRb7_CFR_PEX_Ch_K_AttdB = [1 20];
                tableB7_FIR_PEX_Ch = table(FIRb7_CFR_PEX_Ch_Wtype, FIRb7_CFR_PEX_Ch_Order, FIRb7_CFR_PEX_Ch_fTolerance)
                
            elseif strcmp(FIRb7_CFR_PEX_Ch_Wtype,'Kaiser')
                FIRb7_CFR_PEX_Ch_Wtype = 'Kaiser';
                FIRb7_CFR_PEX_Ch_K_AttdB = [10];
            end
            
            % generate FIRb7_CFR_PEX_Ch
            [FIRb7_CFR_PEX_Ch, ~, tableB7_PEXChFIR] = SYM_FIRApp(FIRb7_CFR_PEX_Ch_Wtype,FIRb7_CFR_PEX_Ch_Ftype,...
                FIRb7_CFR_PEX_Ch_Order,FIRb7_CFR_PEX_Ch_K_AttdB,FIRb7_CFR_PEX_Ch_K_fdelta,FIRb7_CFR_PEX_Ch_fTolerance,...
                FIRb7_CFR_PEX_Ch_fcutoffL,FIRb7_CFR_PEX_Ch_fcutoffH,df,PEX_bwInband,PEX_NCarriers,PEX_Ch_fs,fnum_FIRb7);
            
        end
        
        fnum=fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagB7_Clipping},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag];
        
        % CFR App
        [waveformB7cell_CFRFIR, PARdB_B7CFR, EVMb7_CFR, tableB7_CFR] = SYM_CFRApp(waveformCell, thresholddB_B7, bwInband, fs, FIRb7_CFR_PEX_Ch, bwInband, bwACLROffsetB7, fnum, NameCelltitle, fnum_dir, ...
            PEX_DmcItpRatio, PEX_NCO, FIRb7_CFR_PEX_DmcIpt);
        
        % export
        waveformCell=waveformB7cell_CFRFIR;
        ACLRdBL_B7CFR_min=fix(tableB7_CFR.tableACLR.ACLRdBLeft_min);
        ACLRdBR_B7CFRFIR_min=fix(tableB7_CFR.tableACLR.ACLRdBRight_min);
        
        % excel
        ExcelWrite(round(PARdB_B7CFR,2),'PARdB_B7CFR',1,'r1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(tableB7_PEXDmcIntFIR.FIR_Output_Order,'FIRb7_CFR_PEX_DmcItp_Order_out',size(tableB7_PEXDmcIntFIR.FIR_Output_Order),'r1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(EVMb7_CFR,2),'EVMb7_CFR',1,'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(ACLRdBL_B7CFR_min,'ACLRdBLeft_B7CFR_min',1,'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(ACLRdBR_B7CFRFIR_min,'ACLRdBRight_B7CFR_min',1,'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        
    end
    
    %% B8. DUC stage B8
    flagB8_DUC = ExcelRead('flagB8_DUC',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagB8_DUC,'off')
        % input
        ratioB8_DUC = ExcelRead('ratioB8_DUC',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        fsB8_DUC = ratioB8_DUC*fs;
        % excel
        ExcelWrite(fsB8_DUC,'fsB8_DUC',size(fsB8_DUC),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        
        %% B8. FIRb8_DUC
        FIRb8_Wtype = ExcelRead('FIRb8_Wtype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRb8_Ftype = ExcelRead('FIRb8_Ftype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRb8_Order = ExcelRead('FIRb8_Order',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRb8_fTolerance = ExcelRead('FIRb8_fTolerance',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRb8_K_AttdB = ExcelRead('FIRb8_K_AttdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRb8_K_fdelta = ExcelRead('FIRb8_K_fdelta',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRb8_fcutoffL = ExcelRead('FIRb8_fcutoffL',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRb8_fcutoffL = diff(bwInband)/2+max(bwInband);
        FIRb8_fcutoffH = ExcelRead('FIRb8_fcutoffH',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        
        fnum_Debug = 0419
        FIRb8_DUC = SYM_FIRApp(FIRb8_Wtype,FIRb8_Ftype,FIRb8_Order,FIRb8_K_AttdB,FIRb8_K_fdelta,FIRb8_fTolerance,FIRb8_fcutoffL,FIRb8_fcutoffH,df,bwInband,NCarriers,fs,fnum_Debug);
        
        fnum=fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagB8_DUC},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        [waveformB8cell_DUC, fsB8_DUC, NsampsB8_DUC, EVMb8_DUC, tableB8_DUC] = SYM_UPSampApp(waveformCell, ratioB8_DUC, fs, [FIRb8_DUC], mehtod_filter, bwInband, [], [], fnum, NameCelltitle, fnum_dir);
        
        % export
        waveformCell = waveformB8cell_DUC;
        fs=fsB8_DUC;
        Nsamps=NsampsB8_DUC;
        
        % excel
        ExcelWrite(tableB8_DUC.tableFilterInput.FIR_Order,'FIRb8_Order_out',size(tableB8_DUC.tableFilterInput.FIR_Order),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(EVMb8_DUC,2),'EVMb8_DUC',size(EVMb8_DUC),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(fs/1e6,2),'fsB8_DUC_MHz',size(fs),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(Nsamps,'NsampsB8_DUC',size(Nsamps),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    %% B9. Add Digital Gain
    flagB9_addDigGain = ExcelRead('flagB9_addDigGain',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagB9_addDigGain,'off')
        % input
        GaindB_B9 = ExcelRead('GaindB_B9',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        %% 2020-05-19, Compare results of different Target power
        IpwrdB_B9Target = ExcelRead('IpwrdB_B9Target',[1],'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        %         bwACLROffsetB9 = ExcelRead('bwACLROffsetB9',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        bwACLROffsetB9 = sum(bwChNC);
        
        
        fnum=fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagB9_addDigGain},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        [waveformB9cell_wiGain, EVMb9_Gain,  tableB9_Gain] = SYM_AMPApp(waveformCell, GaindB_B9, [], [], [], [], [], [], fs, IpwrdB_B9Target, bwInband, bwACLROffsetB9, fnum, NameCelltitle, fnum_dir, 'Digital');
        for idC=1:NCarriers
            [IpwrdB_B8, ~, ~] = Pwr_Inband_g(fft(waveformCell{idC,:}, length(waveformCell{idC,:}), DIMFFT), fs, bwInband, [], 'full', 0);
            [IpwrdB_B9, ~, ~] = Pwr_Inband_g(fft(waveformB9cell_wiGain{idC,:}, length(waveformB9cell_wiGain{idC,:}), DIMFFT), fs, bwInband, [], 'full', 0);
            IpwrdB_B8_mean(idC,:)=round(mean(IpwrdB_B8),2);
            IpwrdB_B9_mean(idC,:)=round(mean(IpwrdB_B9),2);
        end
        
        % export
        waveformCell = waveformB9cell_wiGain;
        ACLRdBLeft_B9Gain_min = fix(tableB9_Gain.tableACLR.ACLRdBLeft_min);
        ACLRdBRight_B9Gain_min = fix(tableB9_Gain.tableACLR.ACLRdBRight_min);
        
        % excel
        ExcelWrite(round(EVMb9_Gain,2),'EVMb9_Gain',size(EVMb9_Gain),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(ACLRdBLeft_B9Gain_min,'ACLRdBLeft_B9Gain_min',size(ACLRdBLeft_B9Gain_min),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(ACLRdBRight_B9Gain_min,'ACLRdBRight_B9Gain_min',size(ACLRdBRight_B9Gain_min),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(IpwrdB_B8_mean,'IpwrdB_B8_mean',size(IpwrdB_B8_mean),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(IpwrdB_B9_mean,'IpwrdB_B9_mean',size(IpwrdB_B9_mean),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    %% B10. Export to waveformBB
    waveformCell_BB = waveformCell;
    fsBB = fs;
    NsampsBB = Nsamps;
    bwInbandBB = bwInband;
    clear NsampsTX NsampsRX NsampsDP
    NameCell = NameC6cell;
end

%% Trasmitter Block =============================================================================================================================
if ~exist('NsampsTX','var')
    disp('Trasmitter Block ====================================================')
    %% T0. Baseband Signal import to Transmitter Block
    waveformCell = waveformCell_BB;
    fs = fsBB;
    Nsamps = NsampsBB;
    bwInband = bwInbandBB;
    
    %% T1. DUC stage T1
    flagT1_DUC = ExcelRead('flagT1_DUC',[1],'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagT1_DUC,'off')
        % input
        ratioT1_DUC = ExcelRead('ratioT1_DUC',[1],'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        fsT1_DUC = ratioT1_DUC*fs;
        % excel
        ExcelWrite(round(fsT1_DUC,2),'fsT1_DUC',size(fsT1_DUC),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        
        %% T1b. FIRt1DUC
        %     FIRt1_DUC={b_THB1,b_THB2};
        FIRt1_DUC = ExcelRead('FIRt1_DUC',[1],'row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        
        fnum=fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagT1_DUC},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        [waveformT1cell_DUC, fsT1_DUC, NsampsT1_DUC, EVMt1_DUC, tableT1_DUC] = SYM_UPSampApp(waveformCell, ratioT1_DUC, fs, [FIRt1_DUC(:,:)], mehtod_filter, bwInband, [], [], fnum, NameCelltitle, fnum_dir);
        
        % export
        waveformCell = waveformT1cell_DUC;
        fs=fsT1_DUC;
        Nsamps=NsampsT1_DUC;
        
        % excel
        ExcelWrite(round(EVMt1_DUC,2),'EVMt1_DUC',size(EVMt1_DUC),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(fs/1e6,2),'fsT1_DUC_MHz',size(fs),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(Nsamps,'NsampsT1_DUC',size(Nsamps),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    %% T2 UpSampling stage T2, to expand spectrum
    flagT2_UPS = ExcelRead('flagT2_UPS',[1],'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagT2_UPS,'off')
        % input
        fsT2_UPS = fix(orderHD*2*mean(fRF)/fs)*fs;
        fsT2_UPS = fix(orderHD*1*mean(fRF)/fs)*fs;
        fsT2_UPS = fix(1*1*mean(fRF)/fs)*fs;
        fsT2_UPS = ExcelRead('fsT2_UPS',[1],'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        ratioT2_UPS = fsT2_UPS/fs;
        
        % excel
        ExcelWrite(round(ratioT2_UPS,2),'ratioT2_UPS',size(ratioT2_UPS),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        
        %% T2a. FIRt2UPS
        FIRt2_Wtype = ExcelRead('FIRt2_Wtype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRt2_Ftype = ExcelRead('FIRt2_Ftype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRt2_Order = ExcelRead('FIRt2_Order',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRt2_fTolerance = ExcelRead('FIRt2_fTolerance',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRt2_K_AttdB = ExcelRead('FIRt2_K_AttdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRt2_K_fdelta = ExcelRead('FIRt2_K_fdelta',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRt2_fcutoffL = ExcelRead('FIRt2_fcutoffL',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRt2_fcutoffH = ExcelRead('FIRt2_fcutoffH',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        
        FIRt2_UPS = SYM_FIRApp(FIRt2_Wtype,FIRt2_Ftype,FIRt2_Order,FIRt2_K_AttdB,FIRt2_K_fdelta,FIRt2_fTolerance,FIRt2_fcutoffL,FIRt2_fcutoffH,df,bwInband,NCarriers);
        
        fnum=fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagT2_UPS},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        [waveformT2cell_UPS, fsT2_UPS, NsampsT2_UPS, EVMt2_UPS, tableT2_UPS] = SYM_UPSampApp(waveformCell, ratioT2_UPS, fs, [FIRt2_UPS], mehtod_filter, bwInband, [], [], fnum, NameCelltitle, fnum_dir);
        
        % export
        waveformCell = waveformT2cell_UPS;
        fs=fsT2_UPS;
        Nsamps=NsampsT2_UPS;
        
        % excel
        ExcelWrite(tableT2_UPS.tableFilterInput.FIR_Order,'FIRt2_Order_out',size(tableT2_UPS.tableFilterInput.FIR_Order),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(EVMt2_UPS,2),'EVMt2_UPS',size(EVMt2_UPS),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(fs/1e6,2),'fsT2_UPS_MHz',size(fs),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(Nsamps,'NsampsT2_UPS',size(Nsamps),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    %% T2b. 2021-04-09, TEST, flagT2b_TEST_ADC, ADC output, Quantization Error/Dither/NoiseShape
    flagT2b_TEST_ADC = 'on';
    if ~strcmpi(flagT2b_TEST_ADC,'off')
        Vref = 1.4;
        nbits = 12;
        LSB = Vref/2^nbits
        flag_RemoveDC = 1
        
        waveformT2b_ADC = ADC_Quantizer_Dither(cell2mat(waveformCell(1)), LSB, [], flag_RemoveDC, fs);
        [EVMt2b_ADC, ~] = dsp_evm_timexcorr_inband_g(cell2mat(waveformCell(1)), waveformT2b_ADC, fs, bwInband, [], 2)
        
        waveformT2b_ADCin = cell2mat(waveformCell(1));
        PLOT_FFT_dB_g(waveformT2b_ADCin(:,1), fs, length(waveformT2b_ADCin), ['waveformT2b ADCin'], 'df', 'full', 'pwr', [fnum], [], []);
        PLOT_FFT_dB_g(waveformT2b_ADC(:,1), fs, length(waveformT2b_ADC), ['waveformT2b ADCout'], 'df', 'full', 'pwr', [fnum], [], []);
        
        waveformCell(1) = {waveformT2b_ADC};
        
    end
    
    %% T3. Add Inband Noise to waveform
    flagT3_addAWGN = ExcelRead('flagT3_addAWGN',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagT3_addAWGN,'off')&&~strcmp(flagT2b_TEST_ADC,'on')
        % input
        IpwrdB1Hz_NoiseFloor_T3 = ExcelRead('IpwrdB1Hz_NoiseFloor_T3',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        %     spec_SNRdB_InbandNoise = 50;
        
        fnum=fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagT3_addAWGN},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        [waveformT3cell_AWGN,SNRdBInband_T3,EVMt3_AWGN,tableT3_AWGN] = SYM_AWGNApp(waveformCell, IpwrdB1Hz_NoiseFloor_T3, bwInband, fs, bwInband, fnum, NameCelltitle, fnum_dir);
        
        % export
        waveformCell = waveformT3cell_AWGN;
        
        % excel
        ExcelWrite(round(SNRdBInband_T3,2),'SNRdBInband_T3',size(SNRdBInband_T3),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(EVMt3_AWGN,2),'EVMt3_AWGN',size(EVMt3_AWGN),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    %% 2021-03-15, Debug AWGN results for AC
    flagT3_DebugAWGN_forAC = 'off';
    if strcmp(flagT3_DebugAWGN_forAC,'on')
        waveformT3_CombACNbr = sum(waveformT3cell_AWGN{1,:},mod(DIMFFT,2)+1); % all branch combination
        Nbr = size(waveformT3_CombACNbr,mod(DIMFFT,2)+1); % combine to One branch
        
        %% 2021-03-15, Check AC Combination: Ipwr of carrier and noise floor
        idBR = 1
        fnum(1) = fnum(1)+1;
        PLOT_FFT_dB_g(waveformT2cell_UPS{:}(:,idBR), fs, Nsamps, ['DLAC No.',num2str(idBR),' br'], 'df', 'full', 'pwr', fnum);
        PLOT_FFT_dB_g(waveformT3cell_AWGN{:}(:,idBR), fs, Nsamps, ['DLAC No. ',num2str(idBR),' br + AWGN'], 'df', 'full', 'pwr', fnum);
        PLOT_FFT_dB_g(waveformT3_CombACNbr, fs, Nsamps, ['DLAC all ',num2str(size(waveformT3_CombACNbr,2)),' brs +AWGN combination'], 'df', 'full', 'pwr', fnum);
        [IpwrdB_T3_1BR, ~, ~] = Pwr_Inband_g(fft(waveformT3cell_AWGN{:}(:,idBR), Nsamps, DIMFFT), fs, bwInbandBB, [], 'full', []);
        [IpwrdB_T3_8BR, ~, ~] = Pwr_Inband_g(fft(waveformT3_CombACNbr, Nsamps, DIMFFT), fs, bwInbandBB, [], 'full', []);
        [IpwrdB_T2_1BR_NoiseFloor, ~, ~] = Pwr_Inband_g(fft(waveformT2cell_UPS{:}(:,idBR), Nsamps, DIMFFT), fs, bwInbandBB+50e6, [], 'full_psd', []);
        [IpwrdB_T3_1BR_NoiseFloor, ~, ~] = Pwr_Inband_g(fft(waveformT3cell_AWGN{:}(:,idBR), Nsamps, DIMFFT), fs, bwInbandBB+50e6, [], 'full_psd', []);
        [IpwrdB_T3_8BR_NoiseFloor, ~, ~] = Pwr_Inband_g(fft(waveformT3_CombACNbr, Nsamps, DIMFFT), fs, bwInbandBB+50e6, [], 'full_psd', []);
        delta_IpwrdB_CombNoiseFloor = IpwrdB_T3_8BR_NoiseFloor - IpwrdB_T3_1BR_NoiseFloor;
    end
    
    %% T4. LO Generator
    % excel import
    LOU_AMtoPMPhsDriftDeg = ExcelRead('LOU_AMtoPMPhsDriftDeg',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    LOU_leveldB = ExcelRead('LOU_leveldB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    LOU_fLO = ExcelRead('LOU_fLO',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    LOU_PN_ThetaDeg = ExcelRead('LOU_PN_ThetaDeg',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    LOU_PN_MagDriftdB1Hz = ExcelRead('LOU_PN_MagDriftdB1Hz',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    LOU_PN_offset = ExcelRead('LOU_PN_offset',2,'row','cell',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    LOU_PN_f_offset_Hz = LOU_PN_offset{1};
    LOU_PN_g_offset_dBc1Hz = LOU_PN_offset{2};
    
    LOU_IMB_PhsDeg = ExcelRead('LOU_IMB_PhsDeg',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    LOU_IMB_MagdB = ExcelRead('LOU_IMB_MagdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    LOU_SPURS = ExcelRead('LOU_SPURS',3,'row','cell',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    LOU_SPURS_foffset_spurs_Hz = LOU_SPURS{1};
    LOU_SPURS_g_spurs_dBc1Hz = LOU_SPURS{2};
    
    flagT4_LOU_PN = ExcelRead('flagT4_LOU_PN',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    flagT4_LOU_IMB = ExcelRead('flagT4_LOU_IMB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    flagT4_LOU_SPURS = ExcelRead('flagT4_LOU_SPURS',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    flagT4_LOU_QEC = ExcelRead('flagT4_LOU_QEC',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    % parameter assignment
    flagT4_LOU.AMtoPMPhsDriftDeg = LOU_AMtoPMPhsDriftDeg; %% 2020-03-14, g3, Add PhaseDriftDeg for AM/PM
    
    % ********** LO Level & Frequency: **********
    LOU.leveldB = LOU_leveldB;
    %     LOU.fLO = 400e6-mean(bwInbandB6(:));
    LOU.fLO = fix(mean(fRF)/df)*df;
    LOU.fLO = 300e6;
    LOU.fLO = LOU_fLO;
    
    % ********** LO Phase Noise input: **********
    LOU.PN = [];
    LOU.PN.ThetaDeg = LOU_PN_ThetaDeg; % %% 2020-3-14, if PN_LOU.PN_ThetaDeg==0, no-assignment
    LOU.PN.MagDriftdB1Hz = LOU_PN_MagDriftdB1Hz;
    
    LOU.PN.f_offset_Hz = LOU_PN_f_offset_Hz; % phase noise spectrum, frequencies
    LOU.PN.g_offset_dBc1Hz = LOU_PN_g_offset_dBc1Hz; % phase noise spectrum, magnitude
    
    % ********** LO IQ Imbalance input: **********
    LOU.IMB = [];
    LOU.IMB.PhsDeg = LOU_IMB_PhsDeg;
    LOU.IMB.MagdB = LOU_IMB_MagdB;
    
    % ********** LO SPURS input: **********
    LOU.SPURS = [];
    LOU.SPURS.foffset_spurs_Hz = LOU_SPURS_foffset_spurs_Hz; % discrete spurs, freq relative to fLO
    LOU.SPURS.g_spurs_dBc1Hz = LOU_SPURS_g_spurs_dBc1Hz; % discrete spurs, power relative to fLO
    
    edit SystemSim_LoadCoefficient.m
    
    %% T4a. LO Generated Perfactly!!
    %% T4b. LO Generated with Phase Noise and IQ Imbalance
    flagT4_LOU.PhsNoise = flagT4_LOU_PN; % LO wo PN
    flagT4_LOU.IMB = flagT4_LOU_IMB; % LO wo IQ imbalance
    flagT4_LOU.SPURS = flagT4_LOU_SPURS; % LO wo SPURS
    %% T4c. Check LO Phase Noise Pwr
    %% T4d. QEC for the Imbalance LO
    flagT4_LOU.QEC = flagT4_LOU_QEC;
    fnum=fnum+1;
    if flagT4_LOU_PN==0
        [loU_ideal,loU_realistic,tableT4_LOUInput,tableT4_LOU] = SYM_LOgenApp(LOU, fs, Nsamps, flagT4_LOU, fnum, 'half', 'LOU', fnum_dir);
    else
        [loU_ideal,loU_realistic,tableT4_LOUInput,tableT4_LOU] = SYM_LOgenApp(LOU, fs, Nsamps, flagT4_LOU, fnum, 'semilogx', 'LOU', fnum_dir);
    end
    
    %% 2020-07-31,
    %% T4e. LO with Unlinearity
    flagT4_LOU_Unlinearity = ExcelRead('flagT4_LOU_UNL',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagT4_LOU_Unlinearity,'off')
        % input
        LOU_GaindB = ExcelRead('LOU_GaindB',1,'colloD_idealumn','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        LOU_NFdB = ExcelRead('LOU_NFdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        LOU_Poip3dB = ExcelRead('LOU_Poip3dB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        LOU_PHD2dB = ExcelRead('LOU_PHD2dB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        LOU_PDCdB = ExcelRead('LOU_PDCdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        LOU_AMtoPMDegDrift = ExcelRead('LOU_AMtoPMDegDrift',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        LOU_OP1dB = ExcelRead('LOU_OP1dB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        LOU_IpwrdB_Target = ExcelRead('LOU_IpwrdB_Target',[1],'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        LOU_bwACLROffset = ExcelRead('LOU_bwACLROffset',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        
        fnum = fnum+1;
        [loU_realistic_UNL, ~,  tableT4_LOU_UNL] = SYM_AMPApp(loU_realistic, LOU_GaindB, LOU_NFdB, LOU_Poip3dB, LOU_PHD2dB, LOU_PDCdB, LOU_AMtoPMDegDrift, LOU_OP1dB, fs, LOU_IpwrdB_Target, [], LOU_bwACLROffset, fnum, {flagT4_LOU_Unlinearity}, fnum_dir);
        % export
        loU_output = loU_realistic_UNL;
    else
        % export
        loU_output = loU_realistic;
    end
    
    %% T5. Upconversion with LO
    flagT5_UpConv = ExcelRead('flagT5_UpConv',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagT5_UpConv,'off')
        %% T5a. QEC Estimation from TXCarrier
        flagT5_TXQECEstimation = ExcelRead('flagT5_TXQECEstimation',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        %         bwACLROffsetT5 = ExcelRead('bwACLROffsetT5',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        bwACLROffsetT5 = sum(bwChNC);
        
        if strcmp(flagT4_LOU.QEC,'on')||(flagT4_LOU.IMB==0)
            flagT5_TXQECEstimation='off';
        end
        
        flagT5_EVMcalc = bwInband+LOU.fLO;
        fnum=fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{'T5'},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        [waveformT5cell_Upconv,bwInbandT5,EVMt5_UpConv,tableT5_UpConv] = SYM_MixerApp(waveformCell, loU_output, loU_ideal, LOU.fLO, bwInband, fs, flagT5_UpConv, flagT5_EVMcalc, flagT5_TXQECEstimation, bwACLROffsetT5, fnum, NameCelltitle, fnum_dir);
        
        %% 2020-08-01, TEST N*InterModulation, 'off', default
        flagT5_TEST_nIM='off';
        if strcmp(flagT5_TEST_nIM,'on')
            fcw = 1e9;
            t = (0:Nsamps-1)/fs;
            CW_tone = sin(2*pi*fcw*t).';
            [wf_Mixer, wf_I, wf_Q] = Mixer_Up_Down_Convert_g(CW_tone, loU_output(1,:), loU_output(2,:), fs, ['TEST n*InterModulation'], 'Up',[],fnum*101);
        end
        
        % export
        waveformCell = waveformT5cell_Upconv;
        bwInband = bwInbandT5;
        ACLRdBLeft_T5UpConv_min = fix(tableT5_UpConv.tableACLRorImage.ACLRdBLeft_min);
        ACLRdBRight_T5UpConv_min = fix(tableT5_UpConv.tableACLRorImage.ACLRdBRight_min);
        flagT5_UpConv = [NameCelltitle{:}];
        
        % excel
        if strcmp(flagT5_TXQECEstimation,'TXQEC')
            ExcelWrite(round(tableT5_UpConv.tableQECest.QECest_MagdB_mean,2),'QECest_MagdB_T5',size(tableT5_UpConv.tableQECest.QECest_MagdB_mean),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
            ExcelWrite(round(tableT5_UpConv.tableQECest.QECest_PhsDeg_mean,2),'QECest_PhsDeg_T5',size(tableT5_UpConv.tableQECest.QECest_PhsDeg_mean),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        end
        if exist('bwACLROffsetT5','var')&&~isempty(bwACLROffsetT5)
            ExcelWrite(ACLRdBLeft_T5UpConv_min,'ACLRdBLeft_T5UpConv_min',size(ACLRdBLeft_T5UpConv_min),'c1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
            ExcelWrite(ACLRdBRight_T5UpConv_min,'ACLRdBRight_T5UpConv_min',size(ACLRdBRight_T5UpConv_min),'c1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        end
        ExcelWrite(round(EVMt5_UpConv,2),'EVMt5_UpConv',size(EVMt5_UpConv),'c1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    %% T5b. Apply BPF after Upconversion and before PA to filter the Image
    flagT5b_MixerBPF = ExcelRead('flagT5b_MixerBPF',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagT5b_MixerBPF,'off')
        % input
        FIRt5b_Wtype = ExcelRead('FIRt5b_Wtype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRt5b_Ftype = ExcelRead('FIRt5b_Ftype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRt5b_Order = ExcelRead('FIRt5b_Order',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRt5b_fTolerance = ExcelRead('FIRt5b_fTolerance',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRt5b_K_AttdB = ExcelRead('FIRt5b_K_AttdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRt5b_K_fdelta = ExcelRead('FIRt5b_K_fdelta',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRt5b_fcutoffL = ExcelRead('FIRt5b_fcutoffL',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRt5b_fcutoffH = ExcelRead('FIRt5b_fcutoffH',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        
        FIRt5b_Wtype = 'HAN'
        FIRt5b_Order = 101
        FIRt5b_fTolerance = 60e6;
        FIRt5b_K_fdelta = 5e6;
        FIRt5b_fcutoffL = min(bwInband)
        FIRt5b_fcutoffH = max(bwInband)
        
        FIRt5b_MixBPF = SYM_FIRApp(FIRt5b_Wtype,FIRt5b_Ftype,FIRt5b_Order,FIRt5b_K_AttdB,FIRt5b_K_fdelta,FIRt5b_fTolerance,FIRt5b_fcutoffL,FIRt5b_fcutoffH,df,bwInband,NCarriers,fs,fnum_Debug);
        
        fnum = fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagT5b_MixerBPF},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        [waveformT5Acell_MixerBPF, b_t5bMixerBPF, EVMt5b_MixerBPF, tableT5b_MixerBPF] = SYM_FilterApp(waveformCell, FIRt5b_MixBPF, fs, mehtod_filter, bwInband, fnum, NameCelltitle, fnum_dir);
        [ACLRdB_T5bMixerBPF, IpwrdB_T5bMixerBPF] = ACLR_calc_g(waveformT5Acell_MixerBPF{1,:}, fs, bwInband, bwACLROffsetT5, fnum*200, NameCelltitle, fnum_dir); hold on
        
        % export
        waveformCell= waveformT5Acell_MixerBPF;
        ACLRdBLeft_T5aMixerBPF_min = fix(min(ACLRdB_T5bMixerBPF.ACLRdBLeft));
        ACLRdBRight_T5aMixerBPF_min = fix(min(ACLRdB_T5bMixerBPF.ACLRdBRight));
        
        % excel
        ExcelWrite(tableT5b_MixerBPF.tableInput.FIR_Order,'FIRt5b_Order_out',size(tableT5b_MixerBPF.tableInput.FIR_Order),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(EVMt5b_MixerBPF,2),'EVMt5b_MixerBPF',size(EVMt5b_MixerBPF),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(ACLRdBLeft_T5aMixerBPF_min,'ACLRdBLeft_T5bMixerBPF_min',1,'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(ACLRdBRight_T5aMixerBPF_min,'ACLRdBRight_T5bMixerBPF_min',1,'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    %% T6. Add Gain/IM3/IM2orDCoffset/NF
    flagT6_addGainIM3NF = ExcelRead('flagT6_addGainIM3NF',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagT6_addGainIM3NF,'off')
        % input
        GaindB_T6 = ExcelRead('GaindB_T6',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        Poip3dB_T6 = ExcelRead('Poip3dB_T6',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        Pim2dB_T6 = ExcelRead('Pim2dB_T6',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        PDCdB_T6 = ExcelRead('PDCdB_T6',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        AMtoPMDegDrift_T6 = ExcelRead('AMtoPMDegDrift_T6',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        NFdB_T6 = ExcelRead('NFdB_T6',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        %% 2020-05-19, Compare results of different Target power
        IpwrdB_T6Target = ExcelRead('IpwrdB_T6Target',[1],'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        %         bwACLROffsetT6 = ExcelRead('bwACLROffsetT6',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        bwACLROffsetT6 = sum(bwChNC);
        
        fnum = fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagT6_addGainIM3NF},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        %         [waveformT6cell_wiGainOIP3NF, EVMt6_GainOIP3NF,  tableT6_GainOIP3NF] = SYM_AMPApp(waveformCell, GaindB_T6, NFdB_T6, Poip3dB_T6, Pim2dB_T6, PDCdB_T6, AMtoPMDegDrift_T6, [], fs, IpwrdB_T6Target, bwInband, bwACLROffsetT6, fnum, NameCelltitle, fnum_dir);
        [waveformT6cell_wiGainOIP3NF, EVMt6_GainOIP3NF,  tableT6_GainOIP3NF] = SYM_AMPApp(waveformCell, GaindB_T6, NFdB_T6, Poip3dB_T6, Pim2dB_T6, PDCdB_T6, AMtoPMDegDrift_T6, [], fs, IpwrdB_T6Target, bwInband, bwACLROffsetT6, fnum, NameCelltitle, fnum_dir);
        
        for idC=1:NCarriers
            [IpwrdB_T5, ~, ~] = Pwr_Inband_g(fft(waveformCell{idC,:}, length(waveformCell{idC,:}), DIMFFT), fs, bwInband, [], 'full', 0);
            [IpwrdB_T6, ~, ~] = Pwr_Inband_g(fft(waveformT6cell_wiGainOIP3NF{idC,:}, length(waveformT6cell_wiGainOIP3NF{idC,:}), DIMFFT), fs, bwInband, [], 'full', 0);
            IpwrdB1Hz_Noise_T5 = Pwr_Inband_g(fft(waveformCell{:}, length(waveformCell{:}), DIMFFT), fs, fs/2-10e6+[0 2*df], [], 'full', 0)-10*log10(2*df);
            IpwrdB1Hz_Noise_T6 = Pwr_Inband_g(fft(waveformT6cell_wiGainOIP3NF{:}, length(waveformT6cell_wiGainOIP3NF{:}), DIMFFT), fs, fs/2-10e6+[0 2*df], [], 'full', 0)-10*log10(2*df);
            IpwrdB_T5_mean(idC,:)=round(mean(IpwrdB_T5),2);
            IpwrdB_T6_mean(idC,:)=round(mean(IpwrdB_T6),2);
            IpwrdB1Hz_Noise_T5_mean(idC,:)=round(mean(IpwrdB1Hz_Noise_T5),2);
            IpwrdB1Hz_Noise_T6_mean(idC,:)=round(mean(IpwrdB1Hz_Noise_T6),2);
            
        end
        
        % export
        waveformCell = waveformT6cell_wiGainOIP3NF;
        ACLRdBLeft_T6_min = fix(tableT6_GainOIP3NF.tableACLR.ACLRdBLeft_min);
        ACLRdBRight_T6_min = fix(tableT6_GainOIP3NF.tableACLR.ACLRdBRight_min);
        
        % excel
        ExcelWrite(round(EVMt6_GainOIP3NF,2),'EVMt6_GainOIP3NF',size(EVMt6_GainOIP3NF),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(ACLRdBLeft_T6_min,'ACLRdBLeft_T6GainOIP3NF_min',[],'c1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(ACLRdBRight_T6_min,'ACLRdBRight_T6GainOIP3NF_min',1,'r1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(IpwrdB_T5_mean,'IpwrdB_T5_mean',1,'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(IpwrdB_T6_mean,'IpwrdB_T6_mean',1,'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    
    %% 2021-03-16, PA block
    flagT6b_PAFinal = ExcelRead('flagT6b_PAFinal',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagT6b_PAFinal,'off')
        % input
        GaindB_T6b = ExcelRead('GaindB_T6b',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        Poip3dB_T6b = ExcelRead('Poip3dB_T6b',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        Pim2dB_T6b = ExcelRead('Pim2dB_T6b',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        PDCdB_T6b = ExcelRead('PDCdB_T6b',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        AMtoPMDegDrift_T6b = ExcelRead('AMtoPMDegDrift_T6b',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        NFdB_T6b = ExcelRead('NFdB_T6b',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        %% 2020-05-19, Compare results of different Target power
        IpwrdB_T6bTarget = ExcelRead('IpwrdB_T6bTarget',[1],'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        %         bwACLROffsetT6 = ExcelRead('bwACLROffsetT6',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        bwACLROffsetT6b = bwACLROffsetT6;
        
        IpwrdB_T6bTarget
        if IpwrdB_T6_mean>IpwrdB_T6bTarget
            error('Pwr of stage Driver out of Target!')
        end
        
        fnum = fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagT6b_PAFinal},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        %         [waveformT6cell_wiGainOIP3NF, EVMt6_GainOIP3NF,  tableT6_GainOIP3NF] = SYM_AMPApp(waveformCell, GaindB_T6, NFdB_T6, Poip3dB_T6, Pim2dB_T6, PDCdB_T6, AMtoPMDegDrift_T6, [], fs, IpwrdB_T6Target, bwInband, bwACLROffsetT6, fnum, NameCelltitle, fnum_dir);
        [waveformT6bcell_PA, EVMt6b_PA,  tableT6b_PA] = SYM_AMPApp(waveformCell, GaindB_T6b, NFdB_T6b, Poip3dB_T6b, Pim2dB_T6b, PDCdB_T6b, AMtoPMDegDrift_T6b, [], fs, IpwrdB_T6bTarget, bwInband, bwACLROffsetT6b, fnum, NameCelltitle, fnum_dir);
        
        for idC=1:NCarriers
            [IpwrdB_T6, ~, ~] = Pwr_Inband_g(fft(waveformCell{idC,:}, length(waveformCell{idC,:}), DIMFFT), fs, bwInband, [], 'half', 0);
            [IpwrdB_T6b, ~, ~] = Pwr_Inband_g(fft(waveformT6bcell_PA{idC,:}, length(waveformT6bcell_PA{idC,:}), DIMFFT), fs, bwInband, [], 'half', 0);
            %             IpwrdB1Hz_Noise_T5 = Pwr_Inband_g(fft(waveformCell{:}, length(waveformCell{:}), DIMFFT), fs, fs/2-10e6+[0 2*df], [], 'full', 0)-10*log10(2*df);
            %             IpwrdB1Hz_Noise_T6 = Pwr_Inband_g(fft(waveformT6cell_wiGainOIP3NF{:}, length(waveformT6cell_wiGainOIP3NF{:}), DIMFFT), fs, fs/2-10e6+[0 2*df], [], 'full', 0)-10*log10(2*df);
            %             IpwrdB_T5_mean(idC,:)=round(mean(IpwrdB_T5),2);
            IpwrdB_T6b_mean(idC,:)=round(mean(IpwrdB_T6b),2);
            %             IpwrdB1Hz_Noise_T5_mean(idC,:)=round(mean(IpwrdB1Hz_Noise_T5),2);
            %             IpwrdB1Hz_Noise_T6_mean(idC,:)=round(mean(IpwrdB1Hz_Noise_T6),2);
        end
        
        % export
        waveformCell = waveformT6bcell_PA;
        ACLRdBLeft_T6b_min = fix(tableT6b_PA.tableACLR.ACLRdBLeft_min);
        ACLRdBRight_T6b_min = fix(tableT6b_PA.tableACLR.ACLRdBRight_min);
        
        % excel
        ExcelWrite(round(EVMt6b_PA,2),'EVMt6b_PA',size(EVMt6b_PA),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(ACLRdBLeft_T6b_min,'ACLRdBLeft_T6bPA_min',[],'c1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(ACLRdBRight_T6b_min,'ACLRdBRight_T6bPA_min',1,'r1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(IpwrdB_T6b_mean,'IpwrdB_T6b_mean',1,'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    %% T7. Add PhaseShift for each branch
    flagT7_addPhaseShift = ExcelRead('flagT7_addPhaseShift',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if (strcmp(flagC0_AnalysisWF,'AC')||strcmp(flagC0_AnalysisWF,'DLAC'))&&~strcmp(flagT7_addPhaseShift,'off') % all of PhsShiftNbrDeg are ZERO
        % input
        ACPhsShiftNbrDeg_T7 = ExcelRead('ACPhsShiftNbrDeg_T7',4,'row','vector',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        
        waveformT7_addPhsShift = exp(1i*ACPhsShiftNbrDeg_T7./180*pi).*waveformCell{:};
        % export
        waveformCell = {waveformT7_addPhsShift};
    end
    
    %% T8. Combine all branches for AC
    flagT8_CombACNbr = ExcelRead('flagT8_CombACNbr',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if (strcmp(flagC0_AnalysisWF,'AC')||strcmp(flagC0_AnalysisWF,'DLAC'))&&~strcmp(flagT8_CombACNbr,'off')
        waveformT8_CombACNbr = sum(waveformCell{1,:},mod(DIMFFT,2)+1); % all branch combination
        Nbr = size(waveformT8_CombACNbr,mod(DIMFFT,2)+1); % combine to One branch
        if strcmp(flagC0_AnalysisWF,'DLAC')
            % export to BRComb and BRwoComb, BRwoComb for DL demodulation
            waveformT8_DL = cell2mat(waveformCell);
            waveformT8_DL_NbrDL = waveformT8_DL(:,1:Nbr_DL);
            waveformT8Cell = [{waveformT8_CombACNbr};{waveformT8_DL_NbrDL}];
            NameCellT8={'DLAC' newline 'CombACNbr';'DLAC' newline 'NoCombNbr'};
            
        else
            % export
            waveformT8Cell = {waveformT8_CombACNbr};
            NameCellT8={'ACBRComb'};
        end
        
        %% 2021-03-15, Check AC Combination: Ipwr of carrier and noise floor
        idBR = 1;
        fnum(1) = fnum(1)+1;
        PLOT_FFT_dB_g(waveformT7_addPhsShift(:,idBR), fs, Nsamps, ['DLAC no.',num2str(idBR),' br'], 'df', 'full', 'pwr', fnum);
        PLOT_FFT_dB_g(waveformT8_CombACNbr, fs, Nsamps, ['DLAC all ',num2str(size(waveformT7_addPhsShift,2)),' brs combination'], 'df', 'full', 'pwr', fnum);
        [IpwrdB_T8_1BR, ~, ~] = Pwr_Inband_g(fft(waveformT7_addPhsShift(:,idBR), Nsamps, DIMFFT), fs, bwInband, [], 'full', []);
        [IpwrdB_T8_8BR, ~, ~] = Pwr_Inband_g(fft(waveformT8_CombACNbr, Nsamps, DIMFFT), fs, bwInband, [], 'full', []);
        [IpwrdB_T8_1BR_NoiseFloor, ~, ~] = Pwr_Inband_g(fft(waveformT7_addPhsShift(:,idBR), Nsamps, DIMFFT), fs, bwInband+50e6, [], 'full', []);
        [IpwrdB_T8_8BR_NoiseFloor, ~, ~] = Pwr_Inband_g(fft(waveformT8_CombACNbr, Nsamps, DIMFFT), fs, bwInband+50e6, [], 'full', []);
        
        waveformCell = waveformT8Cell;
        NameCell = NameCellT8;
        NCarriers = size(waveformCell,1);
        % excel
        NameCell = ExcelWrite(NameCell,'CarriersT8_CombACNbr',size(NameCell),'r1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(NCarriers,'NCarriersT8_CombACNbr',size(NCarriers),'r1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    else
        NameCellT8=NameCell;
    end
    
    %% T9. Export to waveformTX
    waveformCell_TX = waveformCell;
    fsTX = fs;
    NsampsTX = Nsamps;
    bwInbandTX = bwInband;
    NameCellTX = NameCellT8;
    clear NsampsRX NsampsDP waveformD7_ACPhsAlign
    
end

%% Receiver Block =============================================================================================================================
if ~exist('NsampsRX','var')
    disp('Receiver Block ====================================================')
    %% R0. Transmitter Signal import to Receiver Block
    waveformCell = waveformCell_TX;
    fs = fsTX;
    Nsamps = NsampsTX;
    bwInband = bwInbandTX;
    NameCell = NameCellTX;
    NCarriers = size(waveformCell_TX,1);
    
    %% R0a. Add Loss for Receiver
    flagR0a_addLoss = ExcelRead('flagR01_addLoss',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagR0a_addLoss,'off')
        % input
        GaindB_R0a = ExcelRead('GaindB_R01',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        Poip3dB_R0a = ExcelRead('Poip3dB_R01',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        Pim2dB_R0a = ExcelRead('Pim2dB_R01',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        PDCdB_R0a = ExcelRead('PDCdB_R01',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        AMtoPMDegDrift_R0a = ExcelRead('AMtoPMDegDrift_R01',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        NFdB_R0a = ExcelRead('NFdB_R01',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        %% 2020-05-19, Compare results of different Target power
        IpwrdB_R0aTarget = ExcelRead('IpwrdB_R01Target',[1],'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        %         bwACLROffsetR0a = ExcelRead('bwACLROffsetR01',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        bwACLROffsetR0a = sum(bwChNC);
        
        % generate NameCelltitle
        NameCellflag=repmat([{flagR0a_addLoss},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        fnum=fnum+1;
        [waveformR0Acell_addLoss, EVMr0a_Loss,  tableR0a_Loss] = SYM_AMPApp(waveformCell, GaindB_R0a, NFdB_R0a, Poip3dB_R0a, Pim2dB_R0a, PDCdB_R0a, AMtoPMDegDrift_R0a, [], fs, IpwrdB_R0aTarget, bwInband, bwACLROffsetR0a, fnum, NameCelltitle, fnum_dir);
        for idC=1:NCarriers
            [IpwrdB_R0, ~, ~] = Pwr_Inband_g(fft(waveformCell{idC,:}, length(waveformCell{idC,:}), DIMFFT), fs, bwInband, [], 'half', 0);
            [IpwrdB_R0a, ~, ~] = Pwr_Inband_g(fft(waveformR0Acell_addLoss{idC,:}, length(waveformR0Acell_addLoss{idC,:}), DIMFFT), fs, bwInband, [], 'half', 0);
            IpwrdB1Hz_Noise_R0 = Pwr_Inband_g(fft(waveformCell{idC,:}, length(waveformCell{idC,:}), DIMFFT), fs, fs/2-10e6+[0 2*df], [], 'full', 0)-10*log10(2*df);
            IpwrdB1Hz_Noise_R0a = Pwr_Inband_g(fft(waveformR0Acell_addLoss{idC,:}, length(waveformR0Acell_addLoss{idC,:}), DIMFFT), fs, fs/2-10e6+[0 2*df], [], 'full', 0)-10*log10(2*df);
            IpwrdB_R0_mean(idC,:)=round(mean(IpwrdB_R0),2);
            IpwrdB_R0a_mean(idC,:)=round(mean(IpwrdB_R0a),2);
            IpwrdB1Hz_Noise_R0_mean(idC,:)=round(mean(IpwrdB1Hz_Noise_R0),2);
            IpwrdB1Hz_Noise_R0a_mean(idC,:)=round(mean(IpwrdB1Hz_Noise_R0a),2);
        end
        
        % export
        waveformCell = waveformR0Acell_addLoss;
        if ~isempty(bwACLROffsetR0a)
            ACLRdBLeft_R0a_min = fix(tableR0a_Loss.tableACLR.ACLRdBLeft_min);
            ACLRdBRight_R0a_min = fix(tableR0a_Loss.tableACLR.ACLRdBRight_min);
            ExcelWrite(ACLRdBLeft_R0a_min,'ACLRdBLeft_R01Loss_min',size(ACLRdBLeft_R0a_min),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
            ExcelWrite(ACLRdBRight_R0a_min,'ACLRdBRight_R01Loss_min',size(ACLRdBRight_R0a_min),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        end
        
        % excel
        %         ExcelWrite(round(EVMr01_Loss,2),{'EVMr01_Loss','EVMr01_Loss_BRAll'},size(EVMr01_Loss),'c1row',[],xlsFile,SheetName,RangeArrayInput);
        ExcelWrite(round(EVMr0a_Loss,2),{'EVMr01_Loss'},size(EVMr0a_Loss),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(IpwrdB_R0_mean,'IpwrdB_R0_mean',size(IpwrdB_R0a_mean),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(IpwrdB_R0a_mean,'IpwrdB_R01_mean',size(IpwrdB_R0a_mean),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    %% 2020-07-24, Add Blocking for RX
    %% R0b. Transmitter Signal import to Receiver Block
    flagR0b_addBlocking = ExcelRead('flagR0b_addBlocking',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
%     flagR0b_addBlocking = 'R0bBLKing'
    method_flagR0b_addBlocking = 2;
    
    if ~strcmp(flagR0b_addBlocking,'off')
        
        switch method_flagR0b_addBlocking
            case {1}
                % generate NameCelltitle
                NameCellflag=repmat([{flagR0b_addBlocking},{', '}],size(waveformCell,1),1);
                NameCelltitle=[NameCellflag,NameCell];
                fnum=fnum+1;
                
                flagR0b_BlockingInput = 'excel';
                flagR0b_BlockingInput = 'IM_CW&WB5MHz_Offset7p125&17p5MHz';
                switch flagR0b_BlockingInput
                    case {'excel'}
                        %                 Blocking_Input1='Excel_Blocking_co-lcation.xlsx'; % Excel File
                        Blocking_Input1 = ExcelRead('Blocking_Input1',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
                        %                 Blocking_Input2 = 'B42'; % Excel Sheet Name
                        Blocking_Input2 = ExcelRead('Blocking_Input2',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
                        %                 Blocking_Input3 = 6; % Excel Column
                        Blocking_Input3 = ExcelRead('Blocking_Input3',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
                    case {'~excel'}
                        Blocking_Input1 = [100 200]*1e6; % Blocking Freq. Inband Hz
                        Blocking_Input2 = [30]; % Blocking Level dB
                        Blocking_Input3 = {'Blocking'}; % Blocking Name
                    case {'Inblocking_BW5MHz_Offset7p5MHz'}
                        Blocking_Input1 = bwInband(end)+7.5e6+5e6*[-1 1]/2
                        Blocking_Input2 = -35; % -35/-40/-33(NBIoT)
                        Blocking_Input3 = {'Inblocking_BW5MHz_Offset7p5MHz'}; % Blocking Name
                    case {'Inblocking_BW5MHz_Offset2p5MHz'}
                        Blocking_Input1 = bwInband(end)+2.5e6+5e6*[-1 1]/2
                        Blocking_Input2 = -42; % -42/-47
                        Blocking_Input3 = {'Inblocking_BW5MHz_Offset2p5MHz'}; % Blocking Name
                    case {'Inblocking_BW1p4MHz_Offset2p1MHz'}
                        Blocking_Input1 = bwInband(end)+2.1e6+1.4e6*[-1 1]/2
                        Blocking_Input2 = -35; % -42/-47
                        Blocking_Input3 = {'Inblocking_BW1p4MHz_Offset2p1MHz'}; % Blocking Name
                    case {'NarInblocking_BW1RB180KHz_Offset240KHz_m0'}
                        m=0;
                        Blocking_Input1 = bwInband(end)+(240e3+m*180e3)+180e3*[-1 1]/2
                        Blocking_Input2 = -40; % -40/-49/-39(NBIoT)
                        Blocking_Input3 = {'NarInblocking_BW1RB180KHz_Offset240KHz_m0'}; % Blocking Name
                    case {'IM_CW&WB5MHz_Offset7p125&17p5MHz'}
                        Blocking_Input1(1,:) = bwInband(end)+7.125e6+df*[-1 1]/2 % CW
                        Blocking_Input1(2,:) = bwInband(end)+17.5e6+5e6*[-1 1]/2 % WB5MHz
                        
                        Blocking_Input2(1,:) = -30; % -40/-49/-39(NBIoT)
                        Blocking_Input2(2,:) = -40; % -40/-49/-39(NBIoT)
                        
                        Blocking_Input3 = {'IMblocking'}; % Blocking Name
                end
                
                [waveformR0Bcell_addBlocking, EVMr0b_Blocking, tableR0bInput_Blocking] = SYM_BlockingApp(waveformCell, Blocking_Input1, Blocking_Input2, Blocking_Input3, fs, bwInband, fnum, NameCelltitle, fnum_dir);
                
            case {2}
                InterferType_R0b = ExcelRead('InterferType_R0b',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
                switch InterferType_R0b
                    case {'NarrowBandBlocking'}
                        InterferConfig_R0b.Type = InterferType_R0b;
                        InterferConfig_R0b.bw_Channel = '180kHz';
                        InterferConfig_R0b.NRB = 1;
                        InterferConfig_R0b.Scs_kHz = 15;
                        InterferConfig_R0b.MOD = 'QPSK';
                        InterferConfig_R0b.Carrier_Type = 'LTE';
                        InterferConfig_R0b.df = df
                        InterferConfig_R0b.fs = fs;
                        InterferConfig_R0b.SamplesDecimation = DLconfig(1).SamplesDecimation;
                        InterferOffsetHz_R0b = 240e3;
                        InterferLeveldBm_R0b = -40;
                end
                fnum = fnum+1;
                [waveformR0Bcell_addBlocking,~,BLKConfig_R0b,EVMr0b_Blocking] = SYM_BlockingApp2(waveformCell, InterferType_R0b, InterferConfig_R0b, InterferOffsetHz_R0b, InterferLeveldBm_R0b, bwInband, InterferConfig_R0b.fs, bwInband, fnum, NameCelltitle, fnum_dir);
                
        end
        
        % export
        waveformCell = waveformR0Bcell_addBlocking;
        % excel
        ExcelWrite(round(EVMr0b_Blocking,2),{'EVMr0b_Blocking'},size(EVMr0b_Blocking),'c1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    
    %% R1. Apply RXBPF Filter
    flagR1_RXBPF = ExcelRead('flagR1_RXBPF',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagR1_RXBPF,'off')
        % input
        FIRr1_Wtype = ExcelRead('FIRr1_Wtype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRr1_Ftype = ExcelRead('FIRr1_Ftype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRr1_Order = ExcelRead('FIRr1_Order',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRr1_fTolerance = ExcelRead('FIRr1_fTolerance',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRr1_K_AttdB = ExcelRead('FIRr1_K_AttdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRr1_K_fdelta = ExcelRead('FIRr1_K_fdelta',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRr1_fcutoffL = ExcelRead('FIRr1_fcutoffL',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRr1_fcutoffH = ExcelRead('FIRr1_fcutoffH',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        
        FIRr1_RXBPF = SYM_FIRApp(FIRr1_Wtype,FIRr1_Ftype,FIRr1_Order,FIRr1_K_AttdB,FIRr1_K_fdelta,FIRr1_fTolerance,FIRr1_fcutoffL,FIRr1_fcutoffH,df,bwInband,NCarriers);
        
        fnum = fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagR1_RXBPF},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        [waveformR1cell_RXBPF, b_R1RXBPF, EVMr1_RXBPF, tableR1_RXBPF] = SYM_FilterApp(waveformCell, FIRr1_RXBPF, fs, mehtod_filter, bwInband, fnum, NameCelltitle, fnum_dir);
        
        % export
        waveformCell= waveformR1cell_RXBPF;
        
        % excel
        ExcelWrite(tableR1_RXBPF.tableInput.FIR_Order,'FIRr1_Order_out',size(tableR1_RXBPF.tableInput.FIR_Order),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(EVMr1_RXBPF,2),'EVMr1_RXBPF',size(EVMr1_RXBPF),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    %% R1b Implement LNA, 2020-05-23
    flagR1a_LNA = ExcelRead('flagR1a_LNA',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagR1a_LNA,'off')
        % input
        GaindB_R1a = ExcelRead('GaindB_R1a',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        Poip3dB_R1a = ExcelRead('Poip3dB_R1a',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        Pim2dB_R1a = ExcelRead('Pim2dB_R1a',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        PDCdB_R1a = ExcelRead('PDCdB_R1a',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        AMtoPMDegDrift_R1a = ExcelRead('AMtoPMDegDrift_R1a',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        NFdB_R1a = ExcelRead('NFdB_R1a',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        %% 2020-05-19, Compare results of different Target power
        IpwrdB_R1aTarget = ExcelRead('IpwrdB_R1aTarget',[1],'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        OP1dB_R1a = ExcelRead('OP1dB_R1a',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        %         bwACLROffsetR1a = ExcelRead('bwACLROffsetR01',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        bwACLROffsetR1a = sum(bwChNC);
        
        
        % generate NameCelltitle
        NameCellflag=repmat([{flagR1a_LNA},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        fnum=fnum+1;
        [waveformR1Acell_LNA, EVMr1a_LNA,  tableR1a_LNA] = SYM_AMPApp(waveformCell, GaindB_R1a, NFdB_R1a, Poip3dB_R1a, Pim2dB_R1a, PDCdB_R1a, AMtoPMDegDrift_R1a, OP1dB_R1a, fs, [], bwInband, bwACLROffsetR1a, fnum, NameCelltitle, fnum_dir);
        for idC=1:NCarriers
            [IpwrdB_R1, ~, ~] = Pwr_Inband_g(fft(waveformCell{idC,:}, length(waveformCell{idC,:}), DIMFFT), fs, bwInband, [], 'half', 0);
            [IpwrdB_R1a, ~, ~] = Pwr_Inband_g(fft(waveformR1Acell_LNA{idC,:}, length(waveformR1Acell_LNA{idC,:}), DIMFFT), fs, bwInband, [], 'half', 0);
            IpwrdB1Hz_Noise_R1 = Pwr_Inband_g(fft(waveformCell{idC,:}, length(waveformCell{idC,:}), DIMFFT), fs, fs/2-10e6+[0 2*df], [], 'full', 0)-10*log10(2*df);
            IpwrdB1Hz_Noise_R1a = Pwr_Inband_g(fft(waveformR1Acell_LNA{idC,:}, length(waveformR1Acell_LNA{idC,:}), DIMFFT), fs, fs/2-10e6+[0 2*df], [], 'full', 0)-10*log10(2*df);
            IpwrdB_R1_mean(idC,:)=round(mean(IpwrdB_R1),2);
            IpwrdB_R1a_mean(idC,:)=round(mean(IpwrdB_R1a),2);
            IpwrdB1Hz_Noise_R1_mean(idC,:)=round(mean(IpwrdB1Hz_Noise_R1),2);
            IpwrdB1Hz_Noise_R1a_mean(idC,:)=round(mean(IpwrdB1Hz_Noise_R1a),2);
        end
        % export
        waveformCell = waveformR1Acell_LNA;
        if ~isempty(bwACLROffsetR1a)
            ACLRdBLeft_R1a_min = fix(tableR1a_LNA.tableACLR.ACLRdBLeft_min);
            ACLRdBRight_R1a_min = fix(tableR1a_LNA.tableACLR.ACLRdBRight_min);
            % excel
            ExcelWrite(ACLRdBLeft_R1a_min,'ACLRdBLeft_R1a_min',size(ACLRdBLeft_R1a_min),'c1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
            ExcelWrite(ACLRdBRight_R1a_min,'ACLRdBRight_R1a_min',size(ACLRdBRight_R1a_min),'c1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        end
        
        % excel
        ExcelWrite(round(EVMr1a_LNA,2),{'EVMr1a_LNA'},size(EVMr1a_LNA),'r1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(IpwrdB_R1_mean,'IpwrdB_R1_mean',size(IpwrdB_R1_mean),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(IpwrdB_R1a_mean,'IpwrdB_R1a_mean',size(IpwrdB_R1a_mean),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    
    %% R2. LOD Generator
    edit SystemSim_LoadCoefficient.m
    
    % excel import
    LOD_AMtoPMPhsDriftDeg = ExcelRead('LOD_AMtoPMPhsDriftDeg',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    LOD_leveldB = ExcelRead('LOD_leveldB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    LOD_fLO = ExcelRead('LOD_fLO',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    LOD_PN_ThetaDeg = ExcelRead('LOD_PN_ThetaDeg',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    LOD_PN_MagDriftdB1Hz = ExcelRead('LOD_PN_MagDriftdB1Hz',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    LOD_PN_offset = ExcelRead('LOD_PN_offset',5,'row','cell',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    LOD_PN_f_offset_Hz = LOD_PN_offset{1};
    LOD_PN_g_offset_dBc1Hz = LOD_PN_offset{2};
    
    LOD_IMB_PhsDeg = ExcelRead('LOD_IMB_PhsDeg',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    LOD_IMB_MagdB = ExcelRead('LOD_IMB_MagdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    LOD_SPURS = ExcelRead('LOD_SPURS',3,'row','cell',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    LOD_SPURS_foffset_spurs_Hz = LOD_SPURS{1};
    LOD_SPURS_g_spurs_dBc1Hz = LOD_SPURS{2};
    
    flagR2_LOD_PN = ExcelRead('flagR2_LOD_PN',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    flagR2_LOD_IMB = ExcelRead('flagR2_LOD_IMB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    flagR2_LOD_SPURS = ExcelRead('flagR2_LOD_SPURS',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    flagR2_LOD_QEC = ExcelRead('flagR2_LOD_QEC',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    % parameter assignment
    flagR2_LOD.AMtoPMPhsDriftDeg = LOD_AMtoPMPhsDriftDeg; %% 2020-03-14, g3, Add PhaseDriftDeg for AM/PM
    
    % ********** LO Level & Frequency: **********
    LOD.leveldB = LOD_leveldB;
    %     LOU.fLO = 400e6-mean(bwInbandB6(:));
    LOD.fLO = fix(mean(fRF)/df)*df;
    LOD.fLO = 300e6;
    LOD.fLO = LOD_fLO;
    
    % ********** LO Phase Noise input: **********
    LOD.PN = [];
    LOD.PN.ThetaDeg = LOD_PN_ThetaDeg; % %% 2020-3-14, if PN_LOU.PN_ThetaDeg==0, no-assignment
    LOD.PN.MagDriftdB1Hz = LOD_PN_MagDriftdB1Hz;
    
    LOD.PN.f_offset_Hz = LOD_PN_f_offset_Hz; % phase noise spectrum, frequencies
    LOD.PN.g_offset_dBc1Hz = LOD_PN_g_offset_dBc1Hz; % phase noise spectrum, magnitude
    
    % ********** LO IQ Imbalance input: **********
    LOD.IMB = [];
    LOD.IMB.PhsDeg = LOD_IMB_PhsDeg;
    LOD.IMB.MagdB = LOD_IMB_MagdB;
    
    % ********** LO SPURS input: **********
    LOD.SPURS = [];
    LOD.SPURS.foffset_spurs_Hz = LOD_SPURS_foffset_spurs_Hz; % discrete spurs, freq relative to fLO
    LOD.SPURS.g_spurs_dBc1Hz = LOD_SPURS_g_spurs_dBc1Hz; % discrete spurs, power relative to fLO
    
    %% R2a. LO Generated Perfactly!!
    %% R2b. LO Generated with Phase Noise and IQ Imbalance
    flagR2_LOD.PhsNoise = flagR2_LOD_PN; % LO wo PN
    flagR2_LOD.IMB = flagR2_LOD_IMB; % LO wo IQ imbalance
    flagR2_LOD.SPURS = flagR2_LOD_SPURS; % LO wo SPURS
    %% R2c. Check LO Phase Noise Pwr
    %% R2d. QEC for the Imbalance LO
    flagR2_LOD.QEC = flagR2_LOD_QEC;
    fnum=fnum+1;
    [loD_ideal,loD_realistic,tableR2_LODInput,tableR2_LOD] = SYM_LOgenApp(LOD, fs, Nsamps, flagR2_LOD, fnum, 'semilogx', 'LOD', fnum_dir);
    
    %% R2e. LO with Unlinearity
    flagR2_LOD_Unlinearity = ExcelRead('flagR2_LOD_UNL',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagR2_LOD_Unlinearity,'off')
        % input
        LOD_GaindB = ExcelRead('LOD_GaindB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        LOD_NFdB = ExcelRead('LOD_NFdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        LOD_Poip3dB = ExcelRead('LOD_Poip3dB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        LOD_PHD2dB = ExcelRead('LOD_PHD2dB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        LOD_PDCdB = ExcelRead('LOD_PDCdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        LOD_AMtoPMDegDrift = ExcelRead('LOD_AMtoPMDegDrift',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        LOD_OP1dB = ExcelRead('LOD_OP1dB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        LOD_IpwrdB_Target = ExcelRead('LOD_IpwrdB_Target',[1],'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        LOD_bwACLROffset = ExcelRead('LOD_bwACLROffset',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        
        [loD_realistic_UNL, ~,  tableR4_LOD_UNL] = SYM_AMPApp(loD_realistic, LOD_GaindB, LOD_NFdB, LOD_Poip3dB, LOD_PHD2dB, LOD_PDCdB, LOD_AMtoPMDegDrift, LOD_OP1dB, fs, LOD_IpwrdB_Target, [], LOD_bwACLROffset, fnum, {flagR2_LOD_Unlinearity});
        % export
        loD_output = loD_realistic_UNL;
    else
        % export
        loD_output = loD_realistic;
    end
    
    %% R3. Downconversion with LO
    flagR3_LODConv = ExcelRead('flagR3_LODConv',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    %% 2021-04-01, RXQEC Method Comparsion: trx_QEC_g and dsp_iqc_widely_g
    flagR3c_TEST_RXQEC = 'off';
    if ~strcmpi(flagR3c_TEST_RXQEC, 'off')
        % input
        LOD_IMB_PhsDeg = 3;
        LOD_IMB_MagdB = -7;
        LOD_fLO = LOU_fLO+fNCO(1)
        
        % generate LO
        
        % Downconversion
        bwACLROffsetR3 = sum(bwChNC);
        flagR3_EVMcalc = bwInband-LOD.fLO;
        [waveformR3cell_DnConv_R3c,bwInbandR3c,EVMr3c_DnConv,tableR3c_DnConv] = SYM_MixerApp(waveformCell, loD_output, loD_ideal, LOD.fLO, bwInband, fs, flagR3_LODConv, flagR3_EVMcalc, [], bwACLROffsetR3, fnum, {'R3c DNConv'}, fnum_dir);
        
        % QEC estimation
        bwInband_BBQEC = flagR3_EVMcalc;
        [QECest_MagdB_R3c, QECest_PhsDeg_R3c, waveformR3C_QECCorr] = trx_QEC_g(waveformR3cell_DnConv_R3c{:}, [], [], fs, [], 'RXQEC',[]);
        
        % QEC estimation2
        [waveformR3C_QECCorr2,QECest_MagdB_R3c2,QECest_PhsDeg_R3c2] = dsp_iqc_widely_g(waveformR3cell_DnConv_R3c{:});
        
        PLOT_FFT_dB_g(waveformR3cell_DnConv_R3c{:}, fs, Nsamps,  ['wf DNConv., GaindBError:',num2str(round(LOD_IMB_MagdB,2)), ', PhaseDegError:',num2str(round(LOD_IMB_PhsDeg,2))], 1, 'full', 'pwr', [20210401]);
        PLOT_FFT_dB_g(waveformR3C_QECCorr, fs, Nsamps, ['wf QEC Corr, GaindBError:',num2str(round(QECest_MagdB_R3c,2)), ', PhaseDegError:',num2str(round(QECest_PhsDeg_R3c,2))], 1, 'full', 'pwr', [20210401]);
        PLOT_FFT_dB_g(waveformR3C_QECCorr2, fs, Nsamps, ['wf QEC Corr2, GaindBError:',num2str(round(QECest_MagdB_R3c2,2)), ', PhaseDegError:',num2str(round(QECest_PhsDeg_R3c2,2))], 1, 'full', 'pwr', [20210401]);
        title('QEC Estimation and Correction')
    end
    
    if strcmp(flagR3_LODConv,'Down')
        % input
        %% R3a. QEC Estimation from RXCarrier
        flagR3_RXQECEstimation = ExcelRead('flagR3_RXQECEstimation',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        %         bwACLROffsetR3 = ExcelRead('bwACLROffsetR3',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        bwACLROffsetR3 = sum(bwChNC);
        
        if strcmp(flagR2_LOD.QEC,'on')||(flagR2_LOD.QEC~=0)||strcmp(flagC0_AnalysisWF,'AC')||strcmp(flagC0_AnalysisWF,'DLAC')
            flagR3_RXQECEstimation ='off'; %% 2020-04-14, AC phase accuracy, QECRX will Large impact the AC phase accuracy
        end
        
        flagR3_EVMcalc = bwInband-LOD.fLO;
        % generate NameCelltitle
        NameCellflag=repmat([{'R3'},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        
        %% 2021-04-19, TEST, flagR3o_TEST_ADC, ADC output, Quantization Error/Dither/NoiseShape
        %% 2021-04-19, flagR3o_TEST_ADC, Add ADC Quantization Error to estimate the QEC ability
        flagR3o_TEST_ADC = 'on';
        if ~strcmpi(flagR3o_TEST_ADC,'off')
            Vref = 1.4;
            nbits = 12;
            LSB = Vref/2^nbits
            flag_RemoveDC = 1
            
            for k = 1:size(waveformCell,1)
                waveformR3o_ADCin = cell2mat(waveformCell(k));
                waveformR3o_ADCout = ADC_Quantizer_Dither(waveformR3o_ADCin, LSB, [], flag_RemoveDC, fs);
                [EVMr3o_ADC, ~] = dsp_evm_timexcorr_inband_g(waveformR3o_ADCin, waveformR3o_ADCout, fs, bwInband, [], 2)
                
                PLOT_FFT_dB_g(waveformR3o_ADCin(:,1), fs, length(waveformT2b_ADCin), ['waveformR3o ADCin'], 'df', 'full', 'pwr', [fnum_Debug], [], []);
                PLOT_FFT_dB_g(waveformR3o_ADCout(:,1), fs, length(waveformR3o_ADCout), ['waveformR3o ADCout'], 'df', 'full', 'pwr', [fnum_Debug], [], []);
                
                waveformCell(k) = {waveformR3o_ADCout};
            end
        end
        
        fnum=fnum+1;
        [waveformR3cell_DnConv,bwInbandR3,EVMr3_DnConv,tableR3_DnConv] = SYM_MixerApp(waveformCell, loD_output, loD_ideal, LOD.fLO, bwInband, fs, flagR3_LODConv, flagR3_EVMcalc, flagR3_RXQECEstimation, bwACLROffsetR3, fnum, NameCelltitle, fnum_dir);
        
        % export
        waveformCell = waveformR3cell_DnConv;
        bwInband = bwInbandR3;
        if ~isempty(bwACLROffsetR3)
            ACLRdBLeft_R3_min = fix(tableR3_DnConv.tableACLRorImage.ACLRdBLeft_min);
            ACLRdBRight_R3_min = fix(tableR3_DnConv.tableACLRorImage.ACLRdBRight_min);
            % excel
            ExcelWrite(ACLRdBLeft_R3_min,'ACLRdBLeft_R3DnConv_min',size(ACLRdBLeft_R3_min),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
            ExcelWrite(ACLRdBRight_R3_min,'ACLRdBRight_R3DnConv_min',size(ACLRdBRight_R3_min),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        end
        
        % excel
        ExcelWrite(round(EVMr3_DnConv,2),'EVMr3_DnConv',size(EVMr3_DnConv),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        if strcmp(flagR3_RXQECEstimation,'RXQEC')
            ExcelWrite(round(tableR3_DnConv.tableQECest.QECest_MagdB_mean,2),'QECest_MagdB_R3',size(tableR3_DnConv.tableQECest.QECest_MagdB_mean),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
            ExcelWrite(round(tableR3_DnConv.tableQECest.QECest_PhsDeg_mean,2),'QECest_PhsDeg_R3',size(tableR3_DnConv.tableQECest.QECest_PhsDeg_mean),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
            ExcelWrite(round(tableR3_DnConv.tableQECest.evmQEC,2),'EVMr3_DnConvQEC',size(tableR3_DnConv.tableQECest.evmQEC),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        end
    end
    
    %% R3b. Apply IFBPF Filter after Mixer or before ADC
    flagR3b_IFBPF = ExcelRead('flagR3b_IFBPF',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagR3b_IFBPF,'off')
        % input
        FIRr3b_Wtype = ExcelRead('FIRr3b_Wtype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRr3b_Ftype = ExcelRead('FIRr3b_Ftype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRr3b_Order = ExcelRead('FIRr3b_Order',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRr3b_fTolerance = ExcelRead('FIRr3b_fTolerance',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRr3b_K_AttdB = ExcelRead('FIRr3b_K_AttdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRr3b_K_fdelta = ExcelRead('FIRr3b_K_fdelta',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRr3b_fcutoffL = ExcelRead('FIRr3b_fcutoffL',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRr3b_fcutoffH = ExcelRead('FIRr3b_fcutoffH',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        
        FIRr3b_IFBPF = SYM_FIRApp(FIRr3b_Wtype,FIRr3b_Ftype,FIRr3b_Order,FIRr3b_K_AttdB,FIRr3b_K_fdelta,FIRr3b_fTolerance,FIRr3b_fcutoffL,FIRr3b_fcutoffH,df,bwInband,NCarriers);
        
        fnum = fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagR3b_IFBPF},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        [waveformR3Acell_IFBPF, b_R3bIFBPF, EVMr3b_IFBPF, tableR3b_IFBPF] = SYM_FilterApp(waveformCell, FIRr3b_IFBPF, fs, mehtod_filter, bwInband, fnum, NameCelltitle, fnum_dir);
        
        % export
        waveformCell= waveformR3Acell_IFBPF;
        
        % excel
        ExcelWrite(tableR3b_IFBPF.tableInput.FIR_Order,'FIRr3b_Order_out',size(tableR3b_IFBPF.tableInput.FIR_Order),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(EVMr3b_IFBPF,2),'EVMr3b_RXBPF',size(EVMr3b_IFBPF),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    %% R4 DownSampling stage R4, to recover original spectrum for ADC sampling
    flagR4_DNSamp = ExcelRead('flagR4_DNSamp',[1],'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagR4_DNSamp,'off')
        % input
        %         fsR4_ADCSampling = ExcelRead('fsR4_ADCSampling',[1],'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        fsR4_ADCSampling = fsT2_UPS
        ratioR4_DNS = fs/fsT2_UPS;
        
        % excel
        ExcelWrite(round(ratioR4_DNS,2),'ratioR4_DNS',size(ratioR4_DNS),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        
        %% R4a. FIRr4DNS
        FIRr4_Wtype = ExcelRead('FIRr4_Wtype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRr4_Ftype = ExcelRead('FIRr4_Ftype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRr4_Order = ExcelRead('FIRr4_Order',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRr4_fTolerance = ExcelRead('FIRr4_fTolerance',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRr4_K_AttdB = ExcelRead('FIRr4_K_AttdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRr4_K_fdelta = ExcelRead('FIRr4_K_fdelta',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRr4_fcutoffL = ExcelRead('FIRr4_fcutoffL',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRr4_fcutoffH = ExcelRead('FIRr4_fcutoffH',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        
        FIRr4_DNS = SYM_FIRApp(FIRr4_Wtype,FIRr4_Ftype,FIRr4_Order,FIRr4_K_AttdB,FIRr4_K_fdelta,FIRr4_fTolerance,FIRr4_fcutoffL,FIRr4_fcutoffH,df,bwInband,NCarriers);
        
        fnum = fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagR4_DNSamp},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        [waveformR4cell_DNS, fsR4_ADCSampling, NsampsR4_DNS, EVMr4_DNS, tableR4_DNS] = SYM_DNSampApp(waveformCell, ratioR4_DNS, fs, [FIRr4_DNS], mehtod_filter, bwInband, [], fnum, NameCelltitle, fnum_dir);
        % export
        waveformCell = waveformR4cell_DNS;
        fs=fsR4_ADCSampling;
        Nsamps=NsampsR4_DNS;
        
        % excel
        ExcelWrite(tableR4_DNS.tableFilterInput.FIR_Order,'FIRr4_Order_out',size(tableR4_DNS.tableFilterInput.FIR_Order),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(EVMr4_DNS,2),'EVMr4_DNS',size(EVMr4_DNS),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(fs/1e6,2),'fsR4_ADCSampling_MHz',size(fs),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(Nsamps,'NsampsR4_DNS',size(fs),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    %% R5. DDC stage R5
    flagR5_DDC = ExcelRead('flagR5_DDC',[1],'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagR5_DDC,'off')
        % input
        %         ratioR5_DDC = ExcelRead('ratioR5_DDC',[1],'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        fsR5_DDC = fsB4_DUC;
        ratioR5_DDC = fs/fsR5_DDC;
        
        %% R5a. FIRr5DDC
        FIRr5_DDC = ExcelRead('FIRr5_DDC',[7],'row','cell',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        %     FIRr5_DDC={b_RHB3,b_RHB2,b_RFIR}; % serier or paralle
        
        fnum = fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagR5_DDC},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        [waveformR5cell_DDC, fsR5_DDC, NsampsR5_DDC, EVMr5_DDC, tableR5_DDC] = SYM_DNSampApp(waveformCell, ratioR5_DDC, fs, [FIRr5_DDC], mehtod_filter, bwInband, [], fnum, NameCelltitle, fnum_dir);
        % export
        waveformCell = waveformR5cell_DDC;
        fs=fsR5_DDC;
        Nsamps=NsampsR5_DDC;
        
        % excel
        ExcelWrite(round(EVMr5_DDC,2),'EVMr5_DDC',size(EVMr5_DDC),'r1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(fs/1e6,2),'fsR5_DDC_MHz',size(fs),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(ratioR5_DDC,2),'ratioR5_DDC',size(ratioR5_DDC),'r1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(Nsamps,'NsampsR5_DDC',size(Nsamps),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    
    %% R6. Add Inband Noise to waveform
    flagR6_addAWGN = ExcelRead('flagR6_addAWGN',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagR6_addAWGN,'off')
        % input
        IpwrdB1Hz_NoiseFloor_R6 = ExcelRead('IpwrdB1Hz_NoiseFloor_R6',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        %     spec_SNRdB_InbandNoise = 50;
        
        fnum=fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagR6_addAWGN},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        [waveformR6cell_AWGN,SNRdBInband_R6,EVMr6_AWGN,tableR6_AWGN] = SYM_AWGNApp(waveformCell, IpwrdB1Hz_NoiseFloor_R6, bwInband, fs, bwInband, fnum, NameCelltitle, fnum_dir);
        
        % export
        waveformCell = waveformR6cell_AWGN;
        
        % excel
        ExcelWrite(round(SNRdBInband_R6,2),'SNRdBInband_R6',size(SNRdBInband_R6),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(EVMr6_AWGN,2),'EVMr6_AWGN',size(EVMr6_AWGN),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    %% R7. Export to waveformRX
    waveformCell_RX = waveformCell;
    fsRX = fs;
    NsampsRX = Nsamps;
    bwInbandRX = bwInband;
    NameCellRX = NameCell;
    clear NsampsDP waveformD7_ACPhsAlign
    
end

%% Digital Signal Process Block =============================================================================================================================
if ~exist('NsampsDP','var')
    disp('Digital Signal Process Block ====================================================')
    %% D0. Receiver Signal import to Digital Processing Block
    waveformCell = waveformCell_RX;
    fs = fsRX;
    Nsamps = NsampsRX;
    bwInband = bwInbandRX;
    NameCell = NameCellRX;
    
    %% D1. NCO shift to Zero and MultiCarrier seperation
    flagD1_NCO2Zero_MultiCarrierDivide = ExcelRead('flagD1_NCO2Zero_MultiCarrierDivide',[1],'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if strcmp(flagD1_NCO2Zero_MultiCarrierDivide,'DIV')
        %         fNCO2Zero = -fNCO;
        fNCO2Zero = ExcelRead('fNCO2Zero',[flagB6_addNCO],'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        fNCO2Zero = -fNCO;
        % excel
        ExcelWrite(fNCO2Zero,'fNCO2Zero',size(fNCO2Zero),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        
        bwInband=bwInbandB6;
        fnum = fnum+1;
        if ~strcmp(flagT8_CombACNbr,'off')&&strcmp(flagC0_AnalysisWF,'DLAC')
            idC=1;
            [waveformD1cell_NCO2Zero_AC, ~, bwInbandD1_AC, ~, tableD1_NCO_AC] = SYM_NCOApp(waveformCell(idC,:), fix(fNCO2Zero(idC,:)/df)*df, fs, bwInband(idC,:), flagD1_NCO2Zero_MultiCarrierDivide, fnum, NameCell(1,:), fnum_dir);
            [waveformD1cell_NCO2Zero_DL, ~, bwInbandD1_DL, ~, tableD1_NCO_DL] = SYM_NCOApp(waveformCell(2:end,:), fix(fNCO2Zero(2:end,:)/df)*df, fs, bwInband(2:end,:), flagD1_NCO2Zero_MultiCarrierDivide, fnum, NameCell(2:end,:), fnum_dir);
            waveformD1cell_NCO2Zero = [waveformD1cell_NCO2Zero_AC;waveformD1cell_NCO2Zero_DL];
            bwInbandD1 = [bwInbandD1_AC;bwInbandD1_DL];
            NameD1Cell=[{'ACComb'};NameC0cell_DL];
            
        else
            [waveformD1cell_NCO2Zero, ~, bwInbandD1, ~, tableD1_NCO] = SYM_NCOApp(waveformCell(:,:), fix(fNCO2Zero(:,:)/df)*df, fs, bwInband(:,:), flagD1_NCO2Zero_MultiCarrierDivide, fnum, NameCell(:,:), fnum_dir);
            NameD1Cell=[NameC0cell];
            
        end
        % export
        waveformCell=waveformD1cell_NCO2Zero;
        NameCell = NameD1Cell;
        NCarriers = size(waveformCell,1);
        bwInband = bwInbandD1;
        
        % excel
        ExcelWrite(NCarriers,'NCarriersD1_AfterNCO2Zero',size(NCarriers),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(NameCell,'CarriersD1_AfterNCO2Zero',size(NameD1Cell),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(bwInband,'bwInbandD1',size(bwInband),'r1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    %% D2. Offset Channel and ChannelPhase Correct
    flagD2_corChOffset = ExcelRead('flagD2_corChOffset',[1],'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    %% 2020-05-07, ISSUE, if NCarriers=1, flagD2_corChannelOffset= 'All'
    if ~strcmp(flagD2_corChOffset,'off')&&~strcmp(flagB5_addChOffset,'off')
        % input
        waveformD2cell_ref = waveformB4cell_DUC;
        
        flagD2_corChOffset= 'CH'; %'CH'/'Phs'/'All'/'off'
        if strcmp(flagC0_AnalysisWF,'DLAC')||strcmp(flagC0_AnalysisWF,'AC')
            waveformD2cell_ref{1,:} = sum(waveformD2cell_ref{1,:},mod(DIMFFT,2)+1);
        end
        [waveformD2cell_OffsetCor, tableD2_OffsetCor, EVMd2_offsetCh,offsetCH_D2est,EVMd2_OffsetChPhs] = SYM_CHOffsetCorApp(waveformCell, waveformD2cell_ref, fs, bwInband, flagD2_corChOffset);
        % export
        waveformCell = waveformD2cell_OffsetCor;
        
        % excel
        ExcelWrite(round(EVMd2_offsetCh,2),'EVMd2_offsetCh',size(EVMd2_offsetCh),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(offsetCH_D2est,'offsetCH_D2est',size(offsetCH_D2est),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    %% D3. DDC stage D3 with FIRd3DDC
    flagD3_DDC = ExcelRead('flagD3_DDC',[1],'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagD3_DDC,'off')
        % input
        ratioD3_DDC = ExcelRead('ratioD3_DDC',[1],'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        ratioD3_DDC = fs/min(fsC0)
        if mod(log2(ratioD3_DDC),1)~=0
            error('ratio should be pwr of 2!')
        end
        
        %% D3a. FIRd3DDC
        FIRd3_Wtype = ExcelRead('FIRd3_Wtype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRd3_Ftype = ExcelRead('FIRd3_Ftype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRd3_Order = ExcelRead('FIRd3_Order',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRd3_fTolerance = ExcelRead('FIRd3_fTolerance',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRd3_K_AttdB = ExcelRead('FIRd3_K_AttdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRd3_K_fdelta = ExcelRead('FIRd3_K_fdelta',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRd3_fcutoffL = ExcelRead('FIRd3_fcutoffL',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRd3_fcutoffH = ExcelRead('FIRd3_fcutoffH',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        
        FIRd3_Wtype = 'Kaiser'
        FIRd3_Ftype = 'LPF'
        FIRd3_fTolerance = 2e6
        FIRd3_K_AttdB = 30
        FIRd3_K_fdelta = 1e6
        FIRd3_fcutoffL = fs/4-FIRd3_fTolerance
        FIRd3_fcutoffL = max(bwInband,[],2)-0
        FIRd3_fcutoffH = 0;
        
        FIRd3_DDC = SYM_FIRApp(FIRd3_Wtype,FIRd3_Ftype,FIRd3_Order,FIRd3_K_AttdB,FIRd3_K_fdelta,FIRd3_fTolerance,FIRd3_fcutoffL,FIRd3_fcutoffH,df,bwInband,NCarriers,fs,[]);
        
        fnum = fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagD3_DDC},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        [waveformD3cell_DDC, fsD3_DDC, NsampsD3_DDC, EVMd3_DDC, tableD3_DDC] = SYM_DNSampApp(waveformCell, ratioD3_DDC, fs, [FIRd3_DDC], mehtod_filter, bwInband, [], fnum, NameCelltitle, fnum_dir);
        % export
        waveformCell = waveformD3cell_DDC;
        fs=fsD3_DDC;
        Nsamps=NsampsD3_DDC;
        
        % excel
        ExcelWrite(tableD3_DDC.tableFilterInput.FIR_Order,'FIRd3_Order_out',size(tableD3_DDC.tableFilterInput.FIR_Order),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(EVMd3_DDC,2),'EVMd3_DDC',size(EVMd3_DDC),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(fs/1e6,2),'fsD3_DDC_MHz',size(fs),[],[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(Nsamps,'NsampsD3_DDC',size(Nsamps),[],[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    %% D4. EQ, Magnitude Ripple Compenstation
    flagD4_corMagRipple = ExcelRead('flagD4_corMagRipple',[1],'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    %     flagD4_corMagRipple = flagB2_addEQRipple;
    flagD4_corMagRipple = 'off'
    if ~strcmp(flagD4_corMagRipple,'off')
        % ref
        waveformD4cell_ref = waveformB1cell_CHF;
        
        fnum = fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagD3_DDC},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        if ~strcmp(flagC0_AnalysisWF,'AC')&&~strcmp(flagC0_AnalysisWF,'DLAC')
            [waveformD4_CorMagRipple, EVMd4_CorMagRipple] = SYM_MagRippleCorApp(waveformCell, waveformD4cell_ref, fs, bwInband, fnum,  NameCelltitle, fnum_dir);
            % export
            waveformCell = waveformD4_CorMagRipple;
            
            % excel
            ExcelWrite(round(EVMd4_CorMagRipple,2),'EVMd4_CorMagRipple',size(EVMd4_CorMagRipple),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
            
        elseif strcmp(flagC0_AnalysisWF,'DLAC')
            %% 2020-05-22, DLAC mode, DL signal is NOT for demodulation, the role is only for Interference to AC signal
            % bypass DL signal
            [waveformD4cell_CorMagRipple_DL, EVMd4_CorMagRipple_DL] = SYM_MagRippleCorApp(waveformCell(2:end,:), waveformD4cell_ref(2:end,:), fs, bwInband(2:end,:), fnum, NameCelltitle(2:end,:), fnum_dir);
            flag_D4corMagRipple_AC =0;
            switch flag_D4corMagRipple_AC
                case{0}
                    waveformCell = [waveformCell(1,:);waveformD4cell_CorMagRipple_DL];
                case{1} % export waveformD4_CorMagRipple_AC
                    [waveformD4_CorMagRipple_AC, EVMd4_CorMagRipple_AC] = SYM_MagRippleCorApp(waveformCell(1), sum(waveformD4cell_ref{1,:},mod(DIMFFT,2)+1), fs, bwInband(1,:), fnum, fnum_dir);
                    waveformCell = [waveformD4_CorMagRipple_AC;waveformD4cell_CorMagRipple_DL];
                    % excel
                    ExcelWrite(round(waveformD4_CorMagRipple_AC,2),'waveformD4_CorMagRipple_AC',size(waveformD4_CorMagRipple_AC),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
            end
            
            % excel
            ExcelWrite(round(EVMd4_CorMagRipple_DL,2),'EVMd4_CorMagRipple_DL',size(EVMd4_CorMagRipple_DL),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,[0 1]+xlsWriteShift);
        end
    end
    
    %% D5. Channel Filter for UL, FIRd5CHFUL
    flagD5_CHFUL = ExcelRead('flagD5_CHFUL',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    
    if ~strcmp(flagD5_CHFUL,'off')
        % input
        FIRd5_Wtype = ExcelRead('FIRd5_Wtype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRd5_Ftype = ExcelRead('FIRd5_Ftype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRd5_Order = ExcelRead('FIRd5_Order',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRd5_fTolerance = ExcelRead('FIRd5_fTolerance',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRd5_K_AttdB = ExcelRead('FIRd5_K_AttdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRd5_K_fdelta = ExcelRead('FIRd5_K_fdelta',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRd5_fcutoffL = ExcelRead('FIRd5_fcutoffL',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        FIRd5_fcutoffH = ExcelRead('FIRd5_fcutoffH',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        
        FIRd5_Wtype = 'Kaiser'
        FIRd5_Ftype = 'LPF'
        FIRd5_fTolerance = 1.0e6
        FIRd5_K_AttdB = 20
        FIRd5_K_fdelta = 0.1e6
        FIRd5_fcutoffL = max(bwInband,[],2)
        FIRd5_fcutoffH = 0
        %         FIRd5_Wtype = 'Remez'
        %         FIRd5_Ftype = 'LPF'
        %         FIRd5_K_AttdB = [0.1 40]
        [FIRd5_CHFUL] = SYM_FIRApp(FIRd5_Wtype,FIRd5_Ftype,FIRd5_Order,FIRd5_K_AttdB,FIRd5_K_fdelta,FIRd5_fTolerance,FIRd5_fcutoffL,FIRd5_fcutoffH,df,bwInband,NCarriers,fs,[]);
        
        fnum=fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagD5_CHFUL},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        [waveformD5cell_CHFUL, b_D5CHFUL, EVMd5_CHFUL, tableD5_CHFUL] = SYM_FilterApp(waveformCell, FIRd5_CHFUL, fs, mehtod_filter, bwInband, fnum, NameCelltitle, fnum_dir);
        % export
        waveformCell= waveformD5cell_CHFUL;
        
        % excel
        ExcelWrite(tableD5_CHFUL.tableInput.FIR_Order,'FIRd5_Order_out',size(tableD5_CHFUL.tableInput.FIR_Order),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(EVMd5_CHFUL,2),'EVMd5_CHFUL',size(EVMd5_CHFUL),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    %% D6. Downsampling stage to Demodulation
    flagD6_DNSamp = ExcelRead('flagD6_DNSamp',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
    ratioD6_DNSamp = fs./fsC0;
    
    if strcmp(flagD6_DNSamp,'D6DNS')&&(~any(ratioD6_DNSamp))
        % input
        method_FIRd6_DNS = 2;
        
        switch method_FIRd6_DNS
            case {1}
                FIRd6_DNS = ExcelRead('FIRd6_DNS',[8],'column','cell',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
            case {2}
                FIRd6_Wtype = ExcelRead('FIRd6_Wtype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
                FIRd6_Ftype = ExcelRead('FIRd6_Ftype',1,'column','char',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
                FIRd6_Order = ExcelRead('FIRd6_Order',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
                FIRd6_fTolerance = ExcelRead('FIRd6_fTolerance',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
                FIRd6_K_AttdB = ExcelRead('FIRd6_K_AttdB',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
                FIRd6_K_fdelta = ExcelRead('FIRd6_K_fdelta',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
                FIRd6_fcutoffL = ExcelRead('FIRd6_fcutoffL',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
                FIRd6_fcutoffH = ExcelRead('FIRd6_fcutoffH',1,'row','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
                
                FIRd6_Ftype = 'LPF'
                FIRd6_fTolerance = 1e6
                FIRd6_K_AttdB = 30
                FIRd6_K_fdelta = 1e6
                FIRd6_fcutoffL = max(bwInband)+FIRd6_fTolerance
                FIRd6_fcutoffH = 0;
                
                [FIRd6_DNS] = SYM_FIRApp(FIRd6_Wtype,FIRd6_Ftype,FIRd6_Order,FIRd6_K_AttdB,FIRd6_K_fdelta,FIRd6_fTolerance,FIRd6_fcutoffL,FIRd6_fcutoffH,df,bwInband,NCarriers,fs,[]);
                
        end
        
        fnum=fnum+1;
        % generate NameCelltitle
        NameCellflag=repmat([{flagD6_DNSamp},{', '}],size(waveformCell,1),1);
        NameCelltitle=[NameCellflag,NameCell];
        [waveformD6cell_DNS, fsD6_DNS, NsampsD6_DNS, EVMd6_DNS, tableD6_DNS] = SYM_DNSampApp(waveformCell, ratioD6_DNSamp, fs, [FIRd6_DNS], mehtod_filter, bwInband, [], fnum, NameCelltitle, fnum_dir);
        % export
        waveformCell = waveformD6cell_DNS;
        fs=fsD6_DNS;
        Nsamps=NsampsD6_DNS;
        
        % excel
        if exist('tableD6_DNS.tableFilterInput.FIR_Order','var')
            ExcelWrite(tableD6_DNS.tableFilterInput.FIR_Order,'FIRd6_Order',size(tableD6_DNS.tableFilterInput.FIR_Order),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        end
        ExcelWrite(round(EVMd6_DNS,2),'EVMd6_DNS',size(EVMd6_DNS),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(fs/1e6,2),'fsD6_DNS_MHz',size(fs),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(Nsamps,'NsampsD6_DNS',size(Nsamps),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
    end
    
    %% D7. Channel-Phase Correct after Channel-Filter application and before Demodulation
    if (~strcmp(flagD2_corChOffset,'All')||~strcmp(flagD2_corChOffset,'Phs'))&&~strcmp(flagC0_AnalysisWF,'AC') %'CH'/'Phs'/'All'/'off'
        flagD7_corChPhsOffset = 'Phs';
        waveformD7cell_ref = waveformC0cell_DL;
        
        if strcmp(flagC0_AnalysisWF,'DLAC') % waveform 1st is AC
            idC=2;
        elseif ~strcmp(flagC0_AnalysisWF,'AC')&&~strcmp(flagC0_AnalysisWF,'DLAC')
            idC=1;
        end
        [waveformD7cell_OffsetPhsCor, tableD7_OffsetPhsCor, EVMd7_offsetCh,offsetCH_D7Est,EVMd7_OffsetChPhs] = SYM_CHOffsetCorApp(waveformCell(idC:end,:), waveformD7cell_ref, fs, bwInband(idC:end,:), flagD7_corChPhsOffset);
        
        % export
        waveformD7cell_DL=waveformD7cell_OffsetPhsCor;
    end
    
    %% D7. OFDM DeModulation
    if exist('waveformD7cell_DL','var')
        fnum=fnum+1;
        [rxGridD7cell, EVMd7_DemodulationDL] = SYM_DeMODApp(waveformD7cell_DL, waveformD7cell_ref, DLconfig, fs, bwInbandC0, fnum, fnum_dir);
        if strcmp(flagC0_AnalysisWF,'DLAC')
            xlsWriteShift2 = [0 1];
        else
            xlsWriteShift2 = [0 0];
        end
        ExcelWrite(round(EVMd7_DemodulationDL,2),'EVMd7_DemodulationDL',[size(EVMd7_DemodulationDL)],'r1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift+xlsWriteShift2);
    end
    
    
    %% D7b. AC DeModulation
    if strcmp(flagC0_AnalysisWF,'DLAC')||strcmp(flagC0_AnalysisWF,'AC')
        ACDomod_Nshift = 0;
        
        fnum=fnum+1;
        [ACdmod_t0Mean, ACdmod_p0DegMean, ACdmod_phEstDeg, ACdmod_phEstDegDrift, ACdmod_dataCapCor, ACdmod_SNR, ACdmod_dataCapNbrwoPD] = AntCal_phaseDemodulateApp_g100(waveformCell{1,:},ACconfig,ACDomod_Nshift,fnum, fnum_dir);
        ACPhsEstDegError = ACdmod_p0DegMean-ACPhsShiftNbrDeg_T7;
        
        % excel
        ACPhsEstDegMeanErrorMax = max(abs(ACPhsEstDegError));
        ACPhsEstDegDriftMax = max(ACdmod_phEstDegDrift);
        ACdmod_SNRdBmin = min(ACdmod_SNR);
        
        ExcelWrite(round(ACPhsEstDegMeanErrorMax,2),'ACPhsEstDegMeanErrorMax',size(ACPhsEstDegMeanErrorMax),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(ACPhsEstDegDriftMax,2),'ACPhsEstDegDriftMax',size(ACPhsEstDegDriftMax),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        ExcelWrite(round(ACdmod_SNRdBmin,2),'ACdmod_SNRdBmin',size(ACdmod_SNRdBmin),'1row',[],xlsFile,xlsSheetName,xlsRangeArrayInput,xlsWriteShift);
        
        %% 2020-06-09, TEST, flag_TEST_Reduce_ACnSubcarrier, reduce ACconfig.nSubcarrier
        flag_TEST_Reduce_ACnSubcarrier = 0;
        if flag_TEST_Reduce_ACnSubcarrier==1
            flagD71_ACnSubcarCut_TEST=30;
            ACconfig_TEST = ACconfig;
            ACconfig_TEST.nSubcarrier =ACconfig_TEST.nSubcarrier(1+flagD71_ACnSubcarCut_TEST/2:end-flagD71_ACnSubcarCut_TEST/2);
            
            fnum=fnum+1;
            [ACdmod2_t0Mean, ACdmod2_p0DegMean, ACdmod2_phEstDeg, ACdmod2_phEstDegDrift2, ACdmod2_dataCapCor, ACdmod2_SNR, ACdmod2_dataCapNbrwoPD] = AntCal_phaseDemodulateApp_g100(waveformCell{1,:},ACconfig_TEST,ACDomod_Nshift,fnum, fnum_dir);
        end
        
        %% D7c. AC Calculate EQ coefficents and Phase Alignment
        ACEQ_numFirTaps = ExcelRead('ACEQ_numFirTaps',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        ACEQ_NoOfAlignBr = ExcelRead('ACEQ_NoOfAlignBr',1,'column','scalar',xlsFile,xlsSheetName,xlsRangeArrayInput,xlsReadShift);
        
        fnum=fnum+1;
        [waveformD7_ACPhsAlign,ACEQ_taps,ACEQ_B] = AntCal_EQphsAlignApp(waveformC0_AC,ACconfig,[], ACdmod_phEstDeg,ACEQ_numFirTaps,ACEQ_NoOfAlignBr,fnum, fnum_dir);
        ACPhsShiftNbrDeg_T7_Ref = (ACPhsShiftNbrDeg_T7(ACEQ_NoOfAlignBr)-ACPhsShiftNbrDeg_T7).'.*ones(1,length(ACconfig.nSubcarrier));
        for idBR = 1:size(ACPhsShiftNbrDeg_T7_Ref,1)
            plot(ACconfig.nSubcarrier, ACPhsShiftNbrDeg_T7_Ref(idBR,:),'r','DisplayName',['Ref idBR',num2str(idBR),', [',num2str(mean(ACPhsShiftNbrDeg_T7_Ref(idBR,:))),']deg'])
            legend;
        end
        
    end
    
    %% D7c. 2021-03-10, TEST, flagD7c_TEST_ACBrGainVariation: waveformC0_AC with Gain variations for each Branch
    flagD7c_TEST_ACBrGainVariation = 'off';
    if strcmpi(flagD7c_TEST_ACBrGainVariation,'on')
        % Gain vs Branch
        GaindB_D7c = 40*[0, 1, -1, 0, 1, 0, -1, 0];
        GaindB_D7c = 60*[0, 1, -1, 0, 1, 0, -1, 0];
        GaindB_D7c = 80*[0, 1, -1, 0, 1, 0, -1, 0];
        gb
        % Add Gain variations
        waveformD7c_AC_GainVsBr = waveformC0_AC.*10.^(GaindB_D7c/10);
        PLOT_FFT_dB_g(waveformD7c_AC_GainVsBr, fsAC, NsampsC0_AC*2^4, 'waveformD7c_AC_GainVsBr', 'df', 'full', 'pwr', 20210401);
        
        % Add Phase shifts
        waveformD7c_AC_GainVsBr_PhShift = exp(1i*ACPhsShiftNbrDeg_T7./180*pi).*waveformD7c_AC_GainVsBr;
        
        % Sum of AC waveform
        waveformD7c_AC_GainVsBr_PhShift_Sum = sum(waveformD7c_AC_GainVsBr_PhShift,2);
        PLOT_FFT_dB_g(waveformD7c_AC_GainVsBr_PhShift_Sum, fsAC, NsampsC0_AC*2^4, 'waveformD7c_AC_GainVsBr_PhShift_Sum', 'df', 'full', 'pwr', 20210401);
        
        % AC DeModulation
        fnum=fnum+1;
        [ACdmod_t0Mean_D7c, ACdmod_p0DegMean_D7c, ACdmod_phEstDeg_D7c, ACdmod_phEstDegDrift_D7c, ACdmod_dataCapCor_D7c, ACdmod_SNR_D7c, ACdmod_dataCapNbrwoPD_D7c] = ...
            AntCal_phaseDemodulateApp_g100(waveformD7c_AC_GainVsBr_PhShift_Sum,ACconfig,ACDomod_Nshift,fnum, fnum_dir);
        
        % AC EQ
        [waveformD7_ACPhsAlign_D7c,ACEQ_taps_D7c,ACEQ_B_D7c] = AntCal_EQphsAlignApp(waveformC0_AC,ACconfig,[], ACdmod_phEstDeg,ACEQ_numFirTaps,ACEQ_NoOfAlignBr,fnum, fnum_dir);
        
    end
    
    %% D.8. Export NsampsDP for ACPhsAlign analysis
    
end



