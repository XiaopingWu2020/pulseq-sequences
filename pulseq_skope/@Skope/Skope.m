%% PulSeq interface for Skope
% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

classdef Skope < handle
%%%%%%%%%%%%%%
% Basic class to use PulSeq with Skope® systems
%%%%%%%%%%%%%%
    properties
        seq_name;
        seq_type;
        scanner;
        seq_params;
        
        % PulSeq objects
        sys;
        seq;

        %
        isImagingSeq = false;
    end
    
    methods
        % Constructor
        function obj = Skope(seq_name, seq_type, scanner_type)
            obj.seq_name = seq_name;
            obj.seq_type = seq_type;
            [obj.sys, obj.scanner] = obj.Set_scanner(scanner_type);
            obj.seq = mr.Sequence();
        end

        TestReport(obj);
        SkopeReport(obj);
    end

    methods (Access = private)
        Prepare_2dGre(obj);
        Prepare_2dEpi(obj);
        Prepare_2dSpiral(obj);
        Prepare_girf_calib(obj);
        Prepare_grad_tones(obj);
        Prepare_2dEpi_moco(obj);
        Prepare_2dSpiral_moco(obj);
        Prepare_2dSpiral_stitch(obj);
        Prepare_2dSpiral_stitch_multishot(obj);
        delayTR= calculateTiming(obj, minExcitTR);
        [gx,gy,gz,trigDelay]= createGradientTones(obj, A0,f0,tw,basegrad);
        grad= createTrapGrad(obj,basegrad,sf);
        gos= designArchimedean(obj,deltak,kRadius,kSamples,safety_margin);
        [k,g,s,time] = designVariableDensitySpiral(obj, Nitlv, rv, res, fov, radius, safetyMargin);
    end

end
