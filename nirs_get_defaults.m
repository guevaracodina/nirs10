function varargout = nirs_get_defaults(defstr, varargin)
% Get/set the defaults values associated with an identifier
% FORMAT defval = nirs_get_defaults(defstr)
% Return the defaults value associated with identifier "defstr". 
% Currently, this is a '.' subscript reference into the global  
% "defaults" variable defined in spm_defaults.m.
%
% FORMAT nirs_get_defaults(defstr, defval)
% Sets the nirs10 value associated with identifier "defstr". The new
% vbm8 value applies immediately to:
% * new modules in batch jobs
% * modules in batch jobs that have not been saved yet
% This value will not be saved for future sessions of SPM. To make
% persistent changes, edit cg_vbm8_defaults.m.
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
%
% based on Christian Gaser version of
% cg_vbm8_get_defaults $Id: cg_vbm8_get_defaults.m 2696 2009-02-05
% 20:29:48Z guillaume $
%
% Clément Bonnéry

global nirs10;
if isempty(nirs10)
    %load various default files
    %nirs_defaults;
    nirs_boxy_defaults;
    nirs_criugm_defaults;
    nirs_paces_defaults;
    nirs_testStimuli_defaults;
    nirs_buildroi_defaults;
    nirs_MCsegment_defaults;
    nirs_coreg_defaults;
    nirs_view3d_defaults;
    nirs_resize_defaults;
    nirs_remove_chn_stdev_defaults;
    nirs_mark_negative_defaults;
    nirs_mark_movement_defaults;
    nirs_normalize_baseline_defaults;
    nirs_ODtoHbOHbR_defaults;
    nirs_configMC_defaults;
    %nirs_MC_defaults; %File missing?
    nirs_NIRS_SPM_specify_defaults;
    nirs_NIRS_SPM_estimate_defaults;
    nirs_LIOM_estimate_defaults;
    nirs_LIOM_specify_defaults;  
end

% construct subscript reference struct from dot delimited tag string
tags = textscan(defstr,'%s', 'delimiter','.');
subs = struct('type','.','subs',tags{1}');
%Added to ensure that default value of first entry of first module loaded
%is found; traced problem to inconsistency/unpredictability in whether 
%nirs10 was included or not in subs, and it should always be excluded
if strcmp(subs(1).subs, 'nirs10')
    subs = subs(2:end);
end
   
if nargin == 1
    varargout{1} = subsref(nirs10, subs);
else
    nirs10 = subsasgn(nirs10, subs, varargin{1});
end
