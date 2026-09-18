function [dpdCfg, compatible] = dpdRefreshCfrCompatibility(dpdCfg, cfrCfg)
    if nargin < 1 || ~isstruct(dpdCfg) || ~isscalar(dpdCfg)
        error('dpdRefreshCfrCompatibility:InvalidDpdConfiguration', 'dpdCfg must be a scalar structure.');
    end
    if nargin < 2 || isempty(cfrCfg)
        cfrCfg = cfrConfig(false);
    end

    currentSignature = cfrSignature(cfrCfg);
    if isfield(dpdCfg, 'cfrSignature') && ~isempty(dpdCfg.cfrSignature)
        previousSignature = dpdCfg.cfrSignature;
    else
        previousSignature = cfrSignature(cfrConfig(false));
    end
    compatible = isequaln(previousSignature, currentSignature);

    if ~compatible
        dpdCfg.modelLoaded = false;
        dpdCfg.loadedModelArchitecture = '';
        dpdCfg.modelStale = true;
        dpdCfg.modelStaleReason = 'The DPD coefficients belong to a different CFR configuration and must be retrained.';
        dpdCfg.trained = false;
        dpdCfg.modelVerified = false;
        if isfield(dpdCfg, 'savedTrainingInfo')
            dpdCfg = rmfield(dpdCfg, 'savedTrainingInfo');
        end
    end

    dpdCfg.cfrSignature = currentSignature;
    if isfield(dpdCfg, 'learningArchitecture') && ischar(dpdCfg.learningArchitecture)
        if ~isfield(dpdCfg, 'modelFile') || isDefaultArtifactFile(dpdCfg.modelFile)
            dpdCfg.modelFile = dpdArchitectureArtifactFile(dpdCfg.learningArchitecture, currentSignature, 'model');
        end
    end
end

function matched = isDefaultArtifactFile(modelFile)
    matched = false;
    if ~ischar(modelFile)
        return
    end
    [~, modelName, modelExtension] = fileparts(modelFile);
    defaultNames = {'dpdModel', 'dpdModelIndirect', 'dpdModelDirect', 'dpdModelIndirectCfr', 'dpdModelDirectCfr'};
    matched = strcmpi(modelExtension, '.mat') && any(strcmpi(modelName, defaultNames));
end
