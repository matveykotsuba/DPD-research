function artifactFile = dpdArchitectureArtifactFile(learningArchitecture, cfrCfgOrSignature, artifactType)
    if nargin < 2 || isempty(cfrCfgOrSignature)
        cfrCfgOrSignature = cfrSignature(cfrConfig(false));
    end
    if nargin < 3 || isempty(artifactType)
        artifactType = 'model';
    end
    if ~ischar(learningArchitecture)
        error('dpdArchitectureArtifactFile:InvalidArchitecture', 'learningArchitecture must be a character vector.');
    end
    if ~isstruct(cfrCfgOrSignature) || ~isscalar(cfrCfgOrSignature) || ~isfield(cfrCfgOrSignature, 'enabled')
        error('dpdArchitectureArtifactFile:InvalidCfrConfiguration', 'The CFR configuration or signature must contain enabled.');
    end
    if ~ischar(artifactType)
        error('dpdArchitectureArtifactFile:InvalidArtifactType', 'artifactType must be a character vector.');
    end

    learningArchitecture = lower(strtrim(learningArchitecture));
    if ~strcmp(learningArchitecture, 'direct') && ~strcmp(learningArchitecture, 'indirect')
        error('dpdArchitectureArtifactFile:UnknownArchitecture', 'learningArchitecture must be direct or indirect.');
    end
    enabled = cfrCfgOrSignature.enabled;
    if ~(islogical(enabled) || isnumeric(enabled)) || ~isscalar(enabled) || ~isreal(enabled) || ~isfinite(enabled) || (enabled ~= 0 && enabled ~= 1)
        error('dpdArchitectureArtifactFile:InvalidCfrEnabled', 'The CFR enabled value must be a logical scalar.');
    end
    if logical(enabled)
        cfrSuffix = 'Cfr';
    else
        cfrSuffix = '';
    end

    architectureName = [upper(learningArchitecture(1)) learningArchitecture(2:end)];
    artifactType = lower(strtrim(artifactType));
    switch artifactType
        case 'model'
            artifactName = ['dpdModel' architectureName cfrSuffix '.mat'];
        case 'comparison'
            artifactName = ['dpdComparison' architectureName cfrSuffix '.mat'];
        otherwise
            error('dpdArchitectureArtifactFile:UnknownArtifactType', 'artifactType must be model or comparison.');
    end
    artifactFile = fullfile(fileparts(mfilename('fullpath')), artifactName);
end
