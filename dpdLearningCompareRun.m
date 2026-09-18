function [results, summaryTable] = dpdLearningCompareRun( cfg, saveModels)

    if nargin < 1 || isempty(cfg)
        cfg = config();
    end
    if nargin < 2 || isempty(saveModels)
        saveModels = true;
    end

    validateInputs(cfg, saveModels);
    saveModels = logical(saveModels);

    indirectRunCfg = architectureConfiguration( cfg, 'indirect', saveModels);
    directRunCfg = architectureConfiguration( cfg, 'direct', saveModels);

    fprintf('\nDPD direct/indirect learning comparison.\n');
    fprintf(['Common training data: 300000 samples in 10 blocks of ' '30000 samples.\n']);
    fprintf('Common train/validation/test seeds: %d/%d/%d.\n', indirectRunCfg.dpd.trainingSeed, indirectRunCfg.dpd.validationSeed, indirectRunCfg.dpd.testSeed);

    [indirectDpdCfg, indirectTrainingInfo] = dpdTrainRun(indirectRunCfg);
    [directDpdCfg, directTrainingInfo] = dpdTrainRun(directRunCfg);

    assertComparableRuns(indirectTrainingInfo, directTrainingInfo);

    results.indirect.dpdConfig = indirectDpdCfg;
    results.indirect.trainingInfo = indirectTrainingInfo;
    results.direct.dpdConfig = directDpdCfg;
    results.direct.trainingInfo = directTrainingInfo;
    results.common.trainingSampleCount = 300000;
    results.common.blockLength = 30000;
    results.common.blockCount = 10;
    results.common.trainingSeed = indirectTrainingInfo.options.trainingSeed;
    results.common.validationSeed = indirectTrainingInfo.options.validationSeed;
    results.common.testSeed = indirectTrainingInfo.options.testSeed;
    results.common.referenceGain = indirectTrainingInfo.referenceGain;
    results.common.modelsSaved = saveModels;

    summaryRows(1) = summaryRow( 'indirect', indirectTrainingInfo);
    summaryRows(2) = summaryRow( 'direct', directTrainingInfo);
    summaryTable = summaryRowsToTable(summaryRows);
    results.summary = summaryRows;
    results.summaryTable = summaryTable;

    fprintf('\nFixed-G validation/test comparison:\n');
    disp(summaryTable);
end

function runCfg = architectureConfiguration( baseCfg, learningArchitecture, saveModels)
    runCfg = baseCfg;
    cfrCfg = systemCfrConfiguration(baseCfg);
    if ~isfield(baseCfg, 'dpd') || isempty(baseCfg.dpd)
        runCfg.dpd = dpdConfig(baseCfg.pa, false, learningArchitecture, cfrCfg);
    else
        runCfg.dpd = baseCfg.dpd;
    end

    runCfg.dpd.learningArchitecture = learningArchitecture;
    [runCfg.dpd, ~] = dpdRefreshCfrCompatibility(runCfg.dpd, cfrCfg);
    runCfg.dpd.blockLength = 30000;
    runCfg.dpd.adaptationBlockLength = 30000;
    runCfg.dpd.trainingSampleCount = 300000;
    runCfg.dpd.saveModel = saveModels;
    runCfg.dpd.modelLoaded = false;
    runCfg.dpd.learningInitialization = 'identity';

    runCfg.dpd.modelFile = dpdArchitectureArtifactFile(learningArchitecture, runCfg.dpd.cfrSignature, 'model');
end

function row = summaryRow(architecture, trainingInfo)
    validation = trainingInfo.validation;
    test = trainingInfo.test;
    options = trainingInfo.options;

    row.Architecture = architecture;
    row.BlockCount = ceil( options.trainingSampleCount / options.blockLength);
    row.ValidationNmseDb = validation.candidate.nmseDb;
    row.ValidationImprovementDb = validation.improvementDb;
    row.ValidationAclr1Db = validation.candidate.worstAclr1Db;
    row.ValidationAclr2Db = validation.candidate.worstAclr2Db;
    row.ValidationDpdPowerDb = validation.candidate.dpdAveragePowerChangeDb;
    row.TestNmseDb = test.candidate.nmseDb;
    row.TestImprovementDb = test.improvementDb;
    row.TestAclr1Db = test.candidate.worstAclr1Db;
    row.TestAclr2Db = test.candidate.worstAclr2Db;
    row.TestDpdPowerDb = test.candidate.dpdAveragePowerChangeDb;
end

function summaryTable = summaryRowsToTable(rows)
    summaryTable = table( {rows.Architecture}.', [rows.BlockCount].', [rows.ValidationNmseDb].', [rows.ValidationImprovementDb].', [rows.ValidationAclr1Db].', [rows.ValidationAclr2Db].', [rows.ValidationDpdPowerDb].', [rows.TestNmseDb].', [rows.TestImprovementDb].', [rows.TestAclr1Db].', [rows.TestAclr2Db].', [rows.TestDpdPowerDb].', 'VariableNames', { 'Architecture', 'BlockCount', 'ValidationNmseDb', 'ValidationImprovementDb', 'ValidationAclr1Db', 'ValidationAclr2Db', 'ValidationDpdPowerDb', 'TestNmseDb', 'TestImprovementDb', 'TestAclr1Db', 'TestAclr2Db', 'TestDpdPowerDb'});
end

function assertComparableRuns(indirectInfo, directInfo)
    commonOptionFields = {'trainingSeed', 'validationSeed', 'testSeed', 'trainingSampleCount', 'validationSampleCount', 'testSampleCount', 'blockLength'};
    for fieldIndex = 1:length(commonOptionFields)
        fieldName = commonOptionFields{fieldIndex};
        if indirectInfo.options.(fieldName) ~= directInfo.options.(fieldName)
            error('dpdLearningCompareRun:IncomparableOptions', ['Direct and indirect learning used different %s ' 'values.'], fieldName);
        end
    end

    assertClose(indirectInfo.referenceGain, directInfo.referenceGain, 'reference gain');
    assertClose(indirectInfo.validation.baseline.nmseDb, directInfo.validation.baseline.nmseDb, 'validation PA baseline NMSE');
    assertClose(indirectInfo.test.baseline.nmseDb, directInfo.test.baseline.nmseDb, 'test PA baseline NMSE');
end

function assertClose(firstValue, secondValue, description)
    scale = max(1, max(abs(firstValue), abs(secondValue)));
    tolerance = 1e-10 * scale;
    if ~isfinite(firstValue) || ~isfinite(secondValue) || abs(firstValue - secondValue) > tolerance
        error('dpdLearningCompareRun:IncomparableData', 'Direct and indirect learning did not use the same %s.', description);
    end
end

function validateInputs(cfg, saveModels)
    if ~isstruct(cfg) || ~isfield(cfg, 'pa') || ~isstruct(cfg.pa)
        error('dpdLearningCompareRun:MissingPaConfiguration', 'cfg.pa is required.');
    end
    if ~(islogical(saveModels) || isnumeric(saveModels)) || ~isscalar(saveModels) || ~isfinite(saveModels) || (saveModels ~= 0 && saveModels ~= 1)
        error('dpdLearningCompareRun:InvalidSaveModels', 'saveModels must be a logical scalar.');
    end
end

function cfrCfg = systemCfrConfiguration(cfg)
    if isfield(cfg, 'cfr') && isstruct(cfg.cfr) && isscalar(cfg.cfr)
        cfrCfg = cfg.cfr;
    else
        cfrCfg = cfrConfig(false);
    end
end
