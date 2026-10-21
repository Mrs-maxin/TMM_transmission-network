function dependencies = setup_external_dependencies(config)
%SETUP_EXTERNAL_DEPENDENCIES Add MATPOWER and the AQPPF_TMM source to the path.
% AQPPF_TMM is used unchanged as the published distribution-network TMM
% implementation; its numerical code must match the registered commit.
arguments
    config (1,1) struct
end
assert(isfolder(config.matpowerRoot), "TMM_Transmission:dependency:MatpowerMissing", ...
    "MATPOWER folder is missing: %s. Set config.matpowerRoot.", config.matpowerRoot);
assert(isfolder(fullfile(config.aqppfRoot, "+TMM_core")), "TMM_Transmission:dependency:AqppfMissing", ...
    "AQPPF_TMM source is missing: %s. Set config.aqppfRoot.", config.aqppfRoot);
addpath(genpath(char(config.matpowerRoot)));
addpath(char(config.aqppfRoot));
required = ["runpf", "makeYbus", "dSbus_dV", "kmeans", "tcdf"];
for name = required
    assert(exist(name, "file") > 0, "TMM_Transmission:dependency:FunctionMissing", ...
        "Required function %s is not visible on the MATLAB path.", name);
end
assert(~isempty(which("TMM_core.mixture.fit_order_candidates")), ...
    "TMM_Transmission:dependency:AqppfPackage", ...
    "TMM_core.mixture.fit_order_candidates is not visible under %s.", config.aqppfRoot);

[status, commit] = system(sprintf('git -C "%s" rev-parse HEAD', config.aqppfRoot));
commit = strtrim(string(commit));
if status ~= 0
    commit = "unavailable";
end
matpowerVersion = mpver;
dependencies = struct("matpowerRoot", string(config.matpowerRoot), ...
    "matpowerVersion", string(matpowerVersion), ...
    "aqppfRoot", string(config.aqppfRoot), "aqppfCommit", commit, ...
    "aqppfNumericalCodeEquivalentTo", string(config.aqppfReferenceCommit), ...
    "matlabVersion", string(version));
end
