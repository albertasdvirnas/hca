configCluster

% Set ProjectName and WallTime before submitting jobs to AURORA
c = parcluster;
c.AdditionalProperties.ProjectName = 'project-name';
c.AdditionalProperties.WallTime = '01:00:00';
c.saveProfile