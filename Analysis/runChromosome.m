function chrom = runChromosome(saveData)
%RUNCHROMOSOME  Best chromosome of a saved run, in rotation-vector genes.
%
%   chrom = runChromosome(saveData) accepts the saveData struct written by
%   saveResults. Runs saved before the exponential-map change store XYZ
%   Euler genes (Specifications.Parameterisation missing or 'euler'); those
%   are converted so every caller can pass the result to setupCameras.
    chrom = saveData.BestSolution.Chromosome(:).';
    specs = saveData.Specifications;
    if ~isfield(specs, 'Parameterisation') || strcmpi(specs.Parameterisation, 'euler')
        chrom = eulerChromToRotvec(chrom);
    end
end
