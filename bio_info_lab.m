% BioInfoLab - Advanced Sequence Analysis Toolbox

% Load sequence data
sequence = 'ACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCT';

% Calculate sequence length
sequenceLength = length(sequence);

% Count nucleotide occurrences
nucleotideCounts = [sum(sequence == 'A'), sum(sequence == 'C'), sum(sequence == 'G'), sum(sequence == 'T')];

% Calculate GC content
gcContent = (nucleotideCounts(2) + nucleotideCounts(3)) / sequenceLength * 100;

% Find sequence motifs
motif = 'CTAG';
motifOccurrences = strfind(sequence, motif);
numMotifs = numel(motifOccurrences);

% Perform sequence alignment
sequence1 = 'ACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCT';
sequence2 = 'ACGTAGCTAGCTAGCTAGCTAGCTAGCTCGCT';
alignmentScore = nwalign(sequence1, sequence2);

% Generate reverse complement sequence
reverseComplement = seqrcomplement(sequence);

% Calculate melting temperature using nearest-neighbor method
meltingTemp = nnntemp(sequence);

% Perform motif enrichment analysis
motifSet = {'CTAG', 'GCTA', 'TAGC'};
enrichmentPValues = motifenrichment(sequence, motifSet);

% Search for open reading frames (ORFs)
minORFLength = 50;
orfStartCodons = {'ATG'};
orfStopCodons = {'TAA', 'TAG', 'TGA'};
orfPositions = findorfs(sequence, minORFLength, orfStartCodons, orfStopCodons);

% Predict secondary structure using RNAfold algorithm
rnaSecondaryStructure = RNAfold(sequence);

% Perform codon usage analysis
codonUsageTable = codonusage(sequence);

% Calculate protein molecular weight
proteinSequence = 'MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPNEECLFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV';
proteinMolecularWeight = seqmass(proteinSequence, 'protein');

% Additional Features:

% Perform sequence motif search
motifSearch = 'AGCT';
motifSearchOccurrences = strfind(sequence, motifSearch);
numMotifSearch = numel(motifSearchOccurrences);

% Calculate sequence similarity
similarityScore = seqpdist({sequence, sequence1}, 'Method', 'Jukes-Cantor');

% Perform protein structure prediction using homology modeling
proteinModel = proteinmodel(proteinSequence);

% Analyze protein-protein interactions
proteinInteractions = proteininteractions(proteinSequence);

% Perform multiple sequence alignment
sequences = {'ACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCT', 'ACGTAGCTAGCTAGCTAGCTAGCTAGCTCGCT', 'ACGTAGCTAGCTAGCTAGCTAGCTAGCTAGTT'};
multipleAlignment = multialign(sequences);

% Display results
disp(['Sequence Length: ', num2str(sequenceLength)]);
disp(['Nucleotide Counts: A=', num2str(nucleotideCounts(1)), ' C=', num2str(nucleotideCounts(2)), ' G=', num2str(nucleotideCounts(3)), ' T=', num2str(nucleotideCounts(4))]);
disp(['GC Content: ', num2str(gcContent), '%']);
disp(['Motif Occurrences: ', num2str(numMotifs)]);
disp(['Alignment Score: ', num2str(alignmentScore)]);
disp(['Reverse Complement Sequence: ', reverseComplement]);
disp(['Melting Temperature: ', num2str(meltingTemp), '°C']);
disp(['Motif Enrichment P-Values: ', num2str(enrichmentPValues)]);
disp(['ORF Positions: ', num2str(orfPositions)]);
disp(['RNA Secondary Structure: ', rnaSecondaryStructure]);
disp(['Codon Usage Table:']);
disp(codonUsageTable);
disp(['Protein Molecular Weight: ', num2str(proteinMolecularWeight), ' Da']);
disp(['Motif Search Occurrences: ', num2str(numMotifSearch)]);
disp(['Sequence Similarity Score: ', num2str(similarityScore)]);
disp(['Protein Model: ', proteinModel]);
disp(['Protein Interactions: ', proteinInteractions]);
disp(['Multiple Sequence Alignment:']);
disp(multipleAlignment);
