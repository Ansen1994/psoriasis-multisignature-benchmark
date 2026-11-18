import fs from 'node:fs';
import path from 'node:path';
import zlib from 'node:zlib';
import readline from 'node:readline';
import { fileURLToPath } from 'node:url';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const OUT = path.join(ROOT, 'data', 'benchmark');

const signatures = [
  {
    id: 'RRDS',
    label: 'R-loop-related diagnostic signature',
    pmid: '42898510',
    pmcid: 'PMC13677910',
    sourceRecord: 'sources/raw/literature/PMC13677910.fullText.xml',
    genes: ['SRPK1', 'BLM', 'SAMD9', 'TPX2', 'KIF23', 'APOBEC3B', 'CDCA8', 'MKI67'],
    originalEndpoint: 'D2',
    compatibleEndpoints: 'D2',
    evidenceGrade: 'B',
    modelStatus: 'no_recovered_weights',
    comparator: 'equal_weight_zscore',
    derivationCohorts: 'GSE30999;GSE41664',
    reportedValidationCohorts: 'GSE13355;GSE14905;GSE6710',
    notes: 'The source describes a diagnostic model, but its main comparisons are LS versus NL; D1 wording is inconsistent and is not treated as a separate validated endpoint.',
  },
  {
    id: 'RiboSisRG',
    label: 'Ribosome-biogenesis-related six-gene model',
    pmid: '41217566',
    pmcid: 'PMC12605623',
    sourceRecord: 'sources/raw/literature/PMC12605623.fullText.xml',
    genes: ['MXD1', 'SCO2', 'FOSL1', 'STAT1', 'LTF', 'DEPDC1B'],
    originalEndpoint: 'D2',
    compatibleEndpoints: 'D2',
    evidenceGrade: 'B',
    modelStatus: 'no_recovered_weights',
    comparator: 'equal_weight_zscore',
    derivationCohorts: 'GSE30999;GSE13355',
    reportedValidationCohorts: 'GSE14905',
    notes: 'Six-gene diagnostic model; source reports a high-accuracy diagnostic claim, but the bulk endpoint is LS versus NL.',
  },
  {
    id: 'TRGS',
    label: 'Treg-associated five-gene signature',
    pmid: '41572692',
    pmcid: 'PMC13575954',
    sourceRecord: 'sources/raw/literature/PMC13575954.fullText.xml',
    genes: ['CRIP1', 'FBXW11', 'CD47', 'ECH1', 'H3F3A'],
    originalEndpoint: 'D2',
    compatibleEndpoints: 'D2',
    evidenceGrade: 'B',
    modelStatus: 'no_recovered_weights',
    comparator: 'equal_weight_zscore',
    derivationCohorts: 'GSE13355',
    reportedValidationCohorts: 'GSE14905;GSE30999',
    notes: 'The paper reports a diagnostic model and a PASI correlation; only the diagnostic/LS-NL component is used for the D2 benchmark.',
  },
  {
    id: 'IMscore',
    label: 'Immune-microenvironment five-gene score',
    pmid: '41857773',
    pmcid: 'PMC13107425',
    sourceRecord: 'sources/raw/literature/PMC13107425.fullText.xml',
    genes: ['CXCL8', 'CXCL9', 'CCL18', 'RGS1', 'SAMSN1'],
    originalEndpoint: 'D2',
    compatibleEndpoints: 'D2',
    evidenceGrade: 'B',
    modelStatus: 'no_recovered_weights',
    comparator: 'equal_weight_zscore',
    derivationCohorts: 'GSE106992;GSE117239;GSE117468;GSE13355;GSE14905;GSE226244;GSE30999;GSE34248;GSE41662;GSE53552;GSE78097',
    reportedValidationCohorts: '',
    notes: 'Pooled multi-cohort risk model; the authors acknowledge no genuinely independent large external cohort. GSE41662 is a subseries of GSE41664 and cannot be called an unseen GSE41664 holdout.',
  },
  {
    id: 'Ferroptosis8',
    label: 'Ferroptosis-related eight-gene model',
    pmid: '42948308',
    pmcid: 'PMC13712149',
    sourceRecord: 'sources/raw/literature/PMC13712149.fullText.xml',
    genes: ['CD24', 'EFNB2', 'DDX39A', 'NDRG4', 'GPRC5A', 'S100A8', 'PRSS8', 'ITPKC'],
    originalEndpoint: 'D2/S',
    compatibleEndpoints: 'D2;S',
    evidenceGrade: 'B',
    modelStatus: 'no_recovered_weights',
    comparator: 'equal_weight_zscore',
    derivationCohorts: 'GSE30999;GSE78097',
    reportedValidationCohorts: 'GSE13355;GSE14905;GSE41662;GSE54456',
    notes: 'The paper combines diagnostic and severity analyses and uses GSE41662, a subseries of GSE41664. D2 and S must remain separate.',
  },
  {
    id: 'HKG34',
    label: 'Housekeeping-gene avNNet feature set',
    pmid: '40959079',
    pmcid: 'PMC12433973',
    sourceRecord: 'sources/raw/literature/PMC12433973.fullText.xml',
    genes: [
      'ABCF2', 'ALDH3A2', 'APOL2', 'CD81', 'CHPT1', 'DCUN1D3', 'EIF5A', 'FKBPL', 'GFOD2', 'GLRX3',
      'GLTP', 'GTF3C6', 'HEBP1', 'HIGD1A', 'HMOX2', 'LAP3', 'MED8', 'MICALL1', 'MYO1E', 'NOP2', 'NR3C2',
      'NUDT15', 'PARP9', 'PDGFC', 'PINX1', 'PPARD', 'PPIF', 'RAN', 'RUVBL1', 'SCO2', 'STAT3', 'TXNDC17',
      'UBE2F', 'ZNF134',
    ],
    originalEndpoint: 'D1',
    compatibleEndpoints: 'D1;S',
    evidenceGrade: 'B',
    modelStatus: 'no_recovered_weights',
    comparator: 'equal_weight_zscore',
    derivationCohorts: 'GSE54456',
    reportedValidationCohorts: 'GSE121212;GSE109248;GSE66511;GSE80047;GSE53431',
    notes: 'The paper reports an avNNet model; the 34-gene feature set is recoverable, but the trained network weights and intercept are not. Equal weighting is a comparator, not a model reconstruction.',
  },
  {
    id: 'Epidermal10',
    label: 'Epidermal-barrier-related ten-gene model',
    pmid: '42974488',
    pmcid: 'PMC13726643',
    sourceRecord: 'sources/raw/literature/PMC13726643.fullText.xml',
    genes: ['UPP1', 'CLDN1', 'SMPD3', 'PPARD', 'CCNE1', 'CARHSP1', 'GDPD3', 'PARD3', 'EZH2', 'PGAM1'],
    originalEndpoint: 'D1',
    compatibleEndpoints: 'D1',
    evidenceGrade: 'B',
    modelStatus: 'no_recovered_weights',
    comparator: 'equal_weight_zscore',
    derivationCohorts: 'GSE54456',
    reportedValidationCohorts: 'GSE30999',
    notes: 'The paper reports extremely high AUC; static classification and survival-oriented StepCox/CoxBoost terminology require methods review before using the original claim.',
  },
  {
    id: 'Apoptosis5',
    label: 'Apoptosis-proliferation five-gene model',
    pmid: '42353157',
    pmcid: 'PMC13299817',
    sourceRecord: 'sources/raw/literature/PMC13299817.fullText.xml',
    genes: ['CCNB1', 'KIF11', 'HDAC1', 'TPX2', 'MELK'],
    originalEndpoint: 'D1/D2',
    compatibleEndpoints: 'D1;D2',
    evidenceGrade: 'B',
    modelStatus: 'no_recovered_weights',
    comparator: 'equal_weight_zscore',
    derivationCohorts: 'GSE30999;GSE53552',
    reportedValidationCohorts: 'GSE55201',
    notes: 'The external cohort is whole blood, not skin. Cross-tissue performance cannot be treated as equivalent to skin D2 validation.',
  },
  {
    id: 'm6A10',
    label: 'm6A-regulator ten-gene model',
    pmid: '38436011',
    pmcid: 'PMC10909359',
    sourceRecord: 'sources/raw/literature/PMC10909359.fullText.xml',
    genes: ['FTO', 'IGF2BP2', 'METTL3', 'YTHDC1', 'ZC3H13', 'HNRNPC', 'IGF2BP3', 'LRPPRC', 'YTHDC2', 'HNRNPA2B1'],
    originalEndpoint: 'D2',
    compatibleEndpoints: 'D2',
    evidenceGrade: 'B',
    modelStatus: 'no_recovered_weights',
    comparator: 'equal_weight_zscore',
    derivationCohorts: 'GSE30999',
    reportedValidationCohorts: 'GSE13355',
    notes: 'The source explicitly contrasts LS with NL. The reported validation AUC is 0.730 and must not be replaced by a stronger proxy result.',
  },
  {
    id: 'PLOS14',
    label: 'LASSO fourteen-gene lesion/nonlesional model',
    pmid: '40560938',
    pmcid: 'PMC12192183',
    sourceRecord: 'sources/raw/literature/PMC12192183.fullText.xml',
    genes: ['ABCG4', 'CD24', 'DEFB103A', 'EPHB2', 'FABP5', 'HPSE', 'IFI16', 'KYNU', 'OAS3', 'OASL', 'PLA2G4D', 'SERPINB13', 'STAT1', 'TMPRSS11D'],
    originalEndpoint: 'D2',
    compatibleEndpoints: 'D2;T',
    evidenceGrade: 'B',
    modelStatus: 'no_recovered_weights',
    comparator: 'equal_weight_zscore',
    derivationCohorts: 'GSE136757;NCT02310750',
    reportedValidationCohorts: 'GSE30999;GSE6710;GSE14905',
    notes: 'The fourteen-gene LS/NL model is separate from the five-gene treatment-response feature BUB1, HERC6, ISG15, MXD1, NCAPG. The latter is not counted as a fourteenth-gene component of the diagnostic model.',
  },
];

const cohorts = [
  {
    cohort: 'GSE13355', platform: 'GPL570', totalSamples: 180,
    lesional: 58, nonlesional: 58, normal: 64, completePairs: 58,
    nontreatmentPairs: 58, treatmentBaselinePairs: 0, pretreatmentPairs: 58,
    posttreatmentLesionalSamples: 0, baselineRemissionPairs: 0, remissionRelapsePairs: 0,
    d1: 'yes', d2: 'yes', severity: 'no', treatment: 'no',
    role: 'D1/D2 cohort heavily reused across signatures; not a clean independent holdout for most named models.',
    sourceFile: 'sources/raw/GSE13355/GSE13355_series_matrix.txt.gz',
    verification: 'verified_series_matrix',
  },
  {
    cohort: 'GSE30999', platform: 'GPL570', totalSamples: 170,
    lesional: 85, nonlesional: 85, normal: 0, completePairs: 81,
    nontreatmentPairs: 81, treatmentBaselinePairs: 0, pretreatmentPairs: 81,
    posttreatmentLesionalSamples: 0, baselineRemissionPairs: 0, remissionRelapsePairs: 0,
    d1: 'no', d2: 'yes', severity: 'possible_unverified', treatment: 'no',
    role: 'Important D2 discovery cohort; repeated reuse means performance here cannot be called independent across signatures.',
    sourceFile: 'sources/raw/GSE30999/GSE30999_series_matrix.txt.gz',
    verification: 'verified_series_matrix',
  },
  {
    cohort: 'GSE14905', platform: 'GPL570', totalSamples: 82,
    lesional: 33, nonlesional: 28, normal: 21, completePairs: 28,
    nontreatmentPairs: 28, treatmentBaselinePairs: 0, pretreatmentPairs: 28,
    posttreatmentLesionalSamples: 0, baselineRemissionPairs: 0, remissionRelapsePairs: 0,
    d1: 'yes', d2: 'yes', severity: 'no', treatment: 'no',
    role: 'D1/D2 validation cohort; used by several signature papers and by the project as a prior holdout.',
    sourceFile: 'sources/raw/GSE14905/GSE14905_series_matrix.txt.gz',
    verification: 'verified_series_matrix',
  },
  {
    cohort: 'GSE41664', platform: 'GPL570', totalSamples: 157,
    lesional: 104, nonlesional: 53, normal: 0, completePairs: 53,
    nontreatmentPairs: 38, treatmentBaselinePairs: 15, pretreatmentPairs: 53,
    posttreatmentLesionalSamples: 51, baselineRemissionPairs: 0, remissionRelapsePairs: 0,
    d1: 'no', d2: 'yes', severity: 'no', treatment: 'yes',
    role: 'SuperSeries with 53 pre-treatment LS/NL pairs (38 from the non-treatment series and 15 at GSE11903 trial baseline) plus 51 post-treatment lesional samples. Clean unseen D2 holdout for only 7/10 candidates.',
    sourceFile: 'sources/raw/GSE41664/GSE41664_series_matrix.txt.gz',
    verification: 'verified_series_matrix',
  },
  {
    cohort: 'GSE6710', platform: 'GPL96', totalSamples: 26,
    lesional: 13, nonlesional: 13, normal: 0, completePairs: 13,
    nontreatmentPairs: 13, treatmentBaselinePairs: 0, pretreatmentPairs: 13,
    posttreatmentLesionalSamples: 0, baselineRemissionPairs: 0, remissionRelapsePairs: 0,
    d1: 'no', d2: 'yes', severity: 'no', treatment: 'no',
    role: 'Small cross-platform D2 holdout; clean unseen cohort for 8/10 candidates, but underpowered and GPL96 coverage must be checked.',
    sourceFile: 'sources/raw/GSE6710/GSE6710_series_matrix.txt.gz',
    verification: 'verified_series_matrix',
  },
  {
    cohort: 'GSE78097', platform: 'GPL570', totalSamples: 33,
    lesional: 27, nonlesional: 0, normal: 6, completePairs: 0,
    nontreatmentPairs: 0, treatmentBaselinePairs: 0, pretreatmentPairs: 0,
    posttreatmentLesionalSamples: 0, baselineRemissionPairs: 0, remissionRelapsePairs: 0,
    d1: 'not_primary', d2: 'no', severity: 'yes', treatment: 'no',
    role: 'Mild (14) versus severe (13) lesional cohort plus 6 normal samples; use for severity only, not diagnostic AUC.',
    sourceFile: 'sources/raw/GSE78097/GSE78097_series_matrix.txt.gz',
    verification: 'verified_series_matrix',
  },
  {
    cohort: 'GSE212476', platform: 'RNA-seq FPKM', totalSamples: 40,
    lesional: 40, nonlesional: 0, normal: 0, completePairs: 0,
    nontreatmentPairs: 0, treatmentBaselinePairs: 0, pretreatmentPairs: 0,
    posttreatmentLesionalSamples: 0, baselineRemissionPairs: 14, remissionRelapsePairs: 12,
    d1: 'no', d2: 'no', severity: 'no', treatment: 'yes',
    role: 'Baseline/remission/relapse longitudinal samples with 14 baseline-remission and 12 remission-relapse pairs; use for reversibility and treatment boundary, not D2 diagnosis.',
    sourceFile: 'sources/raw/GSE212476/GSE212476_submission_fpkm.csv.gz',
    verification: 'verified_local_processed_matrix',
  },
];

function csvEscape(value) {
  const text = value == null ? '' : String(value);
  return /[",\r\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

function writeCsv(file, rows, columns) {
  const lines = [columns.join(','), ...rows.map((row) => columns.map((column) => csvEscape(row[column])).join(','))];
  fs.mkdirSync(path.dirname(file), { recursive: true });
  fs.writeFileSync(file, `${lines.join('\n')}\n`, 'utf8');
}

function openLines(file) {
  const source = fs.createReadStream(file);
  const input = file.endsWith('.gz') ? source.pipe(zlib.createGunzip()) : source;
  const reader = readline.createInterface({ input, crlfDelay: Infinity });
  return { reader, input };
}

async function parseAnnotation(file, expectedPlatform) {
  const { reader, input } = openLines(file);
  const bySymbol = new Map();
  let inTable = false;
  let symbolIndex = -1;
  let headerSeen = false;
  for await (const rawLine of reader) {
    const line = rawLine.replace(/^\uFEFF/, '');
    if (!inTable) {
      if (line === '!platform_table_begin') inTable = true;
      continue;
    }
    if (!headerSeen) {
      const headers = line.split('\t');
      const index = headers.findIndex((header) => /^gene[ _]symbol$/i.test(header.trim()));
      if (index < 0) throw new Error(`No Gene symbol column in ${file}`);
      symbolIndex = index;
      headerSeen = true;
      continue;
    }
    if (line === '!platform_table_end') break;
    const fields = line.split('\t');
    const probe = fields[0];
    const symbols = (fields[symbolIndex] || '')
      .split('///')
      .map((symbol) => symbol.trim())
      .filter((symbol) => symbol && symbol !== '---');
    for (const symbol of new Set(symbols)) {
      if (!bySymbol.has(symbol)) bySymbol.set(symbol, new Set());
      bySymbol.get(symbol).add(probe);
    }
  }
  input.destroy();
  return { platform: expectedPlatform, bySymbol };
}

function sourceMentions(relativePath) {
  const absolute = path.join(ROOT, relativePath);
  if (!fs.existsSync(absolute)) return [];
  return [...new Set(fs.readFileSync(absolute, 'utf8').match(/GSE\d+/g) || [])].sort();
}

function cohortUsedByOriginal(signature, cohort) {
  const declared = `${signature.derivationCohorts};${signature.reportedValidationCohorts}`.split(';').filter(Boolean);
  if (declared.includes(cohort)) return 'yes';
  if (cohort === 'GSE41664' && (declared.includes('GSE34248') || declared.includes('GSE41662'))) return 'subseries_used';
  return 'no';
}

function gateAuditForCohort(cohortId, coverageRows) {
  const noDeclaredUse = signatures.filter((signature) => cohortUsedByOriginal(signature, cohortId) === 'no');
  const d2CompatibleNoDeclaredUse = noDeclaredUse.filter((signature) => signature.compatibleEndpoints.split(';').includes('D2'));
  const completeGeneD2Comparable = d2CompatibleNoDeclaredUse.filter((signature) => {
    const coverage = coverageRows.find((row) => row.signature_id === signature.id);
    return coverage && coverage.missing_gene_count === 0;
  });
  return {
    noDeclaredUseCount: noDeclaredUse.length,
    noDeclaredUseSignatures: noDeclaredUse.map((signature) => signature.id),
    d2CompatibleNoDeclaredUseCount: d2CompatibleNoDeclaredUse.length,
    d2CompatibleNoDeclaredUseSignatures: d2CompatibleNoDeclaredUse.map((signature) => signature.id),
    completeGeneD2ComparableCount: completeGeneD2Comparable.length,
    completeGeneD2ComparableSignatures: completeGeneD2Comparable.map((signature) => signature.id),
  };
}

function rowsForUsage() {
  const rows = [];
  for (const signature of signatures) {
    const mentions = sourceMentions(signature.sourceRecord);
    for (const cohort of cohorts) {
      rows.push({
        signature_id: signature.id,
        cohort: cohort.cohort,
        declared_in_registry: cohortUsedByOriginal(signature, cohort.cohort),
        mentioned_in_source: mentions.includes(cohort.cohort) ? 'yes' : 'no',
        source_mentions: mentions.join(';'),
      });
    }
  }
  return rows;
}

function rowsForCoverage(annotation) {
  const rows = [];
  for (const signature of signatures) {
    const present = [];
    const missing = [];
    const multi = [];
    const probes = new Set();
    for (const gene of signature.genes) {
      const geneProbes = annotation.bySymbol.get(gene);
      if (!geneProbes || geneProbes.size === 0) {
        missing.push(gene);
      } else {
        present.push(gene);
        if (geneProbes.size > 1) multi.push(gene);
        for (const probe of geneProbes) probes.add(probe);
      }
    }
    rows.push({
      platform: annotation.platform,
      signature_id: signature.id,
      signature_label: signature.label,
      signature_size: signature.genes.length,
      present_gene_count: present.length,
      missing_gene_count: missing.length,
      coverage_percent: ((present.length / signature.genes.length) * 100).toFixed(2),
      multi_probe_gene_count: multi.length,
      total_unique_probes: probes.size,
      present_genes: present.join(';'),
      missing_genes: missing.join(';'),
      multi_probe_genes: multi.join(';'),
    });
  }
  return rows;
}

function rowsForOverlap() {
  const rows = [];
  for (let i = 0; i < signatures.length; i += 1) {
    for (let j = i + 1; j < signatures.length; j += 1) {
      const a = new Set(signatures[i].genes);
      const b = new Set(signatures[j].genes);
      const shared = [...a].filter((gene) => b.has(gene)).sort();
      const union = new Set([...a, ...b]);
      rows.push({
        signature_a: signatures[i].id,
        signature_b: signatures[j].id,
        signature_a_size: a.size,
        signature_b_size: b.size,
        shared_gene_count: shared.length,
        jaccard: union.size ? (shared.length / union.size).toFixed(4) : '0.0000',
        shared_genes: shared.join(';'),
      });
    }
  }
  return rows;
}

function matrixRows() {
  return signatures.map((signature) => {
    const row = { signature_id: signature.id };
    for (const other of signatures) {
      row[other.id] = other.id === signature.id
        ? signature.genes.length
        : [...new Set(signature.genes)].filter((gene) => other.genes.includes(gene)).length;
    }
    return row;
  });
}

async function main() {
  const gpl570 = await parseAnnotation(path.join(ROOT, 'sources', 'raw', 'GPL570', 'GPL570.annot.gz'), 'GPL570');
  const gpl96 = await parseAnnotation(path.join(ROOT, 'sources', 'raw', 'GPL96', 'GPL96_full.txt'), 'GPL96');

  const registryRows = signatures.map((signature) => ({
    signature_id: signature.id,
    label: signature.label,
    pmid: signature.pmid,
    pmcid: signature.pmcid,
    source_record: signature.sourceRecord,
    evidence_grade: signature.evidenceGrade,
    model_status: signature.modelStatus,
    comparator: signature.comparator,
    original_endpoint: signature.originalEndpoint,
    compatible_endpoints: signature.compatibleEndpoints,
    gene_count: signature.genes.length,
    genes: signature.genes.join(';'),
    derivation_cohorts: signature.derivationCohorts,
    reported_validation_cohorts: signature.reportedValidationCohorts,
    source_mentions: sourceMentions(signature.sourceRecord).join(';'),
    notes: signature.notes,
  }));

  const cohortRows = cohorts.map((cohort) => ({
    cohort: cohort.cohort,
    platform: cohort.platform,
    total_samples: cohort.totalSamples,
    lesional_samples: cohort.lesional,
    nonlesional_samples: cohort.nonlesional,
    normal_samples: cohort.normal,
    complete_lesional_nonlesional_pairs: cohort.completePairs,
    nontreatment_lesional_nonlesional_pairs: cohort.nontreatmentPairs,
    treatment_baseline_lesional_nonlesional_pairs: cohort.treatmentBaselinePairs,
    pretreatment_lesional_nonlesional_pairs: cohort.pretreatmentPairs,
    posttreatment_lesional_samples: cohort.posttreatmentLesionalSamples,
    longitudinal_baseline_remission_pairs: cohort.baselineRemissionPairs,
    longitudinal_remission_relapse_pairs: cohort.remissionRelapsePairs,
    d1_lesion_vs_healthy: cohort.d1,
    d2_lesion_vs_nonlesional: cohort.d2,
    s_severity: cohort.severity,
    t_treatment_or_reversibility: cohort.treatment,
    benchmark_role: cohort.role,
    source_file: cohort.sourceFile,
    verification: cohort.verification,
  }));

  writeCsv(path.join(OUT, 'signature_registry.csv'), registryRows, [
    'signature_id', 'label', 'pmid', 'pmcid', 'source_record', 'evidence_grade', 'model_status', 'comparator',
    'original_endpoint', 'compatible_endpoints', 'gene_count', 'genes', 'derivation_cohorts',
    'reported_validation_cohorts', 'source_mentions', 'notes',
  ]);
  writeCsv(path.join(OUT, 'cohort_endpoint_matrix.csv'), cohortRows, [
    'cohort', 'platform', 'total_samples', 'lesional_samples', 'nonlesional_samples', 'normal_samples',
    'complete_lesional_nonlesional_pairs', 'nontreatment_lesional_nonlesional_pairs',
    'treatment_baseline_lesional_nonlesional_pairs', 'pretreatment_lesional_nonlesional_pairs',
    'posttreatment_lesional_samples', 'longitudinal_baseline_remission_pairs',
    'longitudinal_remission_relapse_pairs', 'd1_lesion_vs_healthy', 'd2_lesion_vs_nonlesional',
    's_severity', 't_treatment_or_reversibility', 'benchmark_role', 'source_file', 'verification',
  ]);
  writeCsv(path.join(OUT, 'signature_cohort_usage.csv'), rowsForUsage(), [
    'signature_id', 'cohort', 'declared_in_registry', 'mentioned_in_source', 'source_mentions',
  ]);
  writeCsv(path.join(OUT, 'platform_gene_coverage.csv'), [...rowsForCoverage(gpl570), ...rowsForCoverage(gpl96)], [
    'platform', 'signature_id', 'signature_label', 'signature_size', 'present_gene_count', 'missing_gene_count',
    'coverage_percent', 'multi_probe_gene_count', 'total_unique_probes', 'present_genes', 'missing_genes',
    'multi_probe_genes',
  ]);
  const overlapRows = rowsForOverlap();
  writeCsv(path.join(OUT, 'signature_overlap_long.csv'), overlapRows, [
    'signature_a', 'signature_b', 'signature_a_size', 'signature_b_size', 'shared_gene_count', 'jaccard', 'shared_genes',
  ]);
  const matrix = matrixRows();
  writeCsv(path.join(OUT, 'signature_overlap_matrix.csv'), matrix, ['signature_id', ...signatures.map((signature) => signature.id)]);

  const summary = {
    generatedAt: new Date().toISOString(),
    decision: 'CONDITIONAL_GO',
    signatures: signatures.length,
    exactModelWeightsRecovered: signatures.filter((signature) => signature.modelStatus !== 'no_recovered_weights').length,
    evidenceGrades: {
      A: signatures.filter((signature) => signature.evidenceGrade === 'A').length,
      B: signatures.filter((signature) => signature.evidenceGrade === 'B').length,
      C: signatures.filter((signature) => signature.evidenceGrade === 'C').length,
    },
    leakageAndEndpointAudit: {
      GSE41664: gateAuditForCohort('GSE41664', rowsForCoverage(gpl570)),
      GSE6710: gateAuditForCohort('GSE6710', rowsForCoverage(gpl96)),
    },
    geneOverlapAudit: {
      pairwiseComparisons: overlapRows.length,
      nonzeroSharedGenePairs: overlapRows.filter((row) => row.shared_gene_count > 0).length,
      maximumSharedGenes: Math.max(...overlapRows.map((row) => row.shared_gene_count)),
    },
    claimBoundary: {
      allowed: 'Fixed-gene-list transparency benchmark using a prespecified equal-weight z-score comparator.',
      forbidden: [
        'Do not describe equal-weight scores as exact reconstruction of published models.',
        'Do not pool D1, D2, S, and T tasks.',
        'Do not describe all signatures as independently validated in GSE41664 or GSE6710.',
      ],
    },
    cohorts,
    coverage: {
      GPL570: rowsForCoverage(gpl570),
      GPL96: rowsForCoverage(gpl96),
    },
  };
  fs.writeFileSync(path.join(OUT, 'benchmark_gate_summary.json'), `${JSON.stringify(summary, null, 2)}\n`, 'utf8');
  console.log(JSON.stringify({ outputs: OUT, signatures: signatures.length, gpl570Genes: gpl570.bySymbol.size, gpl96Genes: gpl96.bySymbol.size }, null, 2));
}

main().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});
