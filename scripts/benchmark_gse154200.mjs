import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const STAGE = path.join(ROOT, 'data', 'benchmark', 'score_stage_independent', 'GSE154200');
const OUT = path.join(STAGE, 'increment_stage');
const LAMBDA = 1.0;
const BOOTSTRAP_REPLICATES = 2000;
const MIN_MEANINGFUL_DELTA_AUC = 0.03;

const PRIMARY_SIGNATURES = [
  'RRDS', 'RiboSisRG', 'TRGS', 'IMscore',
  'Ferroptosis8', 'Apoptosis5', 'm6A10', 'PLOS14',
];

const GENE_SETS = [
  { id: 'KC', file: 'sources/raw/msigdb/c5.go.v2025.1.Hs.symbols.gmt', term: 'GOBP_KERATINOCYTE_DIFFERENTIATION' },
  { id: 'Inflammatory', file: 'sources/raw/msigdb/h.all.v2025.1.Hs.symbols.gmt', term: 'HALLMARK_INFLAMMATORY_RESPONSE' },
  { id: 'IFN_gamma', file: 'sources/raw/msigdb/h.all.v2025.1.Hs.symbols.gmt', term: 'HALLMARK_INTERFERON_GAMMA_RESPONSE' },
  { id: 'Proliferation', file: 'sources/raw/msigdb/h.all.v2025.1.Hs.symbols.gmt', term: 'HALLMARK_E2F_TARGETS' },
];

function parseCsv(text) {
  const rows = [];
  let row = [];
  let field = '';
  let quoted = false;
  for (let index = 0; index < text.length; index += 1) {
    const char = text[index];
    if (quoted) {
      if (char === '"' && text[index + 1] === '"') {
        field += '"';
        index += 1;
      } else if (char === '"') {
        quoted = false;
      } else {
        field += char;
      }
    } else if (char === '"') {
      quoted = true;
    } else if (char === ',') {
      row.push(field);
      field = '';
    } else if (char === '\n') {
      row.push(field.replace(/\r$/, ''));
      rows.push(row);
      row = [];
      field = '';
    } else {
      field += char;
    }
  }
  if (field.length > 0 || row.length > 0) {
    row.push(field.replace(/\r$/, ''));
    rows.push(row);
  }
  const headers = rows.shift() || [];
  return rows.filter((values) => values.some((value) => value.length > 0)).map((values) => {
    const record = {};
    headers.forEach((header, index) => { record[header] = values[index] ?? ''; });
    return record;
  });
}

function readCsv(relativePath) {
  return parseCsv(fs.readFileSync(path.join(ROOT, relativePath), 'utf8'));
}

function writeCsv(file, rows, columns) {
  const escape = (value) => {
    const text = value == null ? '' : String(value);
    return /[",\r\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
  };
  fs.mkdirSync(path.dirname(file), { recursive: true });
  const lines = [columns.join(','), ...rows.map((row) => columns.map((column) => escape(row[column])).join(','))];
  fs.writeFileSync(file, `${lines.join('\n')}\n`, 'utf8');
}

function finiteMean(values) {
  let sum = 0;
  let count = 0;
  for (const value of values) {
    if (!Number.isFinite(value)) continue;
    sum += value;
    count += 1;
  }
  return count > 0 ? sum / count : null;
}

function finiteSd(values) {
  const mean = finiteMean(values);
  if (mean == null) return null;
  let sum = 0;
  let count = 0;
  for (const value of values) {
    if (!Number.isFinite(value)) continue;
    sum += (value - mean) ** 2;
    count += 1;
  }
  return count > 1 ? Math.sqrt(sum / (count - 1)) : null;
}

function percentile(values, probability) {
  const sorted = values.filter(Number.isFinite).sort((a, b) => a - b);
  if (sorted.length === 0) return null;
  const position = (sorted.length - 1) * probability;
  const lower = Math.floor(position);
  const upper = Math.ceil(position);
  if (lower === upper) return sorted[lower];
  return sorted[lower] + (sorted[upper] - sorted[lower]) * (position - lower);
}

function median(values) {
  return percentile(values, 0.5);
}

function formatNumber(value, digits = 6) {
  return Number.isFinite(value) ? Number(value).toFixed(digits) : '';
}

function readGmtTerm(relativePath, term) {
  const line = fs.readFileSync(path.join(ROOT, relativePath), 'utf8')
    .split(/\r?\n/)
    .find((candidate) => candidate.startsWith(`${term}\t`));
  if (!line) throw new Error(`GMT term not found: ${term}`);
  return [...new Set(line.split('\t').slice(2).map((gene) => gene.trim()).filter(Boolean))];
}

function fitGeneStats(genes, geneValues, sampleIndices) {
  const stats = new Map();
  for (const gene of genes) {
    const values = geneValues.get(gene);
    if (!values) continue;
    const reference = sampleIndices.map((index) => values[index]).filter(Number.isFinite);
    const mean = finiteMean(reference);
    const sd = finiteSd(reference);
    if (mean != null && sd != null && sd > 0) stats.set(gene, { mean, sd });
  }
  return stats;
}

function scoreGeneSet(genes, geneValues, stats, sampleIndices) {
  return sampleIndices.map((sampleIndex) => {
    const z = genes.map((gene) => {
      const value = geneValues.get(gene)?.[sampleIndex];
      const stat = stats.get(gene);
      return stat && Number.isFinite(value) ? (value - stat.mean) / stat.sd : null;
    }).filter(Number.isFinite);
    return z.length ? z.reduce((sum, value) => sum + value, 0) / z.length : null;
  });
}

function aucFromRows(rows, predictionKey) {
  const usable = rows.filter((row) => row.label != null && Number.isFinite(Number(row[predictionKey])));
  const positives = usable.filter((row) => row.label === 1);
  const negatives = usable.filter((row) => row.label === 0);
  if (!positives.length || !negatives.length) return null;
  const ordered = usable.map((row) => ({ score: Number(row[predictionKey]), label: row.label })).sort((a, b) => a.score - b.score);
  let index = 0;
  while (index < ordered.length) {
    let end = index + 1;
    while (end < ordered.length && ordered[end].score === ordered[index].score) end += 1;
    const rank = ((index + 1) + end) / 2;
    for (let rankIndex = index; rankIndex < end; rankIndex += 1) ordered[rankIndex].rank = rank;
    index = end;
  }
  const rankSum = ordered.filter((row) => row.label === 1).reduce((sum, row) => sum + row.rank, 0);
  return (rankSum - positives.length * (positives.length + 1) / 2) / (positives.length * negatives.length);
}

function stableSigmoid(value) {
  if (value >= 0) return 1 / (1 + Math.exp(-value));
  const exp = Math.exp(value);
  return exp / (1 + exp);
}

function fitNonnegativeRidgeLogistic(xRows, labels, lambda) {
  const featureCount = xRows[0].length;
  const positiveFraction = (labels.reduce((sum, label) => sum + label, 0) + 0.5) / (labels.length + 1);
  const coefficients = [Math.log(positiveFraction / (1 - positiveFraction)), ...Array(featureCount).fill(0)];
  let converged = false;
  let iterations = 0;
  for (iterations = 1; iterations <= 5000; iterations += 1) {
    let maxChange = 0;
    for (let coefficientIndex = 0; coefficientIndex <= featureCount; coefficientIndex += 1) {
      let gradient = 0;
      let curvature = 0;
      for (let rowIndex = 0; rowIndex < xRows.length; rowIndex += 1) {
        const xRow = xRows[rowIndex];
        const x = coefficientIndex === 0 ? 1 : xRow[coefficientIndex - 1];
        const linear = coefficients[0] + xRow.reduce((sum, value, index) => sum + coefficients[index + 1] * value, 0);
        const probability = stableSigmoid(linear);
        gradient += x * (labels[rowIndex] - probability);
        curvature += x * x * probability * (1 - probability);
      }
      const penalty = coefficientIndex === 0 ? 0 : lambda;
      gradient -= penalty * coefficients[coefficientIndex];
      curvature += penalty;
      if (!(curvature > 1e-12)) continue;
      const newtonStep = Math.max(-10, Math.min(10, gradient / curvature));
      const previous = coefficients[coefficientIndex];
      const proposed = coefficientIndex === 0 ? previous + newtonStep : Math.max(0, previous + newtonStep);
      coefficients[coefficientIndex] = proposed;
      maxChange = Math.max(maxChange, Math.abs(proposed - previous));
    }
    if (maxChange < 1e-8) {
      converged = true;
      break;
    }
  }
  return { coefficients, converged, iterations };
}

function fitStandardization(rows, features) {
  return features.map((feature) => {
    const values = rows.map((row) => Number(row[feature])).filter(Number.isFinite);
    const mean = finiteMean(values);
    const sd = finiteSd(values);
    return { feature, mean, sd: sd != null && sd > 0 ? sd : 1, constant: !(sd > 0) };
  });
}

function standardizedVector(row, stats) {
  return stats.map((stat) => stat.constant ? 0 : (Number(row[stat.feature]) - stat.mean) / stat.sd);
}

function crossValidateByPatient(rows, features) {
  const patientIds = [...new Set(rows.map((row) => row.patient_id))].sort();
  const predictions = [];
  const foldCoefficients = [];
  for (const patientId of patientIds) {
    const trainRows = rows.filter((row) => row.patient_id !== patientId);
    const testRows = rows.filter((row) => row.patient_id === patientId);
    const stats = fitStandardization(trainRows, features);
    const xTrain = trainRows.map((row) => standardizedVector(row, stats));
    const labels = trainRows.map((row) => row.label);
    const fit = fitNonnegativeRidgeLogistic(xTrain, labels, LAMBDA);
    foldCoefficients.push({
      held_out_patient: patientId,
      intercept: fit.coefficients[0],
      coefficients: Object.fromEntries(features.map((feature, index) => [feature, fit.coefficients[index + 1]])),
      converged: fit.converged,
      iterations: fit.iterations,
    });
    for (const row of testRows) {
      const x = standardizedVector(row, stats);
      const linear = fit.coefficients[0] + x.reduce((sum, value, index) => sum + value * fit.coefficients[index + 1], 0);
      predictions.push({ ...row, prediction: stableSigmoid(linear) });
    }
  }
  return { predictions, foldCoefficients };
}

function createSeededRandom(seed) {
  let state = seed >>> 0;
  return () => {
    state ^= state << 13;
    state ^= state >>> 17;
    state ^= state << 5;
    return (state >>> 0) / 4294967296;
  };
}

function patientBootstrap(predictionRows, modelIds, replicates, seed) {
  const byPatient = new Map();
  for (const row of predictionRows) {
    if (!byPatient.has(row.patient_id)) byPatient.set(row.patient_id, []);
    byPatient.get(row.patient_id).push(row);
  }
  const patientIds = [...byPatient.keys()].sort();
  const random = createSeededRandom(seed);
  const aucDraws = Object.fromEntries(modelIds.map((modelId) => [modelId, []]));
  const deltaDraws = [];
  for (let replicate = 0; replicate < replicates; replicate += 1) {
    const sampledRows = [];
    for (let index = 0; index < patientIds.length; index += 1) {
      const patientId = patientIds[Math.floor(random() * patientIds.length)];
      sampledRows.push(...byPatient.get(patientId));
    }
    const aucs = {};
    for (const modelId of modelIds) {
      aucs[modelId] = aucFromRows(sampledRows, modelId);
      if (Number.isFinite(aucs[modelId])) aucDraws[modelId].push(aucs[modelId]);
    }
    if (Number.isFinite(aucs.baseline_combined) && Number.isFinite(aucs.signature_plus_baselines)) {
      deltaDraws.push(aucs.signature_plus_baselines - aucs.baseline_combined);
    }
  }
  return { aucDraws, deltaDraws };
}

function bootstrapInterval(values) {
  return { low: percentile(values, 0.025), high: percentile(values, 0.975), replicates: values.length };
}

function coefficientSummary(foldCoefficients, features) {
  return features.map((feature) => {
    const values = foldCoefficients.map((fold) => fold.coefficients[feature]).filter(Number.isFinite);
    return {
      feature,
      folds: values.length,
      positive_folds: values.filter((value) => value > 1e-8).length,
      zero_folds: values.filter((value) => Math.abs(value) <= 1e-8).length,
      median_coefficient: median(values),
      minimum_coefficient: values.length ? Math.min(...values) : null,
      maximum_coefficient: values.length ? Math.max(...values) : null,
    };
  });
}

function validateD2Rows(rows) {
  const byPatient = new Map();
  for (const row of rows) {
    if (!byPatient.has(row.patient_id)) byPatient.set(row.patient_id, []);
    byPatient.get(row.patient_id).push(row);
  }
  if (byPatient.size !== 88) throw new Error(`Expected 88 D2 patients, found ${byPatient.size}.`);
  for (const [patientId, patientRows] of byPatient) {
    const labels = new Set(patientRows.map((row) => row.label));
    if (patientRows.length !== 2 || !labels.has(0) || !labels.has(1)) {
      throw new Error(`Patient ${patientId} does not have exactly one lesional and one nonlesional sample.`);
    }
  }
  return byPatient.size;
}

function main() {
  const startedAt = Date.now();
  const metadataRows = readCsv('data/benchmark/score_stage_independent/GSE154200/D2_sample_metadata.csv');
  const expressionRows = readCsv('data/benchmark/score_stage_independent/GSE154200/target_gene_logcpm.csv');
  const registry = readCsv('data/benchmark/signature_registry.csv');
  if (metadataRows.length !== 176) throw new Error(`Expected 176 metadata rows, found ${metadataRows.length}.`);
  if (new Set(metadataRows.map((row) => row.gsm)).size !== 176) throw new Error('GSM identifiers are not unique.');

  const expressionHeader = Object.keys(expressionRows[0] || {}).filter((key) => key !== 'gene');
  if (expressionHeader.length !== 176 || new Set(expressionHeader).size !== 176) throw new Error('Expression matrix is not 176 unique samples.');
  const expressionIndexByTitle = new Map(expressionHeader.map((title, index) => [title, index]));
  const geneValues = new Map();
  for (const row of expressionRows) {
    const values = expressionHeader.map((title) => Number(row[title]));
    if (values.some((value) => !Number.isFinite(value))) throw new Error(`Non-finite expression value for ${row.gene}.`);
    geneValues.set(row.gene, values);
  }

  const rows = metadataRows.map((row) => {
    const expressionIndex = expressionIndexByTitle.get(row.title);
    if (expressionIndex == null) throw new Error(`Expression column missing for ${row.title}.`);
    if (!['lesional', 'nonlesional'].includes(row.tissue)) throw new Error(`Unexpected tissue: ${row.tissue}.`);
    return { ...row, label: row.tissue === 'lesional' ? 1 : 0, expression_index: expressionIndex };
  });
  const patients = validateD2Rows(rows);
  const allIndices = rows.map((row) => row.expression_index);

  const signatureDefinitions = PRIMARY_SIGNATURES.map((signatureId) => {
    const record = registry.find((row) => row.signature_id === signatureId);
    if (!record) throw new Error(`Signature not found in registry: ${signatureId}`);
    return {
      id: signatureId,
      genes: [...new Set(record.genes.split(';').map((gene) => gene.trim()).filter(Boolean))],
      declared: Number(record.gene_count),
    };
  });
  if (signatureDefinitions.some((definition) => definition.genes.length !== definition.declared)) {
    throw new Error('Registry gene_count does not match the listed genes.');
  }

  const coverageRows = [];
  for (const definition of signatureDefinitions) {
    const used = definition.genes.filter((gene) => geneValues.has(gene));
    const stats = fitGeneStats(used, geneValues, allIndices);
    const scored = scoreGeneSet(used, geneValues, stats, allIndices);
    if (stats.size < 1 || scored.some((value) => !Number.isFinite(value))) throw new Error(`Incomplete score for ${definition.id}.`);
    rows.forEach((row, index) => { row[definition.id] = scored[index]; });
    coverageRows.push({
      feature_type: 'signature', feature_id: definition.id, genes_declared: definition.genes.length,
      genes_used: stats.size, coverage_percent: stats.size / definition.genes.length * 100,
      complete: stats.size === definition.genes.length ? 'yes' : 'no',
      missing_genes: definition.genes.filter((gene) => !stats.has(gene)).join(';'),
    });
  }

  const baselineDefinitions = GENE_SETS.map((definition) => ({ ...definition, genes: readGmtTerm(definition.file, definition.term) }));
  for (const definition of baselineDefinitions) {
    const used = definition.genes.filter((gene) => geneValues.has(gene));
    const stats = fitGeneStats(used, geneValues, allIndices);
    const scored = scoreGeneSet(used, geneValues, stats, allIndices);
    if (used.length < 5 || scored.some((value) => !Number.isFinite(value))) throw new Error(`Incomplete baseline score for ${definition.id}.`);
    rows.forEach((row, index) => { row[definition.id] = scored[index]; });
    coverageRows.push({
      feature_type: 'baseline', feature_id: definition.id, genes_declared: definition.genes.length,
      genes_used: used.length, coverage_percent: used.length / definition.genes.length * 100,
      complete: used.length === definition.genes.length ? 'yes' : 'no',
      missing_genes: definition.genes.filter((gene) => !geneValues.has(gene)).join(';'),
    });
  }

  writeCsv(path.join(OUT, 'D2_sample_scores.csv'), rows, [
    'gsm', 'title', 'patient_id', 'tissue', 'label', 'pair_id', 'batch_id',
    ...PRIMARY_SIGNATURES, ...GENE_SETS.map((definition) => definition.id),
  ]);
  writeCsv(path.join(OUT, 'D2_score_gene_coverage.csv'), coverageRows, [
    'feature_type', 'feature_id', 'genes_declared', 'genes_used', 'coverage_percent', 'complete', 'missing_genes',
  ]);

  const coverageById = Object.fromEntries(coverageRows.map((row) => [row.feature_id, row]));
  const signatureMetrics = [];
  for (const [signatureIndex, signatureId] of PRIMARY_SIGNATURES.entries()) {
    const scoreRows = rows.map((row) => ({
      patient_id: row.patient_id, gsm: row.gsm, label: row.label, score: Number(row[signatureId]),
    }));
    const directAuc = aucFromRows(scoreRows, 'score');
    const bootstrap = patientBootstrap(scoreRows, ['score'], BOOTSTRAP_REPLICATES, 20261118 + signatureIndex);
    const interval = bootstrapInterval(bootstrap.aucDraws.score);
    const coverage = coverageById[signatureId];
    signatureMetrics.push({
      cohort: 'GSE154200', signature_id: signatureId, patients, samples: rows.length,
      genes_declared: coverage.genes_declared, genes_used: coverage.genes_used,
      coverage_percent: formatNumber(coverage.coverage_percent), coverage_status: coverage.complete === 'yes' ? 'complete' : 'partial',
      auc: formatNumber(directAuc), auc_ci_low: formatNumber(interval.low), auc_ci_high: formatNumber(interval.high),
      direction_result: directAuc < 0.5 ? 'inverse_without_reversal' : 'positive_or_neutral',
    });
  }
  writeCsv(path.join(OUT, 'D2_signature_metrics.csv'), signatureMetrics, [
    'cohort', 'signature_id', 'patients', 'samples', 'genes_declared', 'genes_used', 'coverage_percent',
    'coverage_status', 'auc', 'auc_ci_low', 'auc_ci_high', 'direction_result',
  ]);

  const baselineFeatures = GENE_SETS.map((definition) => definition.id);
  const baselineCombined = crossValidateByPatient(rows, baselineFeatures);
  const baselinePredictionByGsm = new Map(baselineCombined.predictions.map((row) => [row.gsm, row.prediction]));
  const modelIds = ['signature_only', ...baselineFeatures.map((feature) => `${feature}_alone`), 'baseline_combined', 'signature_plus_baselines'];
  const metricRows = [];
  const coefficientRows = [];
  const predictionRows = [];
  const incrementBySignature = {};
  let convergenceFailures = baselineCombined.foldCoefficients.filter((fold) => !fold.converged).length;

  for (const [signatureIndex, signatureId] of PRIMARY_SIGNATURES.entries()) {
    const signaturePlusBaselines = crossValidateByPatient(rows, [signatureId, ...baselineFeatures]);
    convergenceFailures += signaturePlusBaselines.foldCoefficients.filter((fold) => !fold.converged).length;
    const combinedPredictionByGsm = new Map(signaturePlusBaselines.predictions.map((row) => [row.gsm, row.prediction]));
    const modelRows = rows.map((row) => {
      const result = {
        cohort: 'GSE154200', signature_id: signatureId, patient_id: row.patient_id, gsm: row.gsm, label: row.label,
        signature_only: Number(row[signatureId]),
        baseline_combined: baselinePredictionByGsm.get(row.gsm),
        signature_plus_baselines: combinedPredictionByGsm.get(row.gsm),
      };
      for (const feature of baselineFeatures) result[`${feature}_alone`] = Number(row[feature]);
      return result;
    });
    if (modelRows.some((row) => modelIds.some((modelId) => !Number.isFinite(row[modelId])))) throw new Error(`Incomplete predictions for ${signatureId}.`);

    const pointAucs = Object.fromEntries(modelIds.map((modelId) => [modelId, aucFromRows(modelRows, modelId)]));
    const bootstrap = patientBootstrap(modelRows, modelIds, BOOTSTRAP_REPLICATES, 2026111800 + signatureIndex);
    const deltaPoint = pointAucs.signature_plus_baselines - pointAucs.baseline_combined;
    const deltaInterval = bootstrapInterval(bootstrap.deltaDraws);
    let incrementalStatus = 'uncertain';
    if (Number.isFinite(deltaInterval.low) && deltaInterval.low >= MIN_MEANINGFUL_DELTA_AUC) incrementalStatus = 'supported';
    else if (Number.isFinite(deltaInterval.high) && deltaInterval.high < MIN_MEANINGFUL_DELTA_AUC) incrementalStatus = 'not_supported';

    const definitions = [
      { id: 'signature_only', label: `${signatureId} fixed positive score`, features: [signatureId] },
      ...baselineFeatures.map((feature) => ({ id: `${feature}_alone`, label: `${feature} fixed positive score`, features: [feature] })),
      { id: 'baseline_combined', label: 'Four simple biology baselines combined', features: baselineFeatures },
      { id: 'signature_plus_baselines', label: `${signatureId} plus four simple biology baselines`, features: [signatureId, ...baselineFeatures] },
    ];
    for (const definition of definitions) {
      const interval = bootstrapInterval(bootstrap.aucDraws[definition.id]);
      metricRows.push({
        cohort: 'GSE154200', signature_id: signatureId, model_id: definition.id, model_label: definition.label,
        features: definition.features.join(';'), patients, samples: rows.length, feature_count: definition.features.length,
        oof_auc: formatNumber(pointAucs[definition.id]), auc_ci_low: formatNumber(interval.low), auc_ci_high: formatNumber(interval.high),
        delta_vs_baseline_combined: definition.id === 'signature_plus_baselines' ? formatNumber(deltaPoint) : '',
        delta_ci_low: definition.id === 'signature_plus_baselines' ? formatNumber(deltaInterval.low) : '',
        delta_ci_high: definition.id === 'signature_plus_baselines' ? formatNumber(deltaInterval.high) : '',
        min_meaningful_delta: definition.id === 'signature_plus_baselines' ? MIN_MEANINGFUL_DELTA_AUC : '',
        incremental_status: definition.id === 'signature_plus_baselines' ? incrementalStatus : '',
      });
    }

    incrementBySignature[signatureId] = {
      baseline_alone_auc: Object.fromEntries(baselineFeatures.map((feature) => [feature, formatNumber(pointAucs[`${feature}_alone`])])),
      baseline_combined_auc: formatNumber(pointAucs.baseline_combined),
      signature_plus_baselines_auc: formatNumber(pointAucs.signature_plus_baselines),
      delta_auc: formatNumber(deltaPoint), delta_ci: [formatNumber(deltaInterval.low), formatNumber(deltaInterval.high)],
      incremental_status: incrementalStatus,
    };

    for (const [modelId, fitted, features] of [
      ['baseline_combined', baselineCombined, baselineFeatures],
      ['signature_plus_baselines', signaturePlusBaselines, [signatureId, ...baselineFeatures]],
    ]) {
      for (const summary of coefficientSummary(fitted.foldCoefficients, features)) {
        coefficientRows.push({
          cohort: 'GSE154200', signature_id: signatureId, model_id: modelId, feature: summary.feature,
          folds: summary.folds, positive_folds: summary.positive_folds, zero_folds: summary.zero_folds,
          median_coefficient: formatNumber(summary.median_coefficient), minimum_coefficient: formatNumber(summary.minimum_coefficient),
          maximum_coefficient: formatNumber(summary.maximum_coefficient),
        });
      }
    }
    for (const row of modelRows) {
      for (const modelId of modelIds) {
        predictionRows.push({
          cohort: 'GSE154200', signature_id: signatureId, patient_id: row.patient_id, gsm: row.gsm,
          label: row.label, model_id: modelId, prediction: formatNumber(row[modelId]),
        });
      }
    }
  }

  if (metricRows.length !== PRIMARY_SIGNATURES.length * 7) throw new Error(`Expected 56 metric rows, found ${metricRows.length}.`);
  writeCsv(path.join(OUT, 'increment_model_metrics.csv'), metricRows, [
    'cohort', 'signature_id', 'model_id', 'model_label', 'features', 'patients', 'samples', 'feature_count',
    'oof_auc', 'auc_ci_low', 'auc_ci_high', 'delta_vs_baseline_combined', 'delta_ci_low', 'delta_ci_high',
    'min_meaningful_delta', 'incremental_status',
  ]);
  writeCsv(path.join(OUT, 'increment_coefficient_summary.csv'), coefficientRows, [
    'cohort', 'signature_id', 'model_id', 'feature', 'folds', 'positive_folds', 'zero_folds',
    'median_coefficient', 'minimum_coefficient', 'maximum_coefficient',
  ]);
  writeCsv(path.join(OUT, 'increment_patient_predictions_long.csv'), predictionRows, [
    'cohort', 'signature_id', 'patient_id', 'gsm', 'label', 'model_id', 'prediction',
  ]);

  const summary = {
    generatedAt: new Date().toISOString(),
    status: 'GSE154200_INDEPENDENT_D2_COMPLETED',
    scope: 'Independent D2 score and incremental benchmark in GSE154200; no manuscript or submission files changed.',
    dataQuality: {
      platform: 'GPL20301', patients, samples: rows.length, lesional: 88, nonlesional: 88,
      treatment: 'No treatment for all 176 included samples', expressionTransform: 'edgeR TMM log-CPM(prior.count=2)',
    },
    rules: {
      endpoint: 'D2: lesional skin = 1, nonlesional skin from the same patient = 0.',
      signatureScore: 'Equal-weight mean of available gene z-scores fitted once across all 176 D2 reference samples; positive direction retained.',
      partialCoverage: 'IMscore uses 4/5 genes and is retained only as a partial-coverage benchmark estimate.',
      crossValidation: 'Leave-one-patient-out with both samples from the held-out patient removed from training.',
      model: 'Nonnegative ridge logistic regression, unpenalized intercept, lambda=1.0.',
      uncertainty: `Patient-level bootstrap with ${BOOTSTRAP_REPLICATES} replicates and percentile 95% confidence intervals.`,
      meaningfulDelta: `Supported requires the lower 95% bootstrap bound of Delta AUC versus four-baseline model to be at least ${MIN_MEANINGFUL_DELTA_AUC}.`,
      interpretationLimit: 'Conditional benchmark, not exact original model reconstruction, external validation of an original model, diagnostic utility, treatment response, or clinical utility.',
    },
    baselineCoverage: Object.fromEntries(baselineDefinitions.map((definition) => [definition.id, `${coverageById[definition.id].genes_used}/${coverageById[definition.id].genes_declared}`])),
    signatureCoverage: Object.fromEntries(signatureDefinitions.map((definition) => [definition.id, `${coverageById[definition.id].genes_used}/${coverageById[definition.id].genes_declared}`])),
    directSignatureMetrics: Object.fromEntries(signatureMetrics.map((row) => [row.signature_id, {
      auc: row.auc, auc_ci: [row.auc_ci_low, row.auc_ci_high], coverage_status: row.coverage_status, direction_result: row.direction_result,
    }])),
    incrementResults: incrementBySignature,
    convergenceFailures,
    limitations: [
      'All signatures lack recovered original weights/intercepts and are represented by equal-weight positive-direction gene scores.',
      'Cell composition is not adjusted; lesion versus nonlesion differences include epithelial, immune, and stromal changes.',
      'Eight signature comparisons are exploratory benchmarks without multiplicity control across signatures.',
      'The affected-unaffected paired endpoint is not a clinical screening endpoint.',
    ],
    outputs: {
      sampleScores: 'data/benchmark/score_stage_independent/GSE154200/increment_stage/D2_sample_scores.csv',
      signatureMetrics: 'data/benchmark/score_stage_independent/GSE154200/increment_stage/D2_signature_metrics.csv',
      incrementMetrics: 'data/benchmark/score_stage_independent/GSE154200/increment_stage/increment_model_metrics.csv',
      coefficientSummary: 'data/benchmark/score_stage_independent/GSE154200/increment_stage/increment_coefficient_summary.csv',
      predictions: 'data/benchmark/score_stage_independent/GSE154200/increment_stage/increment_patient_predictions_long.csv',
    },
    runtimeSeconds: (Date.now() - startedAt) / 1000,
  };
  fs.mkdirSync(OUT, { recursive: true });
  fs.writeFileSync(path.join(OUT, 'GSE154200_benchmark_summary.json'), `${JSON.stringify(summary, null, 2)}\n`, 'utf8');
  console.log(JSON.stringify(summary, null, 2));
}

try {
  main();
} catch (error) {
  console.error(error);
  process.exitCode = 1;
}
