import fs from 'node:fs';
import path from 'node:path';
import zlib from 'node:zlib';
import readline from 'node:readline';
import { fileURLToPath } from 'node:url';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const BENCH = path.join(ROOT, 'data', 'benchmark');
const SCORE = path.join(BENCH, 'score_stage');
const OUT = path.join(SCORE, 'increment_stage');
const LAMBDA = 1.0;
const BOOTSTRAP_REPLICATES = 2000;
const MIN_MEANINGFUL_DELTA_AUC = 0.03;

const COHORTS = [
  { id: 'GSE41664', platform: 'GPL570', matrix: 'sources/raw/GSE41664/GSE41664_series_matrix.txt.gz', annotation: 'sources/raw/GPL570/GPL570.annot.gz' },
  { id: 'GSE6710', platform: 'GPL96', matrix: 'sources/raw/GSE6710/GSE6710_series_matrix.txt.gz', annotation: 'sources/raw/GPL96/GPL96_full.txt' },
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
  fs.writeFileSync(file, `${[columns.join(','), ...rows.map((row) => columns.map((column) => escape(row[column])).join(','))].join('\n')}\n`, 'utf8');
}

function openLines(relativePath) {
  const source = fs.createReadStream(path.join(ROOT, relativePath));
  const input = relativePath.endsWith('.gz') ? source.pipe(zlib.createGunzip()) : source;
  return { reader: readline.createInterface({ input, crlfDelay: Infinity }), input };
}

function tabFields(line) {
  return line.split('\t').map((value) => {
    let clean = value.trim();
    if (clean.startsWith('"') && clean.endsWith('"')) clean = clean.slice(1, -1);
    return clean.replaceAll('""', '"');
  });
}

function finiteNumber(value) {
  if (value == null || value === '') return null;
  const number = Number(value);
  return Number.isFinite(number) ? number : null;
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
  const fields = line.split('\t');
  return { term, url: fields[1], genes: [...new Set(fields.slice(2).filter(Boolean))] };
}

async function parseAnnotation(relativePath, targetGenes) {
  const { reader, input } = openLines(relativePath);
  const bySymbol = new Map();
  let inTable = false;
  let headerSeen = false;
  let symbolIndex = -1;
  try {
    for await (const rawLine of reader) {
      const line = rawLine.replace(/^\uFEFF/, '');
      if (!inTable) {
        if (line === '!platform_table_begin') inTable = true;
        continue;
      }
      if (!headerSeen) {
        const headers = line.split('\t');
        symbolIndex = headers.findIndex((header) => /^gene[ _]symbol$/i.test(header.trim()));
        if (symbolIndex < 0) throw new Error(`No Gene symbol column in ${relativePath}`);
        headerSeen = true;
        continue;
      }
      if (line === '!platform_table_end') break;
      const fields = tabFields(line);
      const probe = fields[0];
      if (!probe) continue;
      const symbols = fields[symbolIndex].split('///').map((symbol) => symbol.trim()).filter((symbol) => targetGenes.has(symbol));
      for (const symbol of new Set(symbols)) {
        if (!bySymbol.has(symbol)) bySymbol.set(symbol, new Set());
        bySymbol.get(symbol).add(probe);
      }
    }
  } finally {
    input.destroy();
  }
  return bySymbol;
}

async function parseExpression(relativePath, wantedProbes) {
  const { reader, input } = openLines(relativePath);
  const rows = new Map();
  let header = null;
  let inTable = false;
  try {
    for await (const rawLine of reader) {
      const line = rawLine.replace(/^\uFEFF/, '');
      if (!inTable) {
        if (line === '!series_matrix_table_begin') inTable = true;
        continue;
      }
      if (line === '!series_matrix_table_end') break;
      if (header == null) {
        header = tabFields(line).slice(1);
        continue;
      }
      const fields = tabFields(line);
      if (!wantedProbes.has(fields[0])) continue;
      rows.set(fields[0], fields.slice(1).map(finiteNumber));
    }
  } finally {
    input.destroy();
  }
  if (!header) throw new Error(`No expression header in ${relativePath}`);
  return { header, rows };
}

function collapseProbes(annotation, expression, genes, referenceIndices) {
  const values = new Map();
  for (const gene of genes) {
    const candidates = [...(annotation.get(gene) || [])].map((probe) => {
      const expressionValues = expression.rows.get(probe);
      if (!expressionValues) return null;
      const reference = referenceIndices.map((index) => expressionValues[index]);
      const finite = reference.filter(Number.isFinite).length;
      return { probe, expressionValues, finite, coverage: finite / referenceIndices.length, mean: finiteMean(reference) };
    }).filter(Boolean);
    const eligible = candidates.filter((candidate) => candidate.coverage >= 0.95 && Number.isFinite(candidate.mean));
    const pool = eligible.length ? eligible : candidates;
    pool.sort((a, b) => b.mean - a.mean || b.coverage - a.coverage || a.probe.localeCompare(b.probe));
    if (pool[0]) values.set(gene, pool[0].expressionValues);
  }
  return values;
}

function fitGeneStats(genes, geneValues, trainIndices) {
  const stats = new Map();
  for (const gene of genes) {
    const values = geneValues.get(gene);
    if (!values) continue;
    const trainValues = trainIndices.map((index) => values[index]).filter(Number.isFinite);
    const mean = finiteMean(trainValues);
    const sd = finiteSd(trainValues);
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
  if (xRows.length !== labels.length || xRows.length === 0) throw new Error('Inconsistent logistic training data.');
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
  return {
    low: percentile(values, 0.025),
    high: percentile(values, 0.975),
    replicates: values.length,
  };
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

function validateD2Rows(rows, expectedPairs, cohort) {
  const byPatient = new Map();
  for (const row of rows) {
    if (!byPatient.has(row.patient_id)) byPatient.set(row.patient_id, []);
    byPatient.get(row.patient_id).push(row);
  }
  if (byPatient.size !== expectedPairs) throw new Error(`${cohort}: expected ${expectedPairs} D2 patients, found ${byPatient.size}.`);
  for (const [patientId, patientRows] of byPatient) {
    const labels = new Set(patientRows.map((row) => row.label));
    if (patientRows.length !== 2 || !labels.has(0) || !labels.has(1)) {
      throw new Error(`${cohort}: patient ${patientId} does not have exactly one lesional and one nonlesional sample.`);
    }
  }
}

async function evaluateCohort(cohort) {
  const expectedPairs = cohort.id === 'GSE41664' ? 53 : 13;
  const primarySignatures = cohort.id === 'GSE41664'
    ? ['RiboSisRG', 'TRGS', 'Apoptosis5', 'm6A10', 'PLOS14']
    : ['TRGS', 'IMscore', 'Ferroptosis8', 'Apoptosis5', 'm6A10'];
  const metadataPath = `data/benchmark/score_stage/${cohort.id}/sample_metadata_and_scores.csv`;
  const metadataRows = readCsv(metadataPath);
  const d2Rows = metadataRows
    .filter((row) => row.analysis_role === 'D2_lesional' || row.analysis_role === 'D2_nonlesional')
    .map((row) => ({ ...row, label: row.tissue === 'lesional' ? 1 : 0 }));
  validateD2Rows(d2Rows, expectedPairs, cohort.id);

  const geneSetDefinitions = GENE_SETS.map((definition) => ({ ...definition, ...readGmtTerm(definition.file, definition.term) }));
  const targetGenes = new Set(geneSetDefinitions.flatMap((definition) => definition.genes));
  const annotation = await parseAnnotation(cohort.annotation, targetGenes);
  const wantedProbes = new Set([...annotation.values()].flatMap((probes) => [...probes]));
  const expression = await parseExpression(cohort.matrix, wantedProbes);
  const expressionIndexByGsm = new Map(expression.header.map((gsm, index) => [gsm, index]));
  const referenceIndices = d2Rows.map((row) => {
    const index = expressionIndexByGsm.get(row.gsm);
    if (index == null) throw new Error(`${cohort.id}: ${row.gsm} is absent from the expression matrix.`);
    return index;
  });
  const platformGenes = new Set(geneSetDefinitions.flatMap((definition) => definition.genes.filter((gene) => annotation.has(gene))));
  const collapsed = collapseProbes(annotation, expression, platformGenes, referenceIndices);
  const enrichedRows = d2Rows.map((row, index) => ({ ...row, expression_index: referenceIndices[index] }));

  const coverageRows = [];
  for (const definition of geneSetDefinitions) {
    const genesOnPlatform = definition.genes.filter((gene) => annotation.has(gene));
    const genesUsed = genesOnPlatform.filter((gene) => collapsed.has(gene));
    const stats = fitGeneStats(genesUsed, collapsed, referenceIndices);
    const scored = scoreGeneSet(genesUsed, collapsed, stats, referenceIndices);
    if (genesUsed.length < 5) throw new Error(`${cohort.id}: baseline ${definition.id} has fewer than five usable genes.`);
    enrichedRows.forEach((row, index) => { row[definition.id] = scored[index]; });
    coverageRows.push({
      cohort: cohort.id,
      gene_set_id: definition.id,
      term: definition.term,
      source_file: definition.file,
      source_url: definition.url,
      genes_in_set: definition.genes.length,
      genes_on_platform: genesOnPlatform.length,
      genes_used: genesUsed.length,
      coverage_percent: definition.genes.length ? (genesUsed.length / definition.genes.length) * 100 : 0,
      finite_scores: scored.filter(Number.isFinite).length,
      missing_genes: definition.genes.filter((gene) => !collapsed.has(gene)).join(';'),
    });
  }

  for (const signature of primarySignatures) {
    for (const row of enrichedRows) {
      if (!Number.isFinite(Number(row[signature]))) throw new Error(`${cohort.id}: missing score for ${signature} in ${row.gsm}.`);
    }
  }
  for (const row of enrichedRows) {
    for (const definition of geneSetDefinitions) {
      if (!Number.isFinite(Number(row[definition.id]))) throw new Error(`${cohort.id}: missing baseline ${definition.id} in ${row.gsm}.`);
    }
  }

  const baselineFeatures = geneSetDefinitions.map((definition) => definition.id);
  const modelIds = ['signature_only', ...baselineFeatures.map((feature) => `${feature}_alone`), 'baseline_combined', 'signature_plus_baselines'];
  const predictionRows = [];
  const metricRows = [];
  const coefficientRows = [];
  const summaryBySignature = {};
  let convergenceFailures = 0;
  for (const signature of primarySignatures) {
    const baselineCombined = crossValidateByPatient(enrichedRows, baselineFeatures);
    const signaturePlusBaselines = crossValidateByPatient(enrichedRows, [signature, ...baselineFeatures]);
    convergenceFailures += [...baselineCombined.foldCoefficients, ...signaturePlusBaselines.foldCoefficients]
      .filter((fold) => !fold.converged).length;
    const baselinePredictionByGsm = new Map(baselineCombined.predictions.map((row) => [row.gsm, row.prediction]));
    const combinedPredictionByGsm = new Map(signaturePlusBaselines.predictions.map((row) => [row.gsm, row.prediction]));
    const signatureRowsForBootstrap = enrichedRows.map((row) => {
      const scoredRow = {
        cohort: cohort.id,
        signature_id: signature,
        patient_id: row.patient_id,
        gsm: row.gsm,
        label: row.label,
        signature_only: Number(row[signature]),
        baseline_combined: baselinePredictionByGsm.get(row.gsm),
        signature_plus_baselines: combinedPredictionByGsm.get(row.gsm),
      };
      for (const feature of baselineFeatures) scoredRow[`${feature}_alone`] = Number(row[feature]);
      return scoredRow;
    });
    for (const row of signatureRowsForBootstrap) {
      for (const modelId of modelIds) {
        if (!Number.isFinite(row[modelId])) throw new Error(`${cohort.id}: incomplete prediction for ${signature}, ${row.gsm}, ${modelId}.`);
      }
    }

    const bootstrap = patientBootstrap(signatureRowsForBootstrap, modelIds, BOOTSTRAP_REPLICATES, 20261118 + primarySignatures.indexOf(signature));
    const pointAucs = Object.fromEntries(modelIds.map((modelId) => [modelId, aucFromRows(signatureRowsForBootstrap, modelId)]));
    const deltaPoint = pointAucs.signature_plus_baselines - pointAucs.baseline_combined;
    const deltaInterval = bootstrapInterval(bootstrap.deltaDraws);
    let incrementalStatus = 'uncertain';
    if (Number.isFinite(deltaInterval.low) && deltaInterval.low >= MIN_MEANINGFUL_DELTA_AUC) incrementalStatus = 'supported';
    else if (Number.isFinite(deltaInterval.high) && deltaInterval.high < MIN_MEANINGFUL_DELTA_AUC) incrementalStatus = 'not_supported';

    const definitions = [
      { id: 'signature_only', label: `${signature} fixed positive score`, features: [signature] },
      ...baselineFeatures.map((feature) => ({ id: `${feature}_alone`, label: `${feature} fixed positive score`, features: [feature] })),
      { id: 'baseline_combined', label: 'Four simple biology baselines combined', features: baselineFeatures },
      { id: 'signature_plus_baselines', label: `${signature} plus four simple biology baselines`, features: [signature, ...baselineFeatures] },
    ];
    for (const definition of definitions) {
      const interval = bootstrapInterval(bootstrap.aucDraws[definition.id]);
      metricRows.push({
        cohort: cohort.id,
        signature_id: signature,
        model_id: definition.id,
        model_label: definition.label,
        features: definition.features.join(';'),
        patients: expectedPairs,
        samples: enrichedRows.length,
        feature_count: definition.features.length,
        oof_auc: formatNumber(pointAucs[definition.id]),
        auc_ci_low: formatNumber(interval.low),
        auc_ci_high: formatNumber(interval.high),
        delta_vs_baseline_combined: definition.id === 'signature_plus_baselines' ? formatNumber(deltaPoint) : '',
        delta_ci_low: definition.id === 'signature_plus_baselines' ? formatNumber(deltaInterval.low) : '',
        delta_ci_high: definition.id === 'signature_plus_baselines' ? formatNumber(deltaInterval.high) : '',
        min_meaningful_delta: definition.id === 'signature_plus_baselines' ? MIN_MEANINGFUL_DELTA_AUC : '',
        incremental_status: definition.id === 'signature_plus_baselines' ? incrementalStatus : '',
      });
    }

    summaryBySignature[signature] = {
      baseline_alone_auc: Object.fromEntries(baselineFeatures.map((feature) => [feature, formatNumber(pointAucs[`${feature}_alone`])])),
      baseline_combined_auc: formatNumber(pointAucs.baseline_combined),
      signature_plus_baselines_auc: formatNumber(pointAucs.signature_plus_baselines),
      delta_auc: formatNumber(deltaPoint),
      delta_ci: [formatNumber(deltaInterval.low), formatNumber(deltaInterval.high)],
      incremental_status: incrementalStatus,
    };

    for (const [modelId, fitted, features] of [
      ['baseline_combined', baselineCombined, baselineFeatures],
      ['signature_plus_baselines', signaturePlusBaselines, [signature, ...baselineFeatures]],
    ]) {
      for (const summary of coefficientSummary(fitted.foldCoefficients, features)) {
        coefficientRows.push({
          cohort: cohort.id,
          signature_id: signature,
          model_id: modelId,
          feature: summary.feature,
          folds: summary.folds,
          positive_folds: summary.positive_folds,
          zero_folds: summary.zero_folds,
          median_coefficient: formatNumber(summary.median_coefficient),
          minimum_coefficient: formatNumber(summary.minimum_coefficient),
          maximum_coefficient: formatNumber(summary.maximum_coefficient),
        });
      }
    }
    for (const row of signatureRowsForBootstrap) {
      for (const modelId of modelIds) {
        predictionRows.push({
          cohort: cohort.id,
          signature_id: signature,
          patient_id: row.patient_id,
          gsm: row.gsm,
          label: row.label,
          model_id: modelId,
          prediction: formatNumber(row[modelId]),
        });
      }
    }
  }

  return { coverageRows, metricRows, coefficientRows, predictionRows, summaryBySignature, convergenceFailures };
}

async function main() {
  const startedAt = Date.now();
  const cohortResults = [];
  for (const cohort of COHORTS) cohortResults.push({ cohort, result: await evaluateCohort(cohort) });
  fs.mkdirSync(OUT, { recursive: true });
  const coverageRows = cohortResults.flatMap((entry) => entry.result.coverageRows);
  const metricRows = cohortResults.flatMap((entry) => entry.result.metricRows);
  const coefficientRows = cohortResults.flatMap((entry) => entry.result.coefficientRows);
  const predictionRows = cohortResults.flatMap((entry) => entry.result.predictionRows);
  writeCsv(path.join(OUT, 'baseline_gene_coverage.csv'), coverageRows, [
    'cohort', 'gene_set_id', 'term', 'source_file', 'source_url', 'genes_in_set', 'genes_on_platform',
    'genes_used', 'coverage_percent', 'finite_scores', 'missing_genes',
  ]);
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
  const expectedMetricRows = 2 * 5 * (4 + 3);
  if (metricRows.length !== expectedMetricRows) throw new Error(`Expected ${expectedMetricRows} metric rows, found ${metricRows.length}.`);
  const summary = {
    generatedAt: new Date().toISOString(),
    status: 'INCREMENT_STAGE_COMPLETED',
    scope: 'D2 incremental value of prespecified positive-direction signature scores over simple GO/Hallmark biology baselines; no manuscript or submission files changed.',
    rules: {
      endpoint: 'D2: lesional skin = 1, nonlesional skin from the same patient = 0.',
      signatureComparator: 'Published gene lists scored as the equal-weight mean of available gene z-scores; positive direction is preserved and not reversed after seeing results.',
      baselineComparator: 'Four prespecified simple scores: GOBP_KERATINOCYTE_DIFFERENTIATION, HALLMARK_INFLAMMATORY_RESPONSE, HALLMARK_INTERFERON_GAMMA_RESPONSE, and HALLMARK_E2F_TARGETS. Each is the equal-weight mean of available gene z-scores.',
      crossValidation: 'Leave-one-patient-out; both lesional and nonlesional samples from the held-out patient are excluded from training.',
      model: 'Ridge logistic regression on standardized features with an unpenalized intercept, lambda=1.0, and nonnegative feature coefficients to preserve prespecified positive directions.',
      uncertainty: `Patient-level bootstrap with ${BOOTSTRAP_REPLICATES} replicates and percentile 95% confidence intervals.`,
      meaningfulDelta: `Delta AUC versus the four-baseline model is considered supported only when the bootstrap 95% interval lower bound is at least ${MIN_MEANINGFUL_DELTA_AUC}.`,
      interpretationLimit: 'These are conditional cross-cohort benchmark estimates, not external validation, diagnostic utility, or treatment-response evidence.',
    },
    sources: {
      msigdbVersion: '2025.1.Hs',
      msigdbFiles: GENE_SETS.map((definition) => definition.file).filter((value, index, all) => all.indexOf(value) === index),
    },
    cohorts: Object.fromEntries(cohortResults.map(({ cohort, result }) => [cohort.id, {
      pairs: cohort.id === 'GSE41664' ? 53 : 13,
      samples: cohort.id === 'GSE41664' ? 106 : 26,
      primarySignatures: Object.keys(result.summaryBySignature),
      baselineCoverage: Object.fromEntries(result.coverageRows.map((row) => [row.gene_set_id, `${row.genes_used}/${row.genes_in_set}`])),
      convergenceFailures: result.convergenceFailures,
      results: result.summaryBySignature,
    }])),
    outputs: {
      modelMetrics: 'data/benchmark/score_stage/increment_stage/increment_model_metrics.csv',
      coefficientSummary: 'data/benchmark/score_stage/increment_stage/increment_coefficient_summary.csv',
      predictions: 'data/benchmark/score_stage/increment_stage/increment_patient_predictions_long.csv',
      baselineCoverage: 'data/benchmark/score_stage/increment_stage/baseline_gene_coverage.csv',
    },
    runtimeSeconds: (Date.now() - startedAt) / 1000,
  };
  fs.writeFileSync(path.join(OUT, 'increment_run_summary.json'), `${JSON.stringify(summary, null, 2)}\n`, 'utf8');
  console.log(JSON.stringify(summary, null, 2));
}

main().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});
