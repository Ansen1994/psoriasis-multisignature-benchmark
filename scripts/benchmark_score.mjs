import fs from 'node:fs';
import path from 'node:path';
import zlib from 'node:zlib';
import readline from 'node:readline';
import { fileURLToPath } from 'node:url';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const BENCH = path.join(ROOT, 'data', 'benchmark');
const OUT = path.join(BENCH, 'score_stage');

const COHORTS = [
  {
    id: 'GSE41664',
    platform: 'GPL570',
    matrix: 'sources/raw/GSE41664/GSE41664_series_matrix.txt.gz',
    annotation: 'sources/raw/GPL570/GPL570.annot.gz',
  },
  {
    id: 'GSE6710',
    platform: 'GPL96',
    matrix: 'sources/raw/GSE6710/GSE6710_series_matrix.txt.gz',
    annotation: 'sources/raw/GPL96/GPL96_full.txt',
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

function openLines(relativePath) {
  const source = fs.createReadStream(path.join(ROOT, relativePath));
  const input = relativePath.endsWith('.gz') ? source.pipe(zlib.createGunzip()) : source;
  const reader = readline.createInterface({ input, crlfDelay: Infinity });
  return { reader, input };
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
    if (Number.isFinite(value)) {
      sum += value;
      count += 1;
    }
  }
  return count > 0 ? sum / count : null;
}

function finiteSd(values) {
  const mean = finiteMean(values);
  if (mean == null) return null;
  let sum = 0;
  let count = 0;
  for (const value of values) {
    if (Number.isFinite(value)) {
      sum += (value - mean) ** 2;
      count += 1;
    }
  }
  return count > 1 ? Math.sqrt(sum / (count - 1)) : null;
}

function median(values) {
  const sorted = values.filter(Number.isFinite).sort((a, b) => a - b);
  if (sorted.length === 0) return null;
  const middle = Math.floor(sorted.length / 2);
  return sorted.length % 2 ? sorted[middle] : (sorted[middle - 1] + sorted[middle]) / 2;
}

function slug(value) {
  return String(value).trim().toLowerCase().replace(/[^a-z0-9]+/g, '_').replace(/^_+|_+$/g, '');
}

async function parseAnnotation(relativePath, targetGenes) {
  const { reader, input } = openLines(relativePath);
  const bySymbol = new Map();
  let inTable = false;
  let headerSeen = false;
  let symbolIndex = -1;
  let rows = 0;
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
      const probe = fields[0]?.trim();
      if (!probe) continue;
      rows += 1;
      const symbols = new Set((fields[symbolIndex] || '')
        .split('///')
        .map((symbol) => symbol.trim())
        .filter((symbol) => targetGenes.has(symbol)));
      for (const symbol of symbols) {
        if (!bySymbol.has(symbol)) bySymbol.set(symbol, new Set());
        bySymbol.get(symbol).add(probe);
      }
    }
  } finally {
    input.destroy();
  }
  return { bySymbol, rows };
}

async function parseSeriesMetadata(relativePath) {
  const { reader, input } = openLines(relativePath);
  let titles = [];
  let gsmIds = [];
  let sources = [];
  const characteristicRows = [];
  try {
    for await (const rawLine of reader) {
      const line = rawLine.replace(/^\uFEFF/, '');
      if (line.startsWith('!Sample_title\t')) titles = tabFields(line).slice(1);
      else if (line.startsWith('!Sample_geo_accession\t')) gsmIds = tabFields(line).slice(1);
      else if (line.startsWith('!Sample_source_name_ch1\t')) sources = tabFields(line).slice(1);
      else if (line.startsWith('!Sample_characteristics_ch1\t')) characteristicRows.push(tabFields(line).slice(1));
      else if (line === '!series_matrix_table_begin') break;
    }
  } finally {
    input.destroy();
  }
  if (titles.length === 0 || gsmIds.length === 0 || titles.length !== gsmIds.length) {
    throw new Error(`Incomplete sample metadata in ${relativePath}`);
  }
  return gsmIds.map((gsm, sampleIndex) => {
    const characteristics = {};
    characteristicRows.forEach((values, rowIndex) => {
      const value = (values[sampleIndex] || '').trim();
      if (!value) return;
      const match = value.match(/^([^:]+):\s*(.*)$/);
      if (match) characteristics[slug(match[1])] = match[2].trim();
      else characteristics[`characteristic_${rowIndex + 1}`] = value;
    });
    return {
      gsm,
      title: titles[sampleIndex] || '',
      source_name: sources[sampleIndex] || '',
      characteristics,
    };
  });
}

async function parseExpression(relativePath, metadata, wantedProbes) {
  const { reader, input } = openLines(relativePath);
  const rows = new Map();
  let inTable = false;
  let header = null;
  let dataRows = 0;
  let selectedRows = 0;
  let missingCells = 0;
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
      dataRows += 1;
      const fields = tabFields(line);
      const probe = fields[0]?.trim();
      if (!wantedProbes.has(probe)) continue;
      selectedRows += 1;
      const values = fields.slice(1).map((value) => {
        const parsed = finiteNumber(value.trim());
        if (parsed == null) missingCells += 1;
        return parsed;
      });
      rows.set(probe, values);
    }
  } finally {
    input.destroy();
  }
  if (header == null) throw new Error(`No expression header in ${relativePath}`);
  if (header.length !== metadata.length) {
    throw new Error(`Sample count mismatch in ${relativePath}: header=${header.length}, metadata=${metadata.length}`);
  }
  const metadataByGsm = new Map(metadata.map((sample) => [sample.gsm, sample]));
  const missingMetadata = header.filter((gsm) => !metadataByGsm.has(gsm));
  if (missingMetadata.length > 0) throw new Error(`Missing metadata for ${missingMetadata.join(',')}`);
  return { header, rows, dataRows, selectedRows, missingCells };
}

function collapseProbes(annotation, expression, genes, referenceSampleIndices) {
  const geneValues = new Map();
  const manifest = [];
  for (const gene of genes) {
    const probes = [...(annotation.bySymbol.get(gene) || [])].sort();
    const candidates = probes.map((probe) => {
      const values = expression.rows.get(probe) || null;
      const referenceValues = values ? referenceSampleIndices.map((index) => values[index]) : [];
      const finiteCount = referenceValues.filter(Number.isFinite).length;
      return {
        probe,
        values,
        finiteCount,
        coverage: referenceSampleIndices.length > 0 ? finiteCount / referenceSampleIndices.length : 0,
        mean: finiteMean(referenceValues),
      };
    }).filter((candidate) => candidate.values);
    const eligible = candidates.filter((candidate) => candidate.coverage >= 0.95 && candidate.mean != null);
    const pool = eligible.length > 0 ? eligible : candidates;
    pool.sort((a, b) => {
      if (b.mean !== a.mean) return b.mean - a.mean;
      if (b.coverage !== a.coverage) return b.coverage - a.coverage;
      return a.probe.localeCompare(b.probe);
    });
    const selected = pool[0] || null;
    if (selected) geneValues.set(gene, selected.values);
    manifest.push({
      gene,
      annotation_probe_count: probes.length,
      expression_probe_count: candidates.length,
      selected_probe: selected?.probe || '',
      selected_mean: selected?.mean == null ? '' : selected.mean,
      selected_finite_fraction: selected ? selected.coverage.toFixed(4) : '0.0000',
      reference_sample_count: referenceSampleIndices.length,
      selection_rule: eligible.length > 0 ? 'highest_mean_ge_95pct_complete' : 'fallback_highest_mean_any_coverage',
    });
  }
  return { geneValues, manifest };
}

function zScoreFromReference(values, referenceValues) {
  const mean = finiteMean(referenceValues);
  const sd = finiteSd(referenceValues);
  const transformed = values.map((value) => {
    if (!Number.isFinite(value) || mean == null || sd == null || sd === 0) return null;
    return (value - mean) / sd;
  });
  return { transformed, mean, sd, referenceFiniteCount: referenceValues.filter(Number.isFinite).length };
}

function deriveGse41664Sample(sample) {
  const title = sample.title;
  let patientId = '';
  let tissue = '';
  let timepoint = '';
  let treatmentContext = '';
  let sourceBlock = '';
  let match = title.match(/^Skin-Psoriasis-(\d+)-(LS|NL)$/i);
  if (match) {
    patientId = `skin_psoriasis_${match[1]}`;
    tissue = match[2].toUpperCase() === 'LS' ? 'lesional' : 'nonlesional';
    timepoint = 'baseline';
    treatmentContext = 'untreated_or_unspecified';
    sourceBlock = 'skin_psoriasis_title_encoded';
  } else if ((match = title.match(/^(\d+)\s+(lesional|non-lesional)$/i))) {
    patientId = `numeric_${match[1]}`;
    tissue = match[2].toLowerCase() === 'lesional' ? 'lesional' : 'nonlesional';
    timepoint = 'baseline';
    treatmentContext = 'untreated_or_unspecified';
    sourceBlock = 'numeric_id_title_encoded';
  } else if ((match = title.match(/^PsOEt-(\d+)-(NL|LS|Week\s*\d+)$/i))) {
    patientId = `psoet_${match[1]}`;
    const suffix = match[2];
    if (/^NL$/i.test(suffix)) tissue = 'nonlesional';
    else tissue = 'lesional';
    if (/^Week/i.test(suffix)) {
      timepoint = suffix.replace(/\s+/g, '');
      treatmentContext = 'etanercept';
    } else {
      timepoint = 'Pre-dose';
      treatmentContext = 'etanercept';
    }
    sourceBlock = 'psoet_etanercept_trial';
  }
  return { ...sample, patient_id: patientId, tissue, timepoint, treatment_context: treatmentContext, source_block: sourceBlock };
}

function deriveGse6710Sample(sample) {
  const title = sample.title;
  const match = title.match(/^Patient\s+(.+?),\s*(lesional|uninvolved)\s+skin$/i);
  return {
    ...sample,
    patient_id: sample.characteristics.patient_identifyier ? `patient_${sample.characteristics.patient_identifyier}` : `patient_${slug(match?.[1] || '')}`,
    tissue: match ? (match[2].toLowerCase() === 'lesional' ? 'lesional' : 'nonlesional') : '',
    timepoint: 'baseline',
    treatment_context: 'untreated_or_unspecified',
    source_block: 'paired_biopsies',
  };
}

function assignD2Pairs(samples, cohort) {
  const groups = new Map();
  for (const sample of samples) {
    if (!sample.patient_id || !['lesional', 'nonlesional'].includes(sample.tissue)) continue;
    if (!/pre-?dose|baseline/i.test(sample.timepoint) && !['untreated_or_unspecified'].includes(sample.treatment_context)) continue;
    const key = `${sample.patient_id}|${sample.treatment_context}|${sample.timepoint}`;
    if (!groups.has(key)) groups.set(key, []);
    groups.get(key).push(sample);
  }
  let pairNumber = 0;
  for (const key of [...groups.keys()].sort()) {
    const group = groups.get(key);
    const lesional = group.find((sample) => sample.tissue === 'lesional');
    const nonlesional = group.find((sample) => sample.tissue === 'nonlesional');
    if (!lesional || !nonlesional) continue;
    pairNumber += 1;
    const pairId = `${cohort}_D2_${String(pairNumber).padStart(3, '0')}`;
    lesional.pair_id = pairId;
    nonlesional.pair_id = pairId;
    lesional.analysis_role = 'D2_lesional';
    nonlesional.analysis_role = 'D2_nonlesional';
  }
  return pairNumber;
}

function assignTreatmentPairs(samples, cohort) {
  const baselineByPatient = new Map(samples
    .filter((sample) => sample.tissue === 'lesional' && /pre-?dose|baseline/i.test(sample.timepoint))
    .map((sample) => [sample.patient_id, sample]));
  let count = 0;
  for (const sample of samples) {
    if (!/^Week\d+/i.test(sample.timepoint) || sample.tissue !== 'lesional') continue;
    const baseline = baselineByPatient.get(sample.patient_id);
    if (!baseline) continue;
    count += 1;
    sample.analysis_role = 'T_posttreatment_lesional';
    sample.t_baseline_sample_id = baseline.gsm;
    sample.t_pair_id = `${cohort}_T_${String(count).padStart(3, '0')}`;
  }
  return count;
}

function aucFromLabels(rows, scoreKey) {
  const usable = rows.filter((row) => Number.isFinite(row[scoreKey]));
  const positives = usable.filter((row) => row.label === 1);
  const negatives = usable.filter((row) => row.label === 0);
  if (positives.length === 0 || negatives.length === 0) return null;
  const ordered = usable.map((row) => ({ score: row[scoreKey], label: row.label })).sort((a, b) => a.score - b.score);
  let index = 0;
  while (index < ordered.length) {
    let end = index + 1;
    while (end < ordered.length && ordered[end].score === ordered[index].score) end += 1;
    const averageRank = ((index + 1) + end) / 2;
    for (let rankIndex = index; rankIndex < end; rankIndex += 1) ordered[rankIndex].rank = averageRank;
    index = end;
  }
  const rankSumPositive = ordered.filter((row) => row.label === 1).reduce((sum, row) => sum + row.rank, 0);
  return (rankSumPositive - (positives.length * (positives.length + 1)) / 2) / (positives.length * negatives.length);
}

function averageRanks(values) {
  const indexed = values
    .map((value, index) => ({ value, index }))
    .filter((entry) => Number.isFinite(entry.value))
    .sort((a, b) => a.value - b.value || a.index - b.index);
  const ranks = new Array(values.length).fill(null);
  let start = 0;
  while (start < indexed.length) {
    let end = start + 1;
    while (end < indexed.length && indexed[end].value === indexed[start].value) end += 1;
    const rank = ((start + 1) + end) / 2;
    for (let index = start; index < end; index += 1) ranks[indexed[index].index] = rank;
    start = end;
  }
  return ranks;
}

function spearmanCorrelation(x, y) {
  const paired = x.map((value, index) => [value, y[index]]).filter(([a, b]) => Number.isFinite(a) && Number.isFinite(b));
  if (paired.length < 3) return null;
  const xRanks = averageRanks(paired.map(([a]) => a));
  const yRanks = averageRanks(paired.map(([, b]) => b));
  const xMean = finiteMean(xRanks);
  const yMean = finiteMean(yRanks);
  let numerator = 0;
  let xSum = 0;
  let ySum = 0;
  for (let index = 0; index < paired.length; index += 1) {
    const xDelta = xRanks[index] - xMean;
    const yDelta = yRanks[index] - yMean;
    numerator += xDelta * yDelta;
    xSum += xDelta ** 2;
    ySum += yDelta ** 2;
  }
  return xSum > 0 && ySum > 0 ? numerator / Math.sqrt(xSum * ySum) : null;
}

function symmetricEigenvalues(matrix) {
  const size = matrix.length;
  const work = matrix.map((row) => [...row]);
  for (let sweep = 0; sweep < 100; sweep += 1) {
    let maxOffDiagonal = 0;
    let p = 0;
    let q = 1;
    for (let row = 0; row < size; row += 1) {
      for (let column = row + 1; column < size; column += 1) {
        const magnitude = Math.abs(work[row][column]);
        if (magnitude > maxOffDiagonal) {
          maxOffDiagonal = magnitude;
          p = row;
          q = column;
        }
      }
    }
    if (maxOffDiagonal < 1e-12) break;
    const angle = 0.5 * Math.atan2(2 * work[p][q], work[q][q] - work[p][p]);
    const cosine = Math.cos(angle);
    const sine = Math.sin(angle);
    const app = work[p][p];
    const aqq = work[q][q];
    const apq = work[p][q];
    work[p][p] = cosine ** 2 * app - 2 * sine * cosine * apq + sine ** 2 * aqq;
    work[q][q] = sine ** 2 * app + 2 * sine * cosine * apq + cosine ** 2 * aqq;
    work[p][q] = 0;
    work[q][p] = 0;
    for (let index = 0; index < size; index += 1) {
      if (index === p || index === q) continue;
      const aip = work[index][p];
      const aiq = work[index][q];
      work[index][p] = cosine * aip - sine * aiq;
      work[p][index] = work[index][p];
      work[index][q] = sine * aip + cosine * aiq;
      work[q][index] = work[index][q];
    }
  }
  return work.map((row, index) => row[index]).sort((a, b) => b - a);
}

function redundancyOutputs(result) {
  const primaryIds = result.d2MetricRows
    .filter((row) => row.primary_D2_eligible === 'yes')
    .map((row) => row.signature_id);
  const pairIds = [...new Set(result.pairedLongRows.map((row) => row.pair_id))];
  const contexts = [
    {
      name: 'all_D2_reference',
      unit: 'samples',
      indices: result.d2ReferenceSampleIndices,
      values: (signatureId) => result.scores.map((score) => score[signatureId]),
    },
    {
      name: 'D2_lesional',
      unit: 'samples',
      indices: result.d2ReferenceSampleIndices.filter((index) => result.samples[index].tissue === 'lesional'),
      values: (signatureId) => result.scores.map((score) => score[signatureId]),
    },
    {
      name: 'D2_nonlesional',
      unit: 'samples',
      indices: result.d2ReferenceSampleIndices.filter((index) => result.samples[index].tissue === 'nonlesional'),
      values: (signatureId) => result.scores.map((score) => score[signatureId]),
    },
    {
      name: 'D2_paired_delta',
      unit: 'pairs',
      indices: pairIds.map((_, index) => index),
      values: (signatureId) => {
        const byPair = new Map(result.pairedLongRows
          .filter((row) => row.signature_id === signatureId)
          .map((row) => [row.pair_id, Number(row.delta_lesional_minus_nonlesional)]));
        return pairIds.map((pairId) => byPair.get(pairId));
      },
    },
  ];
  const longRows = [];
  const matrixOutputs = [];
  const dimensionRows = [];
  for (const context of contexts) {
    const vectors = new Map(primaryIds.map((signatureId) => {
      const allValues = context.values(signatureId);
      const values = context.indices.map((index) => allValues[index]).filter(Number.isFinite);
      return [signatureId, values];
    }));
    const matrix = primaryIds.map((rowId) => primaryIds.map((columnId) => (
      rowId === columnId ? 1 : spearmanCorrelation(vectors.get(rowId), vectors.get(columnId))
    )));
    const n = primaryIds.length ? vectors.get(primaryIds[0]).length : 0;
    for (let row = 0; row < primaryIds.length; row += 1) {
      for (let column = row + 1; column < primaryIds.length; column += 1) {
        longRows.push({
          cohort: result.config.id,
          context: context.name,
          unit: context.unit,
          n,
          signature_1: primaryIds[row],
          signature_2: primaryIds[column],
          spearman_rho: formatNumber(matrix[row][column]),
        });
      }
    }
    const matrixRows = primaryIds.map((signatureId, rowIndex) => {
      const row = {
        cohort: result.config.id,
        context: context.name,
        unit: context.unit,
        n,
        signature: signatureId,
      };
      primaryIds.forEach((otherId, columnIndex) => {
        row[otherId] = formatNumber(matrix[rowIndex][columnIndex]);
      });
      return row;
    });
    matrixOutputs.push({ context: context.name, matrixRows, signatureIds: primaryIds });
    const eigenvalues = primaryIds.length > 1 ? symmetricEigenvalues(matrix) : [];
    const positive = eigenvalues.filter((value) => value > 0);
    const positiveSum = positive.reduce((sum, value) => sum + value, 0);
    const squaredSum = positive.reduce((sum, value) => sum + value ** 2, 0);
    const offDiagonal = [];
    for (let row = 0; row < primaryIds.length; row += 1) {
      for (let column = row + 1; column < primaryIds.length; column += 1) {
        if (Number.isFinite(matrix[row][column])) offDiagonal.push(matrix[row][column]);
      }
    }
    dimensionRows.push({
      cohort: result.config.id,
      context: context.name,
      unit: context.unit,
      n,
      signatures: primaryIds.length,
      first_eigenvalue: formatNumber(positive[0]),
      first_eigenvalue_fraction: formatNumber(primaryIds.length ? positive[0] / primaryIds.length : null),
      participation_effective_dimension: formatNumber(squaredSum > 0 ? positiveSum ** 2 / squaredSum : null),
      participation_fraction_of_signatures: formatNumber(squaredSum > 0 ? (positiveSum ** 2 / squaredSum) / primaryIds.length : null),
      median_off_diagonal_rho: formatNumber(median(offDiagonal)),
      min_off_diagonal_rho: formatNumber(offDiagonal.length ? Math.min(...offDiagonal) : null),
      max_off_diagonal_rho: formatNumber(offDiagonal.length ? Math.max(...offDiagonal) : null),
      signatures_included: primaryIds.join(';'),
    });
  }
  return { primaryIds, longRows, matrixOutputs, dimensionRows };
}

function primaryEligibility(signature, cohort, platform, usageMap, coverageMap, selectedGenes) {
  const reasons = [];
  const compatible = signature.compatible_endpoints.split(';').includes('D2');
  const usage = usageMap.get(`${cohort}|${signature.id}`) || 'unknown';
  const coverage = coverageMap.get(`${platform}|${signature.id}`);
  const selected = selectedGenes.get(signature.id) || 0;
  if (!compatible) reasons.push('endpoint_not_D2');
  if (usage !== 'no') reasons.push(`declared_cohort_usage_${usage}`);
  if (!coverage || Number(coverage.missing_gene_count) > 0) reasons.push('platform_gene_missing');
  if (selected < signature.genes.length) reasons.push('selected_expression_gene_missing');
  return { primary: reasons.length === 0, usage, reasons: reasons.join(';') };
}

function formatNumber(value, digits = 6) {
  return Number.isFinite(value) ? Number(value).toFixed(digits) : '';
}

function numericDifference(a, b) {
  return Number.isFinite(a) && Number.isFinite(b) ? a - b : null;
}

async function processCohort(config, signatures, usageMap, coverageMap) {
  const targetGenes = new Set(signatures.flatMap((signature) => signature.genes));
  const annotation = await parseAnnotation(config.annotation, targetGenes);
  const metadata = await parseSeriesMetadata(config.matrix);
  const wantedProbes = new Set([...annotation.bySymbol.values()].flatMap((probes) => [...probes]));
  const expression = await parseExpression(config.matrix, metadata, wantedProbes);
  const metadataByGsm = new Map(metadata.map((sample) => [sample.gsm, sample]));
  const orderedMetadata = expression.header.map((gsm) => metadataByGsm.get(gsm));
  let samples = orderedMetadata.map((sample) => config.id === 'GSE41664' ? deriveGse41664Sample(sample) : deriveGse6710Sample(sample));
  const underived = samples.filter((sample) => !sample.patient_id || !sample.tissue || !sample.timepoint);
  if (underived.length > 0) throw new Error(`Unparsed sample titles in ${config.id}: ${underived.map((sample) => sample.title).join(' | ')}`);
  const d2Pairs = assignD2Pairs(samples, config.id);
  const tPostSamples = assignTreatmentPairs(samples, config.id);
  if (config.id === 'GSE41664' && (d2Pairs !== 53 || tPostSamples !== 51)) {
    throw new Error(`Unexpected GSE41664 task counts: D2 pairs=${d2Pairs}, treatment-post samples=${tPostSamples}`);
  }
  if (config.id === 'GSE6710' && d2Pairs !== 13) {
    throw new Error(`Unexpected GSE6710 D2 pair count: ${d2Pairs}`);
  }
  const d2ReferenceSampleIndices = samples
    .map((sample, index) => sample.analysis_role?.startsWith('D2_') ? index : -1)
    .filter((index) => index >= 0);
  if (d2ReferenceSampleIndices.length !== d2Pairs * 2) {
    throw new Error(`D2 reference sample count mismatch in ${config.id}: samples=${d2ReferenceSampleIndices.length}, pairs=${d2Pairs}`);
  }
  const { geneValues, manifest } = collapseProbes(annotation, expression, [...targetGenes].sort(), d2ReferenceSampleIndices);
  const zValues = new Map();
  const zParameters = [];
  for (const [gene, values] of geneValues.entries()) {
    const referenceValues = d2ReferenceSampleIndices.map((index) => values[index]);
    const fitted = zScoreFromReference(values, referenceValues);
    zValues.set(gene, fitted.transformed);
    zParameters.push({
      gene,
      d2_reference_samples: d2ReferenceSampleIndices.length,
      d2_reference_finite_values: fitted.referenceFiniteCount,
      d2_reference_mean: formatNumber(fitted.mean),
      d2_reference_sd: formatNumber(fitted.sd),
    });
  }
  const scores = samples.map(() => ({}));
  const usedGeneCounts = samples.map(() => ({}));
  for (const signature of signatures) {
    for (let sampleIndex = 0; sampleIndex < samples.length; sampleIndex += 1) {
      const values = signature.genes
        .map((gene) => zValues.get(gene)?.[sampleIndex])
        .filter(Number.isFinite);
      usedGeneCounts[sampleIndex][signature.id] = values.length;
      scores[sampleIndex][signature.id] = values.length > 0 ? values.reduce((sum, value) => sum + value, 0) / values.length : null;
    }
  }

  const selectedGenes = new Map(signatures.map((signature) => [
    signature.id,
    signature.genes.filter((gene) => geneValues.has(gene)).length,
  ]));
  const wideRows = samples.map((sample, sampleIndex) => {
    const row = {
      cohort: config.id,
      platform: config.platform,
      gsm: sample.gsm,
      title: sample.title,
      patient_id: sample.patient_id,
      tissue: sample.tissue,
      treatment_context: sample.treatment_context,
      timepoint: sample.timepoint,
      source_block: sample.source_block,
      analysis_role: sample.analysis_role || 'not_D2_pair',
      pair_id: sample.pair_id || '',
      t_pair_id: sample.t_pair_id || '',
      t_baseline_sample_id: sample.t_baseline_sample_id || '',
    };
    for (const signature of signatures) {
      row[signature.id] = formatNumber(scores[sampleIndex][signature.id]);
      row[`n_genes_${signature.id}`] = usedGeneCounts[sampleIndex][signature.id];
    }
    return row;
  });

  const rawGeneRows = [...targetGenes].sort().map((gene) => {
    const row = { gene };
    samples.forEach((sample, index) => { row[sample.gsm] = formatNumber(geneValues.get(gene)?.[index]); });
    return row;
  });
  const zGeneRows = [...targetGenes].sort().map((gene) => {
    const row = { gene };
    samples.forEach((sample, index) => { row[sample.gsm] = formatNumber(zValues.get(gene)?.[index]); });
    return row;
  });

  const pairs = new Map();
  for (let sampleIndex = 0; sampleIndex < samples.length; sampleIndex += 1) {
    const sample = samples[sampleIndex];
    if (!sample.pair_id) continue;
    if (!pairs.has(sample.pair_id)) pairs.set(sample.pair_id, { pair_id: sample.pair_id, patient_id: sample.patient_id, treatment_context: sample.treatment_context, timepoint: sample.timepoint });
    pairs.get(sample.pair_id)[sample.tissue] = {
      gsm: sample.gsm,
      title: sample.title,
      scores: scores[sampleIndex],
      usedGeneCounts: usedGeneCounts[sampleIndex],
    };
  }
  const d2MetricRows = [];
  for (const signature of signatures) {
    const labelRows = [];
    const deltas = [];
    const primary = primaryEligibility(signature, config.id, config.platform, usageMap, coverageMap, selectedGenes);
    for (const pair of pairs.values()) {
      const ls = pair.lesional?.scores[signature.id];
      const nl = pair.nonlesional?.scores[signature.id];
      if (Number.isFinite(ls)) labelRows.push({ label: 1, [signature.id]: ls });
      if (Number.isFinite(nl)) labelRows.push({ label: 0, [signature.id]: nl });
      if (Number.isFinite(ls) && Number.isFinite(nl)) deltas.push(ls - nl);
    }
    const meanDelta = deltas.length ? deltas.reduce((sum, value) => sum + value, 0) / deltas.length : null;
    const sdDelta = finiteSd(deltas);
    d2MetricRows.push({
      cohort: config.id,
      signature_id: signature.id,
      d2_pairs: pairs.size,
      scored_pairs: deltas.length,
      auc: formatNumber(aucFromLabels(labelRows, signature.id)),
      mean_paired_delta_lesional_minus_nonlesional: formatNumber(meanDelta),
      median_paired_delta: formatNumber(median(deltas)),
      positive_pairs: deltas.filter((value) => value > 0).length,
      negative_pairs: deltas.filter((value) => value < 0).length,
      paired_effect_dz: formatNumber(sdDelta ? meanDelta / sdDelta : null),
      selected_gene_count: selectedGenes.get(signature.id),
      signature_gene_count: signature.genes.length,
      primary_D2_eligible: primary.primary ? 'yes' : 'no',
      declared_cohort_usage: primary.usage,
      exclusion_reason: primary.reasons,
    });
  }

  const pairedLongRows = [];
  for (const pair of pairs.values()) {
    for (const signature of signatures) {
      pairedLongRows.push({
        cohort: config.id,
        pair_id: pair.pair_id,
        patient_id: pair.patient_id,
        treatment_context: pair.treatment_context,
        timepoint: pair.timepoint,
        signature_id: signature.id,
        lesional_gsm: pair.lesional?.gsm || '',
        nonlesional_gsm: pair.nonlesional?.gsm || '',
        lesional_score: formatNumber(pair.lesional?.scores[signature.id]),
        nonlesional_score: formatNumber(pair.nonlesional?.scores[signature.id]),
        delta_lesional_minus_nonlesional: formatNumber(numericDifference(
          pair.lesional?.scores[signature.id],
          pair.nonlesional?.scores[signature.id],
        )),
      });
    }
  }

  return {
    config,
    annotation,
    expression,
    samples,
    scores,
    selectedGenes,
    manifest,
    zParameters,
    rawGeneRows,
    zGeneRows,
    d2ReferenceSampleIndices,
    d2ReferenceGsms: d2ReferenceSampleIndices.map((index) => samples[index].gsm),
    wideRows,
    d2Pairs,
    tPostSamples,
    d2MetricRows,
    pairedLongRows,
  };
}

function treatmentChangeRows(result, signatures) {
  const sampleByGsm = new Map(result.samples.map((sample, index) => [sample.gsm, { sample, score: result.scores[index] }]));
  const rows = [];
  for (const sample of result.samples) {
    if (sample.analysis_role !== 'T_posttreatment_lesional') continue;
    const baseline = sampleByGsm.get(sample.t_baseline_sample_id);
    if (!baseline) continue;
    for (const signature of signatures) {
      const baselineScore = baseline.score[signature.id];
      const postScore = result.scores[result.samples.indexOf(sample)][signature.id];
      rows.push({
        cohort: result.config.id,
        patient_id: sample.patient_id,
        timepoint: sample.timepoint,
        baseline_gsm: baseline.sample.gsm,
        posttreatment_gsm: sample.gsm,
        signature_id: signature.id,
        baseline_score: formatNumber(baselineScore),
        posttreatment_score: formatNumber(postScore),
        delta_post_minus_baseline: formatNumber(numericDifference(postScore, baselineScore)),
      });
    }
  }
  return rows;
}

async function main() {
  const registry = readCsv('data/benchmark/signature_registry.csv');
  const signatures = registry.map((row) => ({
    id: row.signature_id,
    genes: row.genes.split(';').filter(Boolean),
    compatible_endpoints: row.compatible_endpoints,
  }));
  const usageMap = new Map(readCsv('data/benchmark/signature_cohort_usage.csv')
    .map((row) => [`${row.cohort}|${row.signature_id}`, row.declared_in_registry]));
  const coverageMap = new Map(readCsv('data/benchmark/platform_gene_coverage.csv')
    .map((row) => [`${row.platform}|${row.signature_id}`, row]));
  fs.mkdirSync(OUT, { recursive: true });

  const results = [];
  for (const config of COHORTS) results.push(await processCohort(config, signatures, usageMap, coverageMap));

  const allManifestRows = [];
  const allRedundancyLongRows = [];
  const allDimensionRows = [];
  const summaryCohorts = {};
  for (const result of results) {
    const cohortOut = path.join(OUT, result.config.id);
    fs.mkdirSync(cohortOut, { recursive: true });
    const scoreColumns = [
      'cohort', 'platform', 'gsm', 'title', 'patient_id', 'tissue', 'treatment_context', 'timepoint',
      'source_block', 'analysis_role', 'pair_id', 't_pair_id', 't_baseline_sample_id',
      ...signatures.flatMap((signature) => [signature.id, `n_genes_${signature.id}`]),
    ];
    writeCsv(path.join(cohortOut, 'sample_metadata_and_scores.csv'), result.wideRows, scoreColumns);
    writeCsv(path.join(cohortOut, 'signature_gene_expression_raw.csv'), result.rawGeneRows, ['gene', ...result.samples.map((sample) => sample.gsm)]);
    writeCsv(path.join(cohortOut, 'signature_gene_expression_zscore.csv'), result.zGeneRows, ['gene', ...result.samples.map((sample) => sample.gsm)]);
    writeCsv(path.join(cohortOut, 'D2_gene_zscore_parameters.csv'), result.zParameters, [
      'gene', 'd2_reference_samples', 'd2_reference_finite_values', 'd2_reference_mean', 'd2_reference_sd',
    ]);
    writeCsv(path.join(cohortOut, 'probe_collapse_manifest.csv'), result.manifest, [
      'gene', 'annotation_probe_count', 'expression_probe_count', 'selected_probe', 'selected_mean',
      'selected_finite_fraction', 'reference_sample_count', 'selection_rule',
    ]);
    writeCsv(path.join(cohortOut, 'D2_signature_metrics.csv'), result.d2MetricRows, [
      'cohort', 'signature_id', 'd2_pairs', 'scored_pairs', 'auc', 'mean_paired_delta_lesional_minus_nonlesional',
      'median_paired_delta', 'positive_pairs', 'negative_pairs', 'paired_effect_dz', 'selected_gene_count',
      'signature_gene_count', 'primary_D2_eligible', 'declared_cohort_usage', 'exclusion_reason',
    ]);
    writeCsv(path.join(cohortOut, 'D2_paired_scores.csv'), result.pairedLongRows, [
      'cohort', 'pair_id', 'patient_id', 'treatment_context', 'timepoint', 'signature_id', 'lesional_gsm',
      'nonlesional_gsm', 'lesional_score', 'nonlesional_score', 'delta_lesional_minus_nonlesional',
    ]);
    const redundancy = redundancyOutputs(result);
    writeCsv(path.join(cohortOut, 'D2_primary_spearman_long.csv'), redundancy.longRows, [
      'cohort', 'context', 'unit', 'n', 'signature_1', 'signature_2', 'spearman_rho',
    ]);
    for (const matrixOutput of redundancy.matrixOutputs) {
      writeCsv(path.join(cohortOut, `D2_primary_spearman_${matrixOutput.context}.csv`), matrixOutput.matrixRows, [
        'cohort', 'context', 'unit', 'n', 'signature', ...matrixOutput.signatureIds,
      ]);
    }
    writeCsv(path.join(cohortOut, 'D2_primary_effective_dimension.csv'), redundancy.dimensionRows, [
      'cohort', 'context', 'unit', 'n', 'signatures', 'first_eigenvalue', 'first_eigenvalue_fraction',
      'participation_effective_dimension', 'participation_fraction_of_signatures', 'median_off_diagonal_rho',
      'min_off_diagonal_rho', 'max_off_diagonal_rho', 'signatures_included',
    ]);
    allRedundancyLongRows.push(...redundancy.longRows);
    allDimensionRows.push(...redundancy.dimensionRows);
    summaryCohorts[result.config.id] = {
      platform: result.config.platform,
      samples: result.samples.length,
      d2Pairs: result.d2Pairs,
      d2ReferenceSamples: result.d2ReferenceSampleIndices.length,
      treatmentPostSamples: result.tPostSamples,
      d2PrimaryEligibleSignatures: result.d2MetricRows.filter((row) => row.primary_D2_eligible === 'yes').map((row) => row.signature_id),
      redundancyEffectiveDimensions: Object.fromEntries(redundancy.dimensionRows.map((row) => [row.context, {
        n: row.n,
        first_eigenvalue_fraction: row.first_eigenvalue_fraction,
        participation_effective_dimension: row.participation_effective_dimension,
        median_off_diagonal_rho: row.median_off_diagonal_rho,
      }])),
      annotationRows: result.annotation.rows,
      expressionRows: result.expression.dataRows,
      selectedExpressionProbes: result.expression.selectedRows,
      missingExpressionCellsInSelectedProbes: result.expression.missingCells,
    };
    allManifestRows.push(...result.manifest.map((row) => ({ cohort: result.config.id, platform: result.config.platform, ...row })));
  }

  writeCsv(path.join(OUT, 'D2_primary_spearman_long_all.csv'), allRedundancyLongRows, [
    'cohort', 'context', 'unit', 'n', 'signature_1', 'signature_2', 'spearman_rho',
  ]);
  const spearmanSummaryRows = [];
  for (const cohort of [...new Set(allRedundancyLongRows.map((row) => row.cohort))]) {
    for (const context of [...new Set(allRedundancyLongRows.filter((row) => row.cohort === cohort).map((row) => row.context))]) {
      const rows = allRedundancyLongRows.filter((row) => row.cohort === cohort && row.context === context);
      const ordered = rows.slice().sort((a, b) => Number(a.spearman_rho) - Number(b.spearman_rho));
      const minimum = ordered[0];
      const maximum = ordered[ordered.length - 1];
      const rhos = ordered.map((row) => Number(row.spearman_rho)).filter(Number.isFinite);
      spearmanSummaryRows.push({
        cohort,
        context,
        unit: rows[0].unit,
        n: rows[0].n,
        signature_pairs: rows.length,
        median_rho: formatNumber(median(rhos)),
        minimum_rho: formatNumber(Number(minimum.spearman_rho)),
        minimum_pair: `${minimum.signature_1}-${minimum.signature_2}`,
        maximum_rho: formatNumber(Number(maximum.spearman_rho)),
        maximum_pair: `${maximum.signature_1}-${maximum.signature_2}`,
        pairs_rho_at_least_0_75: rhos.filter((rho) => rho >= 0.75).length,
        positive_pairs: rhos.filter((rho) => rho > 0).length,
        negative_pairs: rhos.filter((rho) => rho < 0).length,
      });
    }
  }
  writeCsv(path.join(OUT, 'D2_primary_spearman_pair_summary.csv'), spearmanSummaryRows, [
    'cohort', 'context', 'unit', 'n', 'signature_pairs', 'median_rho', 'minimum_rho', 'minimum_pair',
    'maximum_rho', 'maximum_pair', 'pairs_rho_at_least_0_75', 'positive_pairs', 'negative_pairs',
  ]);
  writeCsv(path.join(OUT, 'D2_primary_effective_dimension_all.csv'), allDimensionRows, [
    'cohort', 'context', 'unit', 'n', 'signatures', 'first_eigenvalue', 'first_eigenvalue_fraction',
    'participation_effective_dimension', 'participation_fraction_of_signatures', 'median_off_diagonal_rho',
    'min_off_diagonal_rho', 'max_off_diagonal_rho', 'signatures_included',
  ]);

  const treatmentResult = results.find((result) => result.config.id === 'GSE41664');
  const tChanges = treatmentChangeRows(treatmentResult, signatures);
  writeCsv(path.join(OUT, 'GSE41664_T_score_changes.csv'), tChanges, [
    'cohort', 'patient_id', 'timepoint', 'baseline_gsm', 'posttreatment_gsm', 'signature_id', 'baseline_score',
    'posttreatment_score', 'delta_post_minus_baseline',
  ]);
  const tSummary = [];
  for (const signature of signatures) {
    for (const timepoint of [...new Set(tChanges.map((row) => row.timepoint))].sort((a, b) => Number(a.replace(/\D+/g, '')) - Number(b.replace(/\D+/g, '')))) {
      const values = tChanges
        .filter((row) => row.signature_id === signature.id && row.timepoint === timepoint)
        .map((row) => Number(row.delta_post_minus_baseline))
        .filter(Number.isFinite);
      tSummary.push({
        cohort: 'GSE41664',
        signature_id: signature.id,
        timepoint,
        patients: values.length,
        mean_delta_post_minus_baseline: formatNumber(values.length ? values.reduce((sum, value) => sum + value, 0) / values.length : null),
        median_delta_post_minus_baseline: formatNumber(median(values)),
        negative_deltas: values.filter((value) => value < 0).length,
        positive_deltas: values.filter((value) => value > 0).length,
        interpretation: 'exploratory_D2_score_reversibility_only',
      });
    }
  }
  writeCsv(path.join(OUT, 'GSE41664_T_score_change_summary.csv'), tSummary, [
    'cohort', 'signature_id', 'timepoint', 'patients', 'mean_delta_post_minus_baseline',
    'median_delta_post_minus_baseline', 'negative_deltas', 'positive_deltas', 'interpretation',
  ]);
  writeCsv(path.join(OUT, 'probe_collapse_manifest_all.csv'), allManifestRows, [
    'cohort', 'platform', 'gene', 'annotation_probe_count', 'expression_probe_count', 'selected_probe',
    'selected_mean', 'selected_finite_fraction', 'reference_sample_count', 'selection_rule',
  ]);
  const availability = [
    { cohort: 'GSE41664', endpoint: 'D1', evaluable: 'no', reason: 'No healthy-skin samples in this cohort.' },
    { cohort: 'GSE41664', endpoint: 'D2', evaluable: 'yes', reason: '53 pretreatment lesional/nonlesional pairs.' },
    { cohort: 'GSE41664', endpoint: 'S', evaluable: 'no', reason: 'No validated severity endpoint in this scoring stage.' },
    { cohort: 'GSE41664', endpoint: 'T', evaluable: 'exploratory_only', reason: '51 posttreatment lesional samples; only D2 score reversibility is examined, not a validated treatment-response model.' },
    { cohort: 'GSE6710', endpoint: 'D1', evaluable: 'no', reason: 'No healthy-skin samples in this cohort.' },
    { cohort: 'GSE6710', endpoint: 'D2', evaluable: 'yes', reason: '13 paired lesional/uninvolved samples.' },
    { cohort: 'GSE6710', endpoint: 'S', evaluable: 'no', reason: 'No validated severity endpoint in this cohort.' },
    { cohort: 'GSE6710', endpoint: 'T', evaluable: 'no', reason: 'No longitudinal treatment samples.' },
  ];
  writeCsv(path.join(OUT, 'endpoint_task_availability.csv'), availability, ['cohort', 'endpoint', 'evaluable', 'reason']);

  const summary = {
    generatedAt: new Date().toISOString(),
    status: 'SCORE_STAGE_COMPLETED',
    scope: 'D2 score construction and descriptive separation in GSE41664/GSE6710; no manuscript or submission files changed.',
    rules: {
      expressionInput: 'Values from GEO series matrix as published; no additional log transform or cross-cohort normalization.',
      probeCollapse: 'For each gene, select the highest-mean probe among probes with >=95% finite values in D2 reference samples only; otherwise use a highest-mean fallback probe.',
      standardization: 'Gene-wise mean/SD is fitted on D2 reference samples only and then frozen for all samples in that cohort. GSE41664 D2 parameters exclude the 51 posttreatment samples.',
      signatureScore: 'Equal-weight mean of available gene z-scores; positive direction applied to every listed gene.',
      directionLimit: 'Published weights, intercepts, and gene-specific directions were not recovered; this is not exact model reconstruction.',
    },
    parameterBoundary: {
      GSE41664: 'D2 parameters fitted on 106 pretreatment paired samples (53 pairs); treatment-post samples use these frozen parameters.',
      GSE6710: 'D2 parameters fitted on 26 paired samples (13 pairs).',
    },
    cohorts: summaryCohorts,
    treatmentStage: {
      posttreatmentSamples: treatmentResult.tPostSamples,
      scoreDeltaRows: tChanges.length,
      claim: 'Exploratory score reversibility only; not a validated treatment-response signature.',
    },
    outputDirectory: path.relative(ROOT, OUT).replaceAll('\\', '/'),
  };
  fs.writeFileSync(path.join(OUT, 'score_run_summary.json'), `${JSON.stringify(summary, null, 2)}\n`, 'utf8');
  console.log(JSON.stringify(summary, null, 2));
}

main().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});
