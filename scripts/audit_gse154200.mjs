import fs from 'node:fs';
import path from 'node:path';
import readline from 'node:readline';
import zlib from 'node:zlib';
import { fileURLToPath } from 'node:url';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const SERIES = path.join(ROOT, 'sources/raw/GSE154200/GSE154200_series_matrix.txt.gz');
const COUNTS = path.join(ROOT, 'sources/raw/GSE154200/GSE154200_countMat.txt.gz');
const SOFT = path.join(ROOT, 'sources/raw/GSE154200/GSE154200_samples_brief.txt');
const OUT = path.join(ROOT, 'data/benchmark');

function tabFields(line) {
  return line.split('\t').map((value) => {
    let clean = value.trim();
    if (clean.startsWith('"') && clean.endsWith('"')) clean = clean.slice(1, -1);
    return clean.replaceAll('""', '"');
  });
}

function csvEscape(value) {
  const text = value == null ? '' : String(value);
  return /[",\r\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

function writeCsv(file, rows, columns) {
  const lines = [columns.join(','), ...rows.map((row) => columns.map((column) => csvEscape(row[column])).join(','))];
  fs.mkdirSync(path.dirname(file), { recursive: true });
  fs.writeFileSync(file, `${lines.join('\n')}\n`, 'utf8');
}

function characteristicKey(value) {
  const match = value.match(/^([^:]+):\s*(.*)$/);
  if (!match) return null;
  return {
    key: match[1].trim().toLowerCase().replace(/[^a-z0-9]+/g, '_').replace(/^_+|_+$/g, ''),
    value: match[2].trim(),
  };
}

async function readSeriesMetadata() {
  const input = fs.createReadStream(SERIES).pipe(zlib.createGunzip());
  const reader = readline.createInterface({ input, crlfDelay: Infinity });
  const fields = new Map();
  for await (const rawLine of reader) {
    const match = rawLine.match(/^!Sample_([^\t]+)/);
    if (!match) continue;
    const key = match[1];
    if (!fields.has(key)) fields.set(key, []);
    fields.get(key).push(tabFields(rawLine).slice(1));
  }
  input.destroy();

  const gsmIds = fields.get('geo_accession')?.[0] || [];
  const titles = fields.get('title')?.[0] || [];
  const sourceNames = fields.get('source_name_ch1')?.[0] || [];
  const characteristicRows = fields.get('characteristics_ch1') || [];
  if (!gsmIds.length || gsmIds.length !== titles.length || gsmIds.length !== sourceNames.length) {
    throw new Error('Incomplete GEO sample metadata.');
  }

  const characteristicKeys = new Set();
  const samples = gsmIds.map((gsm, sampleIndex) => {
    const characteristics = {};
    for (const row of characteristicRows) {
      const parsed = characteristicKey((row[sampleIndex] || '').trim());
      if (!parsed) continue;
      characteristicKeys.add(parsed.key);
      characteristics[parsed.key] = parsed.value;
    }
    return {
      gsm,
      title: titles[sampleIndex] || '',
      source_name: sourceNames[sampleIndex] || '',
      characteristics,
    };
  });
  return { samples, characteristicKeys: [...characteristicKeys].sort() };
}

async function readCountTitles() {
  const input = fs.createReadStream(COUNTS).pipe(zlib.createGunzip());
  const reader = readline.createInterface({ input, crlfDelay: Infinity });
  let titles = null;
  for await (const line of reader) {
    titles = tabFields(line).slice(1);
    break;
  }
  input.destroy();
  if (!titles) throw new Error('No header in count matrix.');
  return titles;
}

async function readSoftSampleFields() {
  const reader = readline.createInterface({
    input: fs.createReadStream(SOFT, 'utf8'),
    crlfDelay: Infinity,
  });
  const samples = new Map();
  let current = null;
  for await (const rawLine of reader) {
    const line = rawLine.replace(/^\uFEFF/, '');
    if (line.startsWith('^SAMPLE = ')) {
      current = line.slice('^SAMPLE = '.length).trim();
      samples.set(current, {});
      continue;
    }
    if (!current || !line.startsWith('!Sample_')) continue;
    const separator = line.indexOf(' = ');
    if (separator < 0) continue;
    const key = line.slice('!Sample_'.length, separator);
    const value = line.slice(separator + 3).trim();
    samples.get(current)[key] = value;
  }
  return samples;
}

function normalizeId(value) {
  return String(value || '').trim().toLowerCase().replace(/[^a-z0-9]+/g, '_').replace(/^_+|_+$/g, '');
}

function tissueFromSource(sourceName) {
  const source = String(sourceName || '').toLowerCase();
  if (/non[- ]?lesional/.test(source)) return 'nonlesional';
  if (/lesional/.test(source)) return 'lesional';
  return '';
}

function treatmentFields(characteristics) {
  return Object.entries(characteristics)
    .filter(([key]) => /treatment|therapy|drug|medication/.test(key))
    .map(([key, value]) => `${key}=${value}`)
    .join(';');
}

function classifySample(sample, softSample) {
  const diagnosis = sample.characteristics.diagnosis || '';
  const psoriasis = /psoriasis/i.test(diagnosis) || /psoriasis/i.test(sample.source_name);
  const tissue = tissueFromSource(sample.source_name);
  const patientId = sample.characteristics.patient_id || sample.characteristics.patient_identifyier || '';
  return {
    ...sample,
    diagnosis,
    psoriasis,
    tissue,
    patient_id: normalizeId(patientId),
    treatment: treatmentFields(sample.characteristics) || softSample?.treatment_protocol_ch1 || '',
    timepoint: sample.characteristics.timepoint || sample.characteristics.time_point || '',
    batch_id: sample.characteristics.batch_id || sample.characteristics.batch || '',
  };
}

function main() {
  return readSeriesMetadata().then(async ({ samples, characteristicKeys }) => {
    const countTitles = await readCountTitles();
    const softSamples = await readSoftSampleFields();
    const titleCounts = new Map();
    for (const title of countTitles) titleCounts.set(title, (titleCounts.get(title) || 0) + 1);
    const duplicateCountTitles = [...titleCounts.entries()].filter(([, count]) => count > 1).map(([title]) => title);
    const classified = samples.map((sample) => classifySample(sample, softSamples.get(sample.gsm)));
    const metadataTitles = new Set(classified.map((sample) => sample.title));
    const countTitleSet = new Set(countTitles);
    const missingFromCounts = classified.filter((sample) => !countTitleSet.has(sample.title)).map((sample) => sample.gsm);
    const extraCountTitles = countTitles.filter((title) => !metadataTitles.has(title));

    const psoriasis = classified.filter((sample) => sample.psoriasis);
    const byPatient = new Map();
    for (const sample of psoriasis) {
      if (!sample.patient_id) continue;
      if (!byPatient.has(sample.patient_id)) byPatient.set(sample.patient_id, []);
      byPatient.get(sample.patient_id).push(sample);
    }

    const pairRows = [];
    const eligiblePatientIds = new Set();
    for (const [patientId, patientSamples] of [...byPatient.entries()].sort(([a], [b]) => a.localeCompare(b))) {
      const lesionals = patientSamples.filter((sample) => sample.tissue === 'lesional');
      const nonlesionals = patientSamples.filter((sample) => sample.tissue === 'nonlesional');
      const treatments = [...new Set(patientSamples.map((sample) => sample.treatment).filter(Boolean))];
      const timepoints = [...new Set(patientSamples.map((sample) => sample.timepoint).filter(Boolean))];
      const pairEligible = patientSamples.length === 2 && lesionals.length === 1 && nonlesionals.length === 1;
      if (pairEligible) eligiblePatientIds.add(patientId);
      pairRows.push({
        patient_id: patientId,
        psoriasis_samples: patientSamples.length,
        lesional_samples: lesionals.length,
        nonlesional_samples: nonlesionals.length,
        other_samples: patientSamples.length - lesionals.length - nonlesionals.length,
        treatment_values: treatments.join(';'),
        timepoint_values: timepoints.join(';'),
        pair_eligible: pairEligible ? 'yes' : 'no',
        gsm_lesional: lesionals.map((sample) => sample.gsm).join(';'),
        gsm_nonlesional: nonlesionals.map((sample) => sample.gsm).join(';'),
        titles: patientSamples.map((sample) => sample.title).join(';'),
      });
    }

    const sampleRows = classified.map((sample) => {
      const pairEligible = eligiblePatientIds.has(sample.patient_id);
      let inclusionStatus = 'excluded';
      let exclusionReason = '';
      if (!sample.psoriasis) exclusionReason = `non_psoriasis:${sample.diagnosis || 'unknown'}`;
      else if (!sample.patient_id) exclusionReason = 'missing_patient_id';
      else if (!sample.tissue) exclusionReason = `unrecognized_source:${sample.source_name}`;
      else if (!pairEligible) exclusionReason = 'psoriasis_patient_without_complete_lesional_nonlesional_pair';
      else {
        inclusionStatus = 'included_D2';
        exclusionReason = '';
      }
      return {
        gsm: sample.gsm,
        title: sample.title,
        source_name: sample.source_name,
        diagnosis: sample.diagnosis,
        patient_id: sample.patient_id,
        tissue: sample.tissue,
        treatment: sample.treatment,
        timepoint: sample.timepoint,
        batch_id: sample.batch_id,
        title_in_count_matrix: countTitleSet.has(sample.title) ? 'yes' : 'no',
        inclusion_status: inclusionStatus,
        exclusion_reason: exclusionReason,
      };
    });

    writeCsv(path.join(OUT, 'GSE154200_cohort_screening.csv'), sampleRows, [
      'gsm', 'title', 'source_name', 'diagnosis', 'patient_id', 'tissue', 'treatment', 'timepoint',
      'batch_id', 'title_in_count_matrix', 'inclusion_status', 'exclusion_reason',
    ]);
    writeCsv(path.join(OUT, 'GSE154200_pair_audit.csv'), pairRows, [
      'patient_id', 'psoriasis_samples', 'lesional_samples', 'nonlesional_samples', 'other_samples',
      'treatment_values', 'timepoint_values', 'pair_eligible', 'gsm_lesional', 'gsm_nonlesional', 'titles',
    ]);

    const included = sampleRows.filter((row) => row.inclusion_status === 'included_D2');
    const summary = {
      generatedAt: new Date().toISOString(),
      accession: 'GSE154200',
      platform: 'GPL20301',
      metadataSamples: classified.length,
      countMatrixSamples: countTitles.length,
      countTitleUnique: duplicateCountTitles.length === 0,
      duplicateCountTitles,
      metadataTitlesMissingFromCounts: missingFromCounts,
      extraCountTitles,
      characteristicKeys,
      psoriasisSamples: psoriasis.length,
      psoriasisPatients: byPatient.size,
      completeD2Pairs: included.length / 2,
      includedD2Samples: included.length,
      includedD2Lesional: included.filter((row) => row.tissue === 'lesional').length,
      includedD2Nonlesional: included.filter((row) => row.tissue === 'nonlesional').length,
      treatmentValueCounts: [...new Set(psoriasis.map((sample) => sample.treatment))].map((value) => ({
        value,
        samples: psoriasis.filter((sample) => sample.treatment === value).length,
      })),
      includedTreatmentValueCounts: [...new Set(included.map((row) => row.treatment))].map((value) => ({
        value,
        samples: included.filter((row) => row.treatment === value).length,
      })),
      timepointValueCounts: [...new Set(psoriasis.map((sample) => sample.timepoint))].map((value) => ({
        value,
        samples: psoriasis.filter((sample) => sample.timepoint === value).length,
      })),
      predeclaredRules: {
        disease: 'diagnosis equals psoriasis or source_name_ch1 contains psoriasis',
        pair: 'exactly two psoriasis samples from one patient: one lesional and one non-lesional',
        exclusions: 'all non-psoriasis samples and all psoriasis patients without a complete pair',
        treatmentAudit: 'all characteristic fields whose key contains treatment, therapy, drug, or medication are reported; no post hoc exclusion by outcome',
        treatmentProtocol: 'when series-matrix characteristics contain no treatment field, !Sample_treatment_protocol_ch1 is read from the immutable GEO SOFT sample snapshot',
      },
      outputs: {
        sampleAudit: 'data/benchmark/GSE154200_cohort_screening.csv',
        pairAudit: 'data/benchmark/GSE154200_pair_audit.csv',
      },
    };
    fs.writeFileSync(path.join(OUT, 'GSE154200_audit_summary.json'), `${JSON.stringify(summary, null, 2)}\n`, 'utf8');
    console.log(JSON.stringify(summary, null, 2));
  });
}

main().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});
