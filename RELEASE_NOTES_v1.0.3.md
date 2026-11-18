# v1.0.3

Archived software release for the PeerJ Research Article:

> Limited cross-cohort added value of published psoriasis transcriptomic signatures: a three-cohort reanalysis

Analysis code and statistical outputs are unchanged from v1.0.0. This release hardens the reproducibility tooling:

- repository checks resolve paths from their own location instead of the caller's working directory;
- the SHA256 manifest is based on tracked or publishable files and excludes local-only `ZENODO_DOI.txt`;
- verification compares files against the committed manifest by default and requires an explicit `-UpdateManifest` flag to rewrite it;
- the one-command workflow checks for Node.js and resolves public inputs from the repository root.

The code remains a transparent comparator benchmark, not an exact reconstruction of the original published models.
