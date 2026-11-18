# Data manifest

The recordings are not part of this repository. They are public under CC0:

- OpenNeuro ds008302, version 1.0.0: doi:10.18112/openneuro.ds008302.v1.0.0
- NEMAR mirror on008302: doi:10.82901/nemar.on008302
- Metadata pinned to GitHub commit af3998a993adb108dc29a58b70080bcc72a42098 of nemarDatasets/on008302

| File | Contents |
| --- | --- |
| `dataset_manifest.tsv` | Every file of version 1.0.0 (4,506 files, 216.9 GB): path, participant, task, source, size, SHA256, S3 version id |
| `download_report_cohort.tsv` | The 742 files of the cohort selection used in the analysis (memory task, 185.9 GB) with size, SHA256 and verification status |

To obtain the files: set the data location (see the main README) and run `scripts\run_analysis.cmd 01_download_metadata`, then `02b_download_cohort`. The downloader checks every file against its SHA256 and resumes interrupted downloads. Nothing is downloaded unless you run these steps.
