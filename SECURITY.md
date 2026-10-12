# Security policy

## Supported code

Security fixes are applied to the current R implementation on `r_version`.
The Python implementation is maintained separately on `py_version`.

## Reporting a vulnerability

Please report suspected vulnerabilities privately through this repository's
GitHub security-advisory interface. Include the affected revision, a minimal
reproduction, and the expected impact. Do not publish sensitive paths, data, or
credentials in a public issue.

## Deployment boundary

The app is designed to read a prepared RDS selected by the server operator. RDS
deserialization can instantiate executable R objects, so only load artifacts
from a trusted analysis pipeline. Do not offer RDS uploads to untrusted users.

`PROXIOME_PRECOMPUTED_RDS` selects the startup artifact. The interactive
server-path field is disabled unless `PROXIOME_ALLOW_LOCAL_LOAD=true`; that flag
must remain disabled on a public deployment because it lets a visitor request
any RDS path readable by the Shiny process.

Bind Shiny to `127.0.0.1` and place it behind the deployment environment's
authenticated TLS proxy. The precomputed app does not need runtime access to
the source PXL files. Prepared sidecars must remain inside the bundle directory;
runtime readers validate paths, dataset identities, source cell mappings, and
asset checksums. Treat the bundle, sidecars, and analysis code as administrator-
managed artifacts, not user uploads.

Keep graph Parquet files and their trusted RDS index read-only in deployment.
The first access validates a sample file's checksum; later cell queries can
reuse that validation for the unchanged file. Checksums are integrity checks,
not authentication. Prepare replacements outside the live dataset directory
and deploy matching files and indexes together; never edit an index checksum
or schema version to bypass validation.

New analyses run in supervised R subprocesses with fixed selection/row/time
limits, one worker per Shiny R process, cancellation, and disconnect cleanup.
Visitors cannot submit code, SQL, source paths, or analysis functions. Deploy
behind authentication and enforce OS-level CPU/RAM limits; application row caps
are not a sandbox or a strict byte-memory limit. Multi-process hosting multiplies
the worker limit and needs an appropriate deployment-level quota.
