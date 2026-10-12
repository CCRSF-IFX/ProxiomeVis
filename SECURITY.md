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
the source PXL files.
