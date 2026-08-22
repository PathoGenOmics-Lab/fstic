# Security

Fstic is a command-line tool that reads local files and writes a matrix. It
opens no network connections and runs no code from its input. The realistic
concerns are a malformed input file causing a crash or unbounded memory use,
and a vulnerability in a dependency.

## Reporting

Please report anything security-relevant privately, through
[GitHub's advisory form](https://github.com/PathoGenOmics-Lab/fstic/security/advisories/new),
rather than as a public issue.

Include the input that triggers it, or a description sufficient to construct
one. Expect an initial reply within a week.

## Supported versions

The latest release. Fixes go into a new release rather than being backported.

## What is already covered

Pull requests are checked by CodeQL, and by a dependency review that refuses a
new dependency carrying a known advisory. `cargo audit` runs on every pull
request against the whole tree, advisory rather than blocking. Dependabot opens
security updates as soon as an advisory lands, independently of the monthly
grouped updates.

A crash on malformed input is a bug worth reporting, but it is an ordinary
issue rather than an advisory unless you can show it does more than crash.
