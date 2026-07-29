# Issue tracker: GitHub

Issues and PRDs for this repository live as GitHub issues. Use the `gh` CLI for
all operations and infer the repository from `git remote -v`.

## Conventions

- Create: `gh issue create --title "..." --body "..."`
- Read: `gh issue view <number> --comments`
- List: `gh issue list --state open --json number,title,body,labels,comments`
- Comment: `gh issue comment <number> --body "..."`
- Add or remove labels: `gh issue edit <number> --add-label "..."` or
  `gh issue edit <number> --remove-label "..."`
- Close: `gh issue close <number> --comment "..."`
- Create a child: `gh issue create --parent <number> ...`
- Add a dependency: `gh issue edit <child> --add-blocked-by <blocker>`

Use structured `--json` output for inspection and re-read every issue after a
mutation. GitHub issue and pull-request numbers share one number space; resolve
ambiguous references with `gh pr view <number>` and then `gh issue view
<number>`.

## Pull requests as a triage surface

**PRs as a request surface: no.**

## When a skill says "publish to the issue tracker"

Create a GitHub issue.

## When a skill says "fetch the relevant ticket"

Run `gh issue view <number> --comments`.

## Wayfinding operations

Use one GitHub issue as the map and native sub-issues as its tickets. Use native
GitHub issue dependencies for blocking relationships. A ticket is available
when it is open, unassigned, and has no open blocker. Claim it by assigning the
driving developer, and resolve it through a recorded answer and verified issue
closure.
