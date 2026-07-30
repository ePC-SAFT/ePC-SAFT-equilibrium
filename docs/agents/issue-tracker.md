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

<!-- project-truss:setup:start -->
## Project Truss tracker contract

Project Truss uses GitHub Issues in `ePC-SAFT/ePC-SAFT-equilibrium` through the authenticated `gh` CLI.

## Project Truss execution issues

- Root issues use Problem Statement, Solution, numbered User Stories, Implementation Decisions, Testing Decisions, Out of Scope, and Further Notes.
- Executable leaves use Parent, What to build, Acceptance criteria, and Blocked by.
- Native sub-issue and blocked-by relationships are authoritative.

## Wayfinder decision issues

Wayfinder maps and decision tickets clarify an uncertain destination before execution shaping. They are not Project Truss execution contracts. Their labels, assignees, comments, and closure never derive Truss readiness or completion. Link a completed map from a new Truss root; do not automatically convert or reuse its issues.

### Wayfinding operations

- Create one map issue with the `wayfinder:map` label and `Destination`, `Notes`, `Decisions so far`, `Not yet specified`, and `Out of scope` sections.
- Create decision tickets as native child issues with a `## Question` body and one `wayfinder:research`, `wayfinder:prototype`, `wayfinder:grilling`, or `wayfinder:task` label.
- Express blocking with native GitHub issue dependencies and derive the frontier from open, unblocked, unassigned children.
- Claim by assignment. Resolve with one answer comment, close the ticket, and append only a linked gist to the map's `Decisions so far`.

## Authority

GitHub, Git, CI, reviews, and current worktrees are lifecycle truth. Labels, milestones, Projects, issue text, and comments cannot grant authority.
<!-- project-truss:setup:end -->
