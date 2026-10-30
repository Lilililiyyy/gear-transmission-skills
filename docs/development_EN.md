# Development, kernel integration, and validation

[简体中文](development.md)

This library maintains agent workflows and access to computational capabilities. Independent projects maintain models, algorithms, complete APIs, and examples. The currently integrated project is [Cylindrical Gear LTCA](https://github.com/Lyy0617/Cylindrical-Gear-LTCA).

## Code ownership

| Location | Responsibility |
|---|---|
| `skills/<skill-name>/` | Task triggers, workflows, resources, and deliverables |
| `src/gear_transmission/` | Explicit unit conversion, input validation, function invocation, reports, and independent reruns |
| `packages/cylindrical-gear-ltca/` | Byte-preserved snapshot of the independent reconstructed computational kernel |
| `tests/` | Automated checks for input adaptation, raw arrays, failure states, and recomputation |
| `evaluations/` | Scenario lists, actual acceptance records, and independent agent-use records |

New code documentation and diagnostics use Chinese; user-facing entry documentation is maintained in Chinese and English. Original functions, protocol fields, units, and indices are retained. Skill scripts and workflows call the kernel without reimplementing its formulas.

## Adding a skill

Define the user task, required inputs, computational entry point, and acceptance criteria before creating `skills/<skill-name>/SKILL.md`. Use lowercase letters, digits, and hyphens in names, with a description that states when the skill applies.

Add resources as required:

```text
skills/<skill-name>/
├─ SKILL.md          # Required: name, trigger description, and essential workflow
├─ agents/           # Optional: agent presentation and invocation metadata
├─ scripts/          # Optional: repeatable executable helpers
├─ references/       # Optional: inputs, outputs, units, and conditional procedures
└─ assets/           # Optional: example inputs and output templates
```

Keep the entry point concise and move conditional details into references. State dependencies, model selection, missing-input handling, and failure behavior. Installing skill files must not be described as installing all computational dependencies.

Validate an explicit request first, then extend evaluation to implicit activation, clarification of missing information, and out-of-scope tasks. State the actual status in the README catalog; empty directories or placeholder skills do not establish an integration.

## Updating a computational kernel

Snapshot identity is stored in [kernel_snapshot.json](../src/gear_transmission/kernel_snapshot.json). Function, file-hash, and line mappings are recorded in [source-map.json](source-map.json). When the source project contains uncommitted changes, its commit identifier is background context; snapshot hashes identify the actual files.

Make algorithm changes in the independent project first, then explicitly update the integrated snapshot, mappings, and acceptance records. Preserve formulas, data layouts, units, signs, and branches. Continued upstream development does not automatically replace an accepted snapshot, and formulas should not be maintained in two locations simultaneously.

See the [kernel snapshot README](../packages/cylindrical-gear-ltca/README.md) for integration details. Consult the independent project for models and complete interfaces. The bundled [model documentation](../packages/cylindrical-gear-ltca/docs/model.md) and [API documentation](../packages/cylindrical-gear-ltca/docs/api.md), currently in Chinese, correspond to the integrated version.

## Local checks

Run in the installed environment:

```bash
python -m unittest discover -s tests -v
```

Existing checks cover unit conversion, duplicate inputs, material and geometry types, internal helical meshing, ISO settings, zero load, output protection, solver failures, raw arrays, and independent reruns. Add behavior-focused checks when changing the input protocol or computational calls; ordinary documentation changes require link and command checks.

Execute new or modified scripts. Complex skills should also be evaluated by an independent agent using realistic requests. [scenarios.json](../evaluations/scenarios.json) lists evaluation scenarios; a listed scenario is not evidence that it has been executed.

## Initial acceptance

Actual records are in [acceptance.json](../evaluations/acceptance.json) and [forward-test.json](../evaluations/forward-test.json):

| Check | Recorded result |
|---|---|
| Computational workflow checks | 14 passed |
| Installation checks | 9 passed; complete resources, updates, read-only verification, edit protection, and Git fast-forward |
| Selected spur benchmark against the original project | 22 raw arrays identical |
| Independent single-condition rerun | 22 raw arrays passed |
| Independent rerun of 3 candidates × 3 conditions | 198 raw arrays passed |
| Independent agent using a different helical input | Actual invocation, plotting, and 22-array rerun passed |
| Execution and recomputation after wheel installation | Passed |
| Remote CI | Four Windows / Linux × Python 3.10 / 3.11 configurations passed; see the [CI record](../evaluations/ci.json) |

Record execution completion, solver termination criteria, numerical reproducibility, spatial / phase convergence, original-component response accuracy, and physical testing separately. Unperformed validation remains unverified. Check whether independent-project model changes or validation conclusions apply to the integrated snapshot and operating conditions.

## Structural references

The library adopts the parts of these projects that suit its current scope:

- [Anthropic Skills](https://github.com/anthropics/skills): skills grouped under `skills/`, each with a clear `SKILL.md` entry point.
- [Vercel Agent Skills](https://github.com/vercel-labs/agent-skills): available skills, installation, and invocation examples near the front of a concise README.
- [Scientific Agent Skills](https://github.com/K-Dense-AI/scientific-agent-skills): separate skill registration from scientific dependencies, load references as needed, and retain workflow validation.

The skill format follows the [Agent Skills specification](https://agentskills.io/specification).
