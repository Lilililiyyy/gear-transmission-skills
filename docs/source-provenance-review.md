# Source provenance review / 源码来源审查

> Review date: 2026-10-31  
> Scope: public repository contents, README and development notes, package metadata, selected source/configuration files, and the visible Git commit history.  
> Status: preliminary repository-level review, not a legal opinion or an exhaustive global code-similarity scan.

## Executive summary / 摘要

| Repository | Finding | Confidence |
|---|---|---|
| [Gearbox-Design-Studio-2D](https://github.com/Lilililiyyy/Gearbox-Design-Studio-2D) | No explicit third-party source attribution was found in the reviewed files. The Vite config retains Gemini environment-variable mappings; this is a template/configuration clue, not proof that the application code came from Gemini or another source. The repository has no LICENSE file. | Uncertain; further provenance evidence needed |
| [Cylindrical-Gear-LTCA](https://github.com/Lilililiyyy/Cylindrical-Gear-LTCA) | The repository presents itself as an independently maintained calculation kernel. MASTA is described as a numerical benchmark, not a source-code dependency. It includes an MIT license. The reviewed files and history do not establish a third-party source-code origin. | No clear external-source signal found in the reviewed material; originality is not proven |
| [gear-transmission-skills](https://github.com/Lilililiyyy/gear-transmission-skills) | The integration explicitly documents reuse of a kernel snapshot from the separate Cylindrical-Gear-LTCA project. Source maps, file hashes, and a source base commit are recorded. The skills/workflow layer is described separately from the reused calculation kernel. | Explicit internal provenance record |

## 1. Gearbox-Design-Studio-2D

### Observations

- The Vite configuration maps `GEMINI_API_KEY` to `process.env.API_KEY` and `process.env.GEMINI_API_KEY`. The reviewed application components do not show a corresponding Gemini API call. This may be leftover scaffolding or a starter-template remnant; it does not identify a specific upstream template.
- The README's old clone URL used the former `Lyy0617` account name. The visible commit history associates that handle with the same GitHub account as the current repository owner; this is not evidence of a third-party author.
- The visible history includes an initial project commit and subsequent documentation changes. Commit authorship is useful context, but does not prove who originally authored every line of code.
- No LICENSE file is present. Public visibility alone does not grant permission to reuse, modify, or redistribute the code.

### Assessment

**Uncertain / 有苗头但不能确认来源。** The Gemini environment mapping is worth tracing against the project's earliest source/template, but there is insufficient evidence to attribute the application itself to Gemini, Google AI Studio, or a named third party.

### Recommended follow-up

- Check the original local project folder, browser-based AI Studio history, exported ZIP, and any starter-template URL used at project creation.
- Compare the initial commit with any known template/export archive available to the owner.
- Decide and add an explicit LICENSE only after confirming that all included code and assets can be licensed by the repository owner.

## 2. Cylindrical-Gear-LTCA

### Observations

- The README describes MASTA as a numerical benchmark. Benchmarking against a commercial engineering application does not by itself indicate that its source code was copied.
- The repository includes an MIT LICENSE and a Python package manifest declaring MIT.
- The visible history shows initialization and subsequent project releases under the repository owner's GitHub account.
- The reviewed core files contain domain-specific Chinese comments and explicit model limitations. These characteristics are consistent with a maintained engineering implementation, but are not proof of authorship or originality.

### Assessment

**No clear third-party source-code signal was found in the reviewed material.** This should not be read as a guarantee that every algorithm, equation, benchmark input, figure, or dataset is original or independently authored. Engineering formulas and standards may be independently implemented from published references; those references should be cited in technical documentation where applicable.

### Recommended follow-up

- Keep references for standards, papers, textbooks, benchmark software, and externally supplied datasets close to the relevant model or case documentation.
- Preserve the distinction between numerical benchmark results and source-code reuse.
- Confirm that the MIT license is appropriate for all included files and figures.

## 3. gear-transmission-skills

### Observations

- `docs/development.md` explicitly says that `packages/cylindrical-gear-ltca/` contains a byte-preserved snapshot of the separate calculation project.
- `src/gear_transmission/kernel_snapshot.json` records the source project, source base commit, snapshot version, and SHA-256 values for the packaged files.
- `docs/source-map.json` maps selected functions and source files to their integrated paths and records hashes.
- The snapshot metadata states that it reflects the source working tree at the time of capture, including uncommitted updates; it is not represented as a clean export of the recorded base commit.
- The skill/workflow adapter under `src/gear_transmission/` is described as newly written integration code, while the computational kernel is maintained in the separate project.
- The repository contains an MIT LICENSE. The packaged kernel also carries its own license file in the snapshot.

### Assessment

**Source reuse is explicit and documented.** This is an internal reuse relationship between two repositories under the same account, not an unexplained third-party code appearance. Keep the source project link, license notices, source map, and hash records intact when updating the snapshot.

### Recommended follow-up

- When updating the kernel snapshot, regenerate the file hashes and source map together.
- Keep the upstream source commit and the actual snapshot state clearly distinguished, especially when the snapshot includes uncommitted upstream changes.
- Review third-party licenses for dependencies and any future imported examples or assets.

## Review limitations / 审查边界

This review inspected repository metadata, visible commit history, README/development documentation, selected configuration and core files, license files, and the explicit kernel provenance records. It did not establish a complete chain of custody for the source code, inspect the owner's private local history or AI-tool conversations, run a comprehensive code clone detector across all public/private repositories, or provide a legal determination. No source should be labelled "copied", "AI-generated", or "fully original" based only on style, comments, commit authorship, or a configuration remnant.

For a stronger conclusion, provide the earliest local project archive, source-template links, AI Studio export/history (if applicable), and references used to derive the engineering models.
