---
name: comm-research-idea-generator
description: Generate communication research ideas from papers, experiment results, advisor feedback, or known gaps. Use when the user asks for research gaps, thesis ideas, innovation points, next directions, or how to turn a limitation into a feasible communication research project.
---

Read before execution:

- `../references/research-principles.md`
- `../references/evidence-policy.md`
- `../references/failure-handling.md`
- `../references/communication-domain.md`
- `../references/judgment-protocol.md`

# Communication Research Idea Generator

Generate research ideas from evidence and constraints, not hype.

## Scope

Responsible for:

- gap extraction
- idea generation
- novelty framing
- feasibility screening
- experiment direction suggestions

Not responsible for:

- full paper review
- detailed implementation plan
- literature survey without a seed topic

## Workflow

1. Identify the research context: current project, paper cluster, system model, and target claim.
2. Extract explicit and implicit gaps: assumptions, missing scenarios, missing baselines, weak metrics, scalability limits, or deployment barriers.
3. Generate candidate ideas in three tiers: quick test, thesis-scale direction, and higher-risk extension.
4. For each idea, state the novelty source: formulation, mechanism, assumption relaxation, complexity tradeoff, evaluation setting, or implementation route.
5. Apply feasibility checks: available codebase, baseline clarity, experiment cost, required theory, and risk of unverifiable claims.
6. Suggest the smallest experiment that can validate or kill the idea.

## Output Style

Prefer a small number of well-constrained ideas over many vague ideas. For weak ideas, say why they are weak and how to narrow them.

