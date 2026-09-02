# Judgment Protocol

Use this file to encode personal research judgment. Update it after real reading, experiments, and advisor feedback.

Paper depth gate:
- Skim if the paper is only background, trend mapping, or loosely related.
- Read carefully if it defines a baseline, introduces a method close to the current project, or changes the feasibility of an active idea.
- Deep read if its assumptions, equations, or experiments directly affect the current thesis claim.

Idea novelty gate:
- A useful idea should change at least one of: problem formulation, assumption set, algorithmic mechanism, complexity tradeoff, evaluation scenario, or implementation path.
- Naming, parameter tuning, or combining fashionable topics is not enough unless it creates a measurable new capability.

Feasibility gate:
- Prefer ideas that can produce an interpretable result within the current codebase or a small extension of it.
- Penalize ideas requiring unavailable datasets, unclear baselines, excessive compute, or major theory development before any experiment is possible.
- Reject or narrow ideas whose success depends on unrealistic CSI, unfair baselines, or unverifiable claims.

Result support gate:
- A result supports a claim only if the metric, baseline, and scenario match the claim.
- A single best-case curve is not enough for a broad robustness claim.
- A method should report the cost paid for performance gains: complexity, overhead, PAPR, BER, or sensitivity.

Advisor feedback loop:
- Convert feedback into one of: revise claim, add baseline, fix assumption, improve explanation, run experiment, narrow scope, or discard direction.
- Keep the next action concrete enough to execute within one work session.

