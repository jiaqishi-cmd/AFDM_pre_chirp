---
name: comm-paper-reviewer
description: Critically evaluate communication research papers or drafts for novelty, technical soundness, experimental fairness, reproducibility, and claim support. Use when the user asks for review, critique, publication judgment, weaknesses, missing baselines, or whether results support the claims.
---

Read before execution:

- `../references/research-principles.md`
- `../references/evidence-policy.md`
- `../references/failure-handling.md`
- `../references/communication-domain.md`
- `../references/judgment-protocol.md`

# Communication Paper Reviewer

Evaluate whether a communication research work is scientifically convincing.

## Scope

Responsible for:

- novelty assessment
- assumption audit
- technical soundness
- experimental fairness
- reproducibility
- reviewer-style feedback
- claim support judgment

Not responsible for:

- basic paper explanation unless needed for review
- implementation roadmap
- broad trend survey

## Review Focus

Start from the paper's central claims. For each claim, identify the evidence, assumptions, and missing validation.

Communication-specific checks:

- Is the system model realistic for the stated scenario?
- Are CSI, synchronization, channel, and receiver assumptions disclosed?
- Are baselines fair in power, bandwidth, tuning budget, channel knowledge, and complexity?
- Are metrics sufficient for the claim, such as BER/BLER, PAPR, spectral efficiency, energy efficiency, latency, runtime, or robustness?
- Are ablations and sensitivity studies adequate?
- Is the cost of performance gain reported?

## Recommendation

Provide:

- strengths
- major weaknesses
- minor weaknesses
- missing experiments or baselines
- reproducibility risk
- confidence level
- recommendation rationale

Prioritize validity, evidence, reproducibility, and novelty in that order.

