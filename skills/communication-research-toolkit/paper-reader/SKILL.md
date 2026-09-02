---
name: comm-paper-reader
description: Read and explain communication research papers, especially waveform, AFDM, channel estimation, detection, PAPR, BER, and receiver-design work. Use when the user wants to understand a paper, formulas, figures, experiments, metadata, code links, or prerequisites. Focus on understanding, not full peer review.
---

Read before execution:

- `../references/research-principles.md`
- `../references/evidence-policy.md`
- `../references/search-policy.md`
- `../references/failure-handling.md`
- `../references/communication-domain.md`

# Communication Paper Reader

Help the user understand a communication paper efficiently and accurately.

## Scope

Responsible for:

- metadata extraction
- contribution summary
- system model explanation
- method and equation explanation
- experiment interpretation
- figure and table explanation
- learning roadmap

Not responsible for:

- full accept/reject style review
- implementation planning
- generating new research proposals

## Workflow

1. Extract metadata: title, authors, affiliation, venue, year, DOI, paper link, code repository, project page.
2. Summarize the core contribution in a few precise sentences.
3. Explain the problem background, motivation, and prior limitations.
4. Identify the communication system model: waveform, channel, mobility, CSI, receiver assumptions, metrics, constraints.
5. Explain the method: inputs, outputs, pipeline, key modules, equations, and intuition.
6. Interpret experiments: datasets or simulations, baselines, metrics, parameter settings, and what each result actually supports.
7. Explain figures and tables one by one when available.
8. Build a concept hierarchy or lightweight knowledge graph.
9. Provide prerequisites and a learning roadmap.
10. Give a brief assessment of novelty, limitations, and possible application, without turning it into a full review.

## Output Style

Use clear structure and plain language. Explain intuition before equations. Mark missing metadata with `[NOT FOUND]` and unsupported claims with `[INSUFFICIENT EVIDENCE]`.

